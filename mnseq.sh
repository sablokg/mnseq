#!/usr/bin/env bash
#===============================================================================
# MNase-seq Analysis Pipeline
#
# Steps:
#   1. Raw read QC              (FastQC)
#   2. Adapter/quality trimming (Trim Galore / cutadapt)
#   3. Alignment                (Bowtie2, paired-end, short-fragment tuned)
#   4. BAM cleanup              (sort, mark dups, filter MAPQ, remove dups)
#   5. Fragment size QC         (nucleosome-sized ~147 bp enrichment check)
#   6. Fragment-size filtering  (keep mono-nucleosome fraction, e.g. 100-200bp)
#   7. Coverage tracks          (bamCoverage -> bigWig, normalized)
#   8. Nucleosome calling       (DANPOS3 dpos)
#   9. Nucleosome occupancy/positioning QC summary
#
# Requirements (conda/mamba recommended):
#   fastqc, trim-galore (cutadapt), bowtie2, samtools, picard,
#   deeptools (bamCoverage, bamPEFragmentSize), danpos3 (or danpos.py), python3, R
#
#   conda create -n mnaseseq -c bioconda -c conda-forge \
#       fastqc trim-galore bowtie2 samtools picard deeptools danpos
#
# Usage:
#   ./mnase_seq_pipeline.sh -1 R1.fastq.gz -2 R2.fastq.gz -x /path/to/bowtie2_index \
#       -g genome.fa -o results -s sample1 -t 8
#===============================================================================

set -euo pipefail
IFS=$'\n\t'

#-------------------------------- Defaults ------------------------------------
THREADS=8
OUTDIR="mnase_results"
SAMPLE="sample"
GENOME_FASTA=""
BT2_INDEX=""
FASTQ_R1=""
FASTQ_R2=""
MIN_INSERT=100      # mono-nucleosome lower bound (bp)
MAX_INSERT=200      # mono-nucleosome upper bound (bp)
MAPQ_MIN=30
KEEP_TMP=false

usage() {
    cat <<EOF
Usage: $0 -1 R1.fastq.gz -2 R2.fastq.gz -x BOWTIE2_INDEX_PREFIX -g GENOME.fa -s SAMPLE_NAME [options]

Required:
  -1  FILE   Forward reads (fastq.gz)
  -2  FILE   Reverse reads (fastq.gz)
  -x  PATH   Bowtie2 index prefix (built ahead of time with bowtie2-build)
  -g  FILE   Reference genome FASTA (needed for chrom sizes / DANPOS)
  -s  NAME   Sample name (used for output prefixes)

Options:
  -o  DIR    Output directory              (default: $OUTDIR)
  -t  INT    Threads                       (default: $THREADS)
  -m  INT    Min mono-nucleosome insert bp (default: $MIN_INSERT)
  -M  INT    Max mono-nucleosome insert bp (default: $MAX_INSERT)
  -q  INT    Minimum MAPQ to keep          (default: $MAPQ_MIN)
  -k         Keep intermediate files
  -h         Show this help
EOF
    exit 1
}

while getopts "1:2:x:g:s:o:t:m:M:q:kh" opt; do
    case $opt in
        1) FASTQ_R1="$OPTARG" ;;
        2) FASTQ_R2="$OPTARG" ;;
        x) BT2_INDEX="$OPTARG" ;;
        g) GENOME_FASTA="$OPTARG" ;;
        s) SAMPLE="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MIN_INSERT="$OPTARG" ;;
        M) MAX_INSERT="$OPTARG" ;;
        q) MAPQ_MIN="$OPTARG" ;;
        k) KEEP_TMP=true ;;
        h|*) usage ;;
    esac
done

[[ -z "$FASTQ_R1" || -z "$FASTQ_R2" || -z "$BT2_INDEX" || -z "$GENOME_FASTA" ]] && usage

for exe in fastqc trim_galore bowtie2 samtools picard bamCoverage bamPEFragmentSize; do
    command -v "$exe" >/dev/null 2>&1 || { echo "ERROR: '$exe' not found in PATH." >&2; exit 1; }
done

#------------------------------- Directories ----------------------------------
QC_RAW="$OUTDIR/01_fastqc_raw"
TRIM_DIR="$OUTDIR/02_trimmed"
ALIGN_DIR="$OUTDIR/03_aligned"
FILT_DIR="$OUTDIR/04_filtered"
FRAG_DIR="$OUTDIR/05_fragment_qc"
NUCFILT_DIR="$OUTDIR/06_nucleosome_bam"
TRACK_DIR="$OUTDIR/07_bigwig"
DANPOS_DIR="$OUTDIR/08_danpos"
LOG_DIR="$OUTDIR/logs"

mkdir -p "$QC_RAW" "$TRIM_DIR" "$ALIGN_DIR" "$FILT_DIR" "$FRAG_DIR" \
         "$NUCFILT_DIR" "$TRACK_DIR" "$DANPOS_DIR" "$LOG_DIR"

log() { echo -e "\n[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_DIR/${SAMPLE}.pipeline.log"; }

#===============================================================================
# 1. Raw QC
#===============================================================================
log "STEP 1/9: FastQC on raw reads"
fastqc -t "$THREADS" -o "$QC_RAW" "$FASTQ_R1" "$FASTQ_R2" \
    &>> "$LOG_DIR/${SAMPLE}.fastqc_raw.log"

#===============================================================================
# 2. Adapter / quality trimming
#    --nextera / --illumina auto-detected by trim_galore; MNase inserts can be
#    short (<read length) so adapter read-through trimming matters a lot.
#===============================================================================
log "STEP 2/9: Trim Galore adapter/quality trimming"
trim_galore --paired --cores "$THREADS" --quality 20 --length 20 \
    --output_dir "$TRIM_DIR" \
    "$FASTQ_R1" "$FASTQ_R2" \
    &>> "$LOG_DIR/${SAMPLE}.trim_galore.log"

TRIM_R1=$(find "$TRIM_DIR" -name "*_val_1.fq.gz")
TRIM_R2=$(find "$TRIM_DIR" -name "*_val_2.fq.gz")

#===============================================================================
# 3. Alignment
#    -X 1000: allow full range of fragment sizes (MNase captures sub-, mono-,
#    di-, tri-nucleosome fragments). --no-mixed/--no-discordant enforce clean
#    proper pairs, which matters for downstream fragment-size-based filtering.
#===============================================================================
log "STEP 3/9: Bowtie2 alignment"
bowtie2 -p "$THREADS" -x "$BT2_INDEX" \
    -1 "$TRIM_R1" -2 "$TRIM_R2" \
    -X 1000 --no-mixed --no-discordant --very-sensitive \
    2> "$LOG_DIR/${SAMPLE}.bowtie2.log" \
    | samtools sort -@ "$THREADS" -o "$ALIGN_DIR/${SAMPLE}.sorted.bam" -
samtools index -@ "$THREADS" "$ALIGN_DIR/${SAMPLE}.sorted.bam"

#===============================================================================
# 4. BAM cleanup: mark/remove duplicates, filter MAPQ, keep proper pairs only,
#    drop mitochondrial reads (common MNase-seq QC step; adjust chrom name).
#===============================================================================
log "STEP 4/9: Mark duplicates (Picard)"
picard MarkDuplicates \
    I="$ALIGN_DIR/${SAMPLE}.sorted.bam" \
    O="$FILT_DIR/${SAMPLE}.markdup.bam" \
    M="$LOG_DIR/${SAMPLE}.dup_metrics.txt" \
    REMOVE_DUPLICATES=false \
    &>> "$LOG_DIR/${SAMPLE}.picard.log"

log "STEP 4/9: Filter MAPQ >= $MAPQ_MIN, proper pairs, remove dups + mito"
samtools view -@ "$THREADS" -b -f 2 -F 1284 -q "$MAPQ_MIN" \
    "$FILT_DIR/${SAMPLE}.markdup.bam" \
    | samtools sort -@ "$THREADS" -o "$FILT_DIR/${SAMPLE}.filtered.bam" -
samtools index -@ "$THREADS" "$FILT_DIR/${SAMPLE}.filtered.bam"

# Optional: drop chrM/MT reads if present
if samtools view -H "$FILT_DIR/${SAMPLE}.filtered.bam" | grep -qE 'SN:(chrM|MT)\b'; then
    MITO=$(samtools view -H "$FILT_DIR/${SAMPLE}.filtered.bam" | grep -oE 'SN:(chrM|MT)\b' | head -1 | cut -d: -f2)
    CHRS=$(samtools view -H "$FILT_DIR/${SAMPLE}.filtered.bam" | grep '^@SQ' | grep -v "SN:${MITO}\b" | sed 's/.*SN://; s/\t.*//')
    samtools view -@ "$THREADS" -b "$FILT_DIR/${SAMPLE}.filtered.bam" $CHRS \
        > "$FILT_DIR/${SAMPLE}.filtered.nomito.bam"
    mv "$FILT_DIR/${SAMPLE}.filtered.nomito.bam" "$FILT_DIR/${SAMPLE}.filtered.bam"
    samtools index -@ "$THREADS" "$FILT_DIR/${SAMPLE}.filtered.bam"
fi

samtools flagstat "$FILT_DIR/${SAMPLE}.filtered.bam" > "$LOG_DIR/${SAMPLE}.flagstat.txt"

#===============================================================================
# 5. Fragment size distribution QC
#    Expect a peak around ~147 bp (mono-nucleosome) and possibly a smaller
#    sub-nucleosomal (~50-100bp) or di-nucleosome (~300bp) shoulder depending
#    on MNase digestion extent. This is the key MNase-seq quality check.
#===============================================================================
log "STEP 5/9: Fragment size distribution (deeptools bamPEFragmentSize)"
bamPEFragmentSize \
    -b "$FILT_DIR/${SAMPLE}.filtered.bam" \
    --histogram "$FRAG_DIR/${SAMPLE}.fragment_size_hist.png" \
    --table "$FRAG_DIR/${SAMPLE}.fragment_size_table.txt" \
    -p "$THREADS" \
    &>> "$LOG_DIR/${SAMPLE}.fragsize.log"

#===============================================================================
# 6. Filter to mono-nucleosome fragment size range for nucleosome positioning
#    analysis (adjust -m/-M via CLI flags if your digestion is over/under).
#===============================================================================
log "STEP 6/9: Filtering to mono-nucleosome fragments (${MIN_INSERT}-${MAX_INSERT} bp)"
samtools view -h -@ "$THREADS" "$FILT_DIR/${SAMPLE}.filtered.bam" \
    | awk -v min="$MIN_INSERT" -v max="$MAX_INSERT" \
        'substr($0,1,1)=="@" || ($9>=min && $9<=max) || ($9<=-min && $9>=-max)' \
    | samtools sort -@ "$THREADS" -o "$NUCFILT_DIR/${SAMPLE}.mononuc.bam" -
samtools index -@ "$THREADS" "$NUCFILT_DIR/${SAMPLE}.mononuc.bam"

#===============================================================================
# 7. Coverage tracks (RPGC-normalized bigWig), centered on fragment midpoints
#    which best represent inferred nucleosome dyad positions.
#===============================================================================
log "STEP 7/9: Generating normalized bigWig coverage track"
GENOME_SIZE=$(samtools view -H "$FILT_DIR/${SAMPLE}.filtered.bam" \
    | awk '/^@SQ/{for(i=1;i<=NF;i++) if($i ~ /^LN:/){gsub("LN:","",$i); sum+=$i}} END{print sum}')

bamCoverage \
    -b "$NUCFILT_DIR/${SAMPLE}.mononuc.bam" \
    -o "$TRACK_DIR/${SAMPLE}.mononuc.bw" \
    --binSize 1 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize "$GENOME_SIZE" \
    --MNase \
    -p "$THREADS" \
    &>> "$LOG_DIR/${SAMPLE}.bamcoverage.log"

#===============================================================================
# 8. Nucleosome calling / positioning with DANPOS3
#    Produces nucleosome positions, occupancy, and fuzziness scores.
#    Skips gracefully if danpos.py is not installed.
#===============================================================================
log "STEP 8/9: Nucleosome calling (DANPOS3)"
if command -v danpos.py >/dev/null 2>&1; then
    samtools view -h "$NUCFILT_DIR/${SAMPLE}.mononuc.bam" > "$DANPOS_DIR/${SAMPLE}.mononuc.sam"
    ( cd "$DANPOS_DIR" && danpos.py dpos "${SAMPLE}.mononuc.sam" \
        -o "danpos_${SAMPLE}" &>> "../../${LOG_DIR}/${SAMPLE}.danpos.log" )
    $KEEP_TMP || rm -f "$DANPOS_DIR/${SAMPLE}.mononuc.sam"
else
    log "  danpos.py not found - skipping nucleosome calling. Install with: pip install danpos"
fi

#===============================================================================
# 9. Summary report
#===============================================================================
log "STEP 9/9: Writing summary report"
{
    echo "MNase-seq Pipeline Summary — Sample: $SAMPLE"
    echo "=============================================="
    echo "Raw reads:            $(zcat "$FASTQ_R1" | wc -l | awk '{print $1/4}')"
    echo "--- Alignment (samtools flagstat) ---"
    cat "$LOG_DIR/${SAMPLE}.flagstat.txt"
    echo "--- Duplicate metrics ---"
    grep -A1 "LIBRARY" "$LOG_DIR/${SAMPLE}.dup_metrics.txt" | head -2
    echo "--- Mono-nucleosome filtered read pairs ---"
    samtools view -c -f 64 "$NUCFILT_DIR/${SAMPLE}.mononuc.bam"
} > "$OUTDIR/${SAMPLE}.summary_report.txt"

$KEEP_TMP || rm -f "$ALIGN_DIR/${SAMPLE}.sorted.bam" "$FILT_DIR/${SAMPLE}.markdup.bam"

log "PIPELINE COMPLETE. Outputs in: $OUTDIR"
log "  Filtered BAM:       $FILT_DIR/${SAMPLE}.filtered.bam"
log "  Mono-nuc BAM:       $NUCFILT_DIR/${SAMPLE}.mononuc.bam"
log "  Fragment size QC:   $FRAG_DIR/${SAMPLE}.fragment_size_hist.png"
log "  BigWig track:       $TRACK_DIR/${SAMPLE}.mononuc.bw"
log "  Summary report:     $OUTDIR/${SAMPLE}.summary_report.txt"
