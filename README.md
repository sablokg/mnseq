# mnseq

- my entire MNaseq pipeline for analysis, developed at IBCH Poland. 

```
FastQC on raw reads
Trim Galore adapter/quality trimming
Bowtie2 alignment (-X 1000, --no-mixed --no-discordant for clean fragment-size analysis)
Picard MarkDuplicates + MAPQ/proper-pair filtering + optional chrM removal
Fragment size distribution QC (bamPEFragmentSize) — the key MNase sanity check, looking for the ~147 bp mono-nucleosome peak
Mono-nucleosome filtering (keeps 100–200 bp fragments by default, tunable via -m/-M)
Normalized bigWig track generation with --MNase mode (centers signal on fragment midpoints = inferred dyads)
DANPOS3 nucleosome calling (positions, occupancy, fuzziness) 
Summary report

``` 

Gaurav Sablok \
gsablok@proton.me
