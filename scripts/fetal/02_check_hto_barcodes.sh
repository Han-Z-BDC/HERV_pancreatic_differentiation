#!/bin/bash
# Which hashtag sequences are in each HTO library? TotalSeq-B on 10x 3' v3: the 15 bp hashtag barcode
# sits at read 2 positions 11-25. Tallies the 15-mers in the first 2 M reads of each HTO run; the top
# ones should be BioLegend TotalSeq-B Hashtags 4, 6 and 8 (paper). Run on the login node after 01.
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
declare -A HT=([AGTAAGTTCAGCGTA]=Hashtag4_B0254 [GGTTGCCAGATGTCA]=Hashtag6_B0256 [CTCCTCTGCAATTAC]=Hashtag8_B0258
               [GTCAACTCTTTAGCG]=Hashtag1 [TGATGGCCTATTGGG]=Hashtag2 [TTCCGCCTCTCTTTG]=Hashtag3
               [AAGTATCGTTTCGCA]=Hashtag5 [TGTCTTTCCTGCCAG]=Hashtag7)
tail -n +2 "$(dirname "$0")/../../configs/fetal_olaniru_samples.tsv" | while IFS=$'\t' read -r task pool pcw gsm gex hto demux; do
  f="$FETAL_OUT/raw/$hto/${hto}_2.fastq.gz"
  echo "== $pool ($hto) read 2 length: $(zcat "$f" | head -2 | tail -1 | tr -d '\n' | wc -c)"
  zcat "$f" | head -8000000 | awk 'NR % 4 == 2 { print substr($0, 11, 15) }' | sort | uniq -c | sort -k1,1nr | head -6 \
    | while read -r n s; do echo "  $n $s ${HT[$s]:-?}"; done
done
