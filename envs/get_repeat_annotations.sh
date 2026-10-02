#!/bin/bash
# Download the hg38 RepeatMasker annotation (UCSC rmsk table) into $TE_ANNOT_DIR and make:
#   rmsk_hg38.bed.gz      all repeats: chrom start end name(subfamily) score strand family class
#   rmsk_hg38_ERV.bed.gz  LTR-class repeats only (ERV1, ERVK, ERVL, ERVL-MaLR, Gypsy, unknown LTR)
# Source: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz (~160 MB)
set -euo pipefail
source "$(dirname "$0")/../config.sh"
mkdir -p "$TE_ANNOT_DIR"; cd "$TE_ANNOT_DIR"
[ -s rmsk.txt.gz ] || wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
# UCSC rmsk columns: 6 genoName 7 genoStart 8 genoEnd 10 strand 11 repName 12 repClass 13 repFamily
zcat rmsk.txt.gz | awk 'BEGIN{OFS="\t"} $6 !~ /_/ {print $6,$7,$8,$11,$2,$10,$13,$12}' \
  | sort -k1,1 -k2,2n | gzip > rmsk_hg38.bed.gz
zcat rmsk_hg38.bed.gz | awk '$8=="LTR"' | gzip > rmsk_hg38_ERV.bed.gz
echo "repeats: $(zcat rmsk_hg38.bed.gz | wc -l); LTR/ERV: $(zcat rmsk_hg38_ERV.bed.gz | wc -l)"
zcat rmsk_hg38_ERV.bed.gz | cut -f7 | sort | uniq -c | sort -k1,1nr
