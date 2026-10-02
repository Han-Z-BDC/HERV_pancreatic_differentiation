#!/bin/bash
# Build the scTE indexes (subfamily level) from RepeatMasker hg38 + the Cell Ranger ARC GRCh38 GTF.
#   TE BED (6 columns, name = subfamily): transposable-element classes only
#     (LTR, LINE, SINE, DNA, Retroposon, RC, incl. "?"-uncertain calls); simple repeats,
#     low complexity, satellites and small-RNA genes are left out.
#   RNA indexes (scTE_build): exclusive (drop reads on exons of coding genes / lncRNAs) and
#     nointron (also drop reads inside introns; intergenic TE reads only) -- comparing the two
#     separates autonomous TE transcription from host-gene read-through in single-nucleus data.
#   ATAC index (scTEATAC_build): TE BED only.
# Outputs: $TE_ANNOT_DIR/scte/  (hg38_TE.bed, hg38.exclusive.idx, hg38.nointron.idx, hg38.te_atac.idx)
set -euo pipefail
source "$(dirname "$0")/../config.sh"
out="$TE_ANNOT_DIR/scte"; mkdir -p "$out"; cd "$out"

zcat "$TE_ANNOT_DIR/rmsk_hg38.bed.gz" \
  | awk 'BEGIN{OFS="\t"} {c=$8; sub(/\?$/,"",c)} c ~ /^(LTR|LINE|SINE|DNA|Retroposon|RC)$/ {print $1,$2,$3,$4,$5,$6}' \
  > hg38_TE.bed
echo "TE loci: $(wc -l < hg38_TE.bed); subfamilies: $(cut -f4 hg38_TE.bed | sort -u | wc -l)"

gtf="$out/genes.gtf"
[ -s "$gtf" ] || zcat "$CELLRANGER_REF/genes/genes.gtf.gz" > "$gtf"

for mode in exclusive nointron; do
  "$TE_ENV/bin/scTE_build" -te hg38_TE.bed -gene "$gtf" -m "$mode" -g other -o "hg38.$mode"
done
"$TE_ENV/bin/scTEATAC_build" -g hg38_TE.bed -o hg38.te_atac
ls -la "$out"
