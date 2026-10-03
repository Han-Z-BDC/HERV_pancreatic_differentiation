#!/bin/bash
set -euo pipefail
cd /pl/active/ZhuLab_BDC/GitHub/HERV_pancreatic_differentiation
source config.sh
for m in rna atac scte_exclusive scte_nointron; do "$R_ENV/bin/Rscript" scripts/08_deseq2_te.R "$TE_OUT/deseq2/$m" 2>&1 | grep -v R_HOME; done
