#!/bin/bash
# single-cell-pipeline config for the fetal pancreas atlas. See README.md.
IN_DIR=/scratch/alpine/hzhu10@xsede.org/fetal_pancreas/olaniru2023/cellranger
OUT_DIR=/scratch/alpine/hzhu10@xsede.org/fetal_pancreas/olaniru2023/qc_filtered
SAMPLE_SHEET=/pl/active/ZhuLab_BDC/GitHub/HERV_pancreatic_differentiation/scp_fetal/samples.tsv
CONFIG_DIR=/pl/active/ZhuLab_BDC/GitHub/HERV_pancreatic_differentiation/scp_fetal/configs
TOOLS=/projects/hzhu10@xsede.org/software/anaconda/envs/single.cell/bin
BLACKLIST_BED=/pl/active/ZhuLab_BDC/genome/GRCh38_unified_blacklist.bed

# the single.cell env ships a newer libstdc++ than the system one, and rpy2/anndata2ri (SoupX,
# scDblFinder) need R_HOME -- same two exports the pipeline's own config.sh sets, for the same reason
_ENV_LIB="/projects/hzhu10@xsede.org/software/anaconda/envs/single.cell/lib"
export LD_LIBRARY_PATH="${_ENV_LIB}:${LD_LIBRARY_PATH:-}"
export R_HOME="${_ENV_LIB}/R"
