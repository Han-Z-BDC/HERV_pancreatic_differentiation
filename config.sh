#!/bin/bash
# Shared paths for HERV_pancreatic_differentiation. Source this from every job script.
# Inputs come from the diff-multiome single-cell-pipeline outputs (read only); outputs go to TE_OUT.

REPO=/pl/active/ZhuLab_BDC/GitHub/HERV_pancreatic_differentiation
SCP_REPO=/pl/active/ZhuLab_BDC/GitHub/single-cell-pipeline          # shared helper scripts
SCP_OUT=/scratch/alpine/hzhu10@xsede.org/sc_pipeline_full/out        # shared intermediate data
SAMPLE_MANIFEST=/pl/active/ZhuLab_BDC/Han_analysis/diff_multiome_2026/renamed/rename_manifest.csv
CELLRANGER_DIRS="/pl/active/ZhuLab_BDC/Han_analysis/cellranger/diff_multiome/batch1 /pl/active/ZhuLab_BDC/Han_analysis/cellranger/diff_multiome/batch2"
ANNOTATION_CSV=$SCP_OUT/annotation/scviB_x_atac_mnc_sample_graph/leiden_r0p8/ATAC_WNN_annotation.csv

GENOME_DIR=/pl/active/ZhuLab_BDC/genome
CELLRANGER_REF=$GENOME_DIR/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
BLACKLIST_BED=$GENOME_DIR/GRCh38_unified_blacklist.bed
TE_ANNOT_DIR=$GENOME_DIR/repeats_hg38                              # RepeatMasker / HERV annotations

TE_OUT=/scratch/alpine/hzhu10@xsede.org/herv_pancreatic_differentiation/out
TE_ENV=/projects/hzhu10@xsede.org/software/anaconda/envs/herv_te    # TE tools (to be built)
SC_ENV=/projects/hzhu10@xsede.org/software/anaconda/envs/single.cell # scanpy / snapatac2 side
R_ENV=/projects/hzhu10@xsede.org/software/anaconda/envs/r_env        # DESeq2

mkdir -p "$TE_OUT"

# Replication data: human fetal pancreas, Olaniru et al. 2023 Cell Metab (GSE197064), on scratch only
FETAL_OUT=/scratch/alpine/hzhu10@xsede.org/fetal_pancreas/olaniru2023
CELLRANGER=/projects/hzhu10@xsede.org/software/cellranger-9.0.1/cellranger
FETAL_ROOT=/scratch/alpine/hzhu10@xsede.org/fetal_pancreas       # all fetal datasets + the atlas (scratch only)
MIGL_OUT=/scratch/alpine/hzhu10@xsede.org/fetal_pancreas/migliorini2024  # Migliorini 2024 snRNA (GSE230403)
SRA_ENV=/projects/hzhu10@xsede.org/software/anaconda/envs/cml-env # sra-tools 3.4.1 (prefetch, fasterq-dump)
