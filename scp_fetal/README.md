# scp_fetal — the lab single-cell-pipeline, instanced for the fetal atlas

The pipeline at `/pl/active/ZhuLab_BDC/GitHub/single-cell-pipeline` is config-driven and meant to be
run from your own directory ("clone it into your own folder ... every path you need is set once, in
your own `config.sh`"). Its own `config.sh` points at the diff-multiome project and is left alone.

This directory is that instance for the fetal data:

    config.sh    fetal paths (IN_DIR / OUT_DIR / SAMPLE_SHEET / CONFIG_DIR / TOOLS)
    samples.tsv  sample_id <TAB> path to raw_feature_bc_matrix.h5 relative to IN_DIR
    configs/     a COPY of the pipeline's configs, so the thresholds actually used on this
                 project are recorded here and not shared with diff-multiome. qc_RNA.yaml
                 carries the fetal choice of n_mads (3, not the repo default of 5).
    scripts ->   symlink to the pipeline's scripts, so there is one copy of the code.
    logs/

Step 1.1 is NOT run from here. `scripts/fetal/09_qc_filter.py` in the HERV repo stands in for it,
because 1.1 has no notion of the erythrocyte removal this data needs (pool4 is 25% fetal
erythrocytes, which corrupt its own MAD statistics). 09 writes the same
`rna_qc/h5ad/<sample>.qc_filtered.h5ad` that 1.1 would, so 1.2 and 1.3 run unmodified.

Run (from this directory, so SLURM_SUBMIT_DIR resolves here):

    sbatch --array=0-3 scripts/1.2_remove_background.sb
    sbatch --array=0-3 scripts/1.3_detect_doublets.sb
