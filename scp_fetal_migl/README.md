# single-cell-pipeline instance — Migliorini 2024 nuclei

Second instance of the lab pipeline, identical in structure to `../scp_fetal/` (Olaniru whole
cells) and pointing at the Migliorini output tree. Two instances rather than one sample sheet
because `OUT_DIR` is per-dataset and steps 1.2/1.3 are per-library anyway.

`scripts` is a symlink to the shared `single-cell-pipeline/scripts`, so the pipeline code is
never forked — only `config.sh`, `samples.tsv` and `configs/` are local.

## What differs from the Olaniru instance

Nothing in `configs/`. SoupX and the doublet callers need no cells-vs-nuclei retuning, and
`qc_RNA.yaml` is unused here: step 1.1 was replaced by `scripts/fetal/09_qc_filter.py`
(`09b_qc_filter_migliorini.sb`), which applies a **fixed 3% %MT cutoff** because the MAD rule
collapses on nuclei. See `docs/fetal_qc_cutoffs.md` for that decision and its alternatives.

## Running

    cd scp_fetal_migl
    scripts/submit_1.2_remove_background.sh     # SoupX, 6 tasks
    scripts/submit_1.3_detect_doublets.sh       # scDblFinder + Scrublet, 6 tasks

Step 1.1 is NOT run here — `09b` stands in for it and already wrote
`rna_qc/h5ad/<sample>.qc_filtered.h5ad` with the `passed_qc_step1_2` column the later steps read.
