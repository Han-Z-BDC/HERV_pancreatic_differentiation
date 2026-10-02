# HERV_pancreatic_differentiation

Are endogenous retroviruses (ERVs) and other viral-origin repeats induced at specific stages or in
specific cell types during stem-cell-derived islet (SC-islet) differentiation?

Data: the diff-multiome 2026 10x multiome time course (D0-D45, 15 samples, 98,629 cells; RNA + ATAC
from the same nuclei). This repo does not copy that data. It reads the shared inputs through
`config.sh` and reuses helpers from the `single-cell-pipeline` repo next to it:

| Input | Where | Made by |
|---|---|---|
| Cell annotation (`ATAC_WNN_annotation`, `include_downstream`) | `$SCP_OUT/annotation/scviB_x_atac_mnc_sample_graph/leiden_r0p8/ATAC_WNN_annotation.csv` | single-cell-pipeline 7.2 |
| scVI latent space / UMAP (98,236 kept cells) | `$SCP_OUT/integration/scvi_filtered/` | single-cell-pipeline 7.5 |
| WNN graph and UMAP | `$SCP_OUT/wnn/scviB_x_atac_mnc_sample_graph/` | single-cell-pipeline 6.3 |
| Per-group ATAC peaks | `$SCP_OUT/peaks/wnn_groups/` | single-cell-pipeline 7.6 |
| Per-sample ATAC fragments (SnapATAC2) | `$SCP_OUT/merged/h5ad/<sample>.atac_final.h5ad` | single-cell-pipeline |
| Cell Ranger ARC outputs (GEX/ATAC BAMs, fragments) | `$CELLRANGER_DIRS` | Cell Ranger ARC |

Outputs go to `$TE_OUT` (scratch), never into the shared pipeline folders.

## Plan

1. **Subfamily level (first).**
   - RNA: TE counts per cell at subfamily level with IRescue (UMI-aware EM for multi-mapping reads),
     cross-checked with scTE, from the Cell Ranger GEX BAMs. Annotation: RepeatMasker hg38.
   - ATAC: TE-subfamily accessibility per cell (chromVAR-style deviations over each subfamily's
     loci) and enrichment of subfamilies in the cell-type peak sets.
   - Statistics: pseudobulk (sample x cell type) DESeq2 with sample as the blocking term, as in
     single-cell-pipeline 7.8.
2. **Locus level (second).** Stellarscope (single-cell Telescope) with the curated HERV/L1
   annotation, for individual induced copies; strand and host-gene checks to separate autonomous LTR
   transcription from intronic read-through (single-nucleus RNA has many intronic reads).

## Known pitfalls
- Young ERVs (HERV-K HML-2, young HERV-H) are multi-mapping; discarding multi-mappers removes them.
- Intronic ERVs pick up host-gene pre-mRNA reads in single-nucleus data.
- 10x 3' chemistry captures the 3' end of polyadenylated transcripts only.
- Stage and differentiation batch coincide for most stages (only D7/D9 have two batches).

Nothing is pushed to a remote yet.
