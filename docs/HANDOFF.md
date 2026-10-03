# Handoff: HERV / TE analysis (state on 2026-10-03)

Read this first when continuing the project in a new chat.

## Question
Are endogenous retroviruses (ERVs) and other viral-origin repeats induced at specific stages or in
specific cell types during SC-islet differentiation? Data: the diff-multiome 10x multiome time course
(D0-D45, 15 samples, 98,236 kept cells, RNA + ATAC from the same nuclei).

## Where things are
- Repo (Alpine): `/pl/active/ZhuLab_BDC/GitHub/HERV_pancreatic_differentiation` (local git, not pushed).
  Local mirror on the Mac: `~/HERV_pancreatic_differentiation` (scripts, docs, small results only).
- SSH: `ssh -i ~/.ssh/id_rsa hzhu10@xsede.org@login.rc.colorado.edu`
- `config.sh` holds every path. Shared inputs from the diff-multiome `single-cell-pipeline` are read only;
  outputs go to `TE_OUT=/scratch/alpine/hzhu10@xsede.org/herv_pancreatic_differentiation/out`.
- Cell annotation: `$ANNOTATION_CSV` (columns `ATAC_WNN_annotation`, `include_downstream`; use only
  `include_downstream == True`). 24 cell groups (hPSC, DE, GT, OTX2_endo, FG, liver, PP1, PP2, ductal,
  exocrine, ENP1-3, early/late SC_beta, EC_beta, SC_alpha, SC_delta, SC_EC, stroma).
- Environments: `herv_te` (IRescue 1.2.0, scTE (patched), pychromvar, scanpy, samtools, bedtools);
  `r_env` (DESeq2 1.46). Build script `envs/build_herv_te_env.sh`.
- Annotation: `/pl/active/ZhuLab_BDC/genome/repeats_hg38/` (UCSC rmsk; `rmsk_hg38.bed.gz`,
  `rmsk_hg38_ERV.bed.gz`); scTE indexes and TE BED (4,516,048 TE loci, 1,090 subfamilies) in `scte/`.
- Samples: `configs/samples.tsv` (task, analysis name, Cell Ranger folder, GEX/ATAC BAM, cells);
  per-sample cell whitelists in `$TE_OUT/whitelists/`.

## Done
| Step | Script | Output |
|---|---|---|
| scTE indexes | `01_build_scte_indexes.sh` | `repeats_hg38/scte/` (exclusive, nointron, ATAC) |
| sample table, whitelists | `02_make_sample_table.py` | `configs/samples.tsv`, `$TE_OUT/whitelists/` |
| RNA subfamily counts (IRescue, forward-stranded, UB) | `03_irescue_subfamily.sb` | `$TE_OUT/irescue_subfamily/<sample>/counts/` (all 15) |
| RNA merge + first look | `06_merge_irescue.py` | `irescue_subfamily_merged.h5ad`, `$TE_OUT/qc/` |
| ATAC TE counts from fragments (cut sites x TE loci) | `05_atac_te_counts.sb` | `$TE_OUT/atac_te_subfamily/` (all 15) |
| pseudobulk (sample x group, >= 20 cells) | `07_pseudobulk_te.py` | `$TE_OUT/deseq2/{rna,atac}/` (119 units each) |
| DESeq2, group vs rest and stage vs rest | `08_deseq2_te.R`, `09_run_deseq2_te.sh` | `$TE_OUT/deseq2/{rna,atac}/results/` |

DESeq2 design: `~ group` (sample x group pseudobulks as replicates); size factors from library size
(gene UMIs / all cut sites). Sample blocking is NOT used because early cell types come from a single
sample each; stage and batch are therefore confounded.

## Results so far (leads, not final)
- Sanity: HERVH-int / LTR7 fall ~40x from hPSC to D9 (RNA), as expected.
- RNA: an hPSC/DE/GT/FG block (HERVH, LTR7, HERV4_I, HERVS71, HERVK11, HERVK, HERVE); HERV15-int
  in PP1/PP2 (log2FC 3.4/4.5); LTR12_/LTR12E/MLT1F2/LTR48B across ENP1-3 and early endocrine groups;
  MER92-int / MER110-int in late beta / EC-beta; HERVIP10B3 / HERVI / PABL_A in liver.
- ATAC: LTR7/LTR7Y open in hPSC and DE; LTR6A/B in DE, GT, OTX2_endo; LTR22/LTR22B in the endocrine
  groups; many LTR10 / LTR19 families in PP, OTX2_endo, liver.

## Running at handoff
- scTE RNA, all samples, `exclusive` and `nointron` indexes: job 33356420 (amem, 4 at a time;
  outputs `$TE_OUT/scte_rna_{exclusive,nointron}/<sample>.h5ad`).

## Next steps
1. Merge scTE; compare `exclusive` vs `nointron` per subfamily to flag host-gene read-through (single-
   nucleus RNA has many intronic reads). Candidates to check first: HERVS71, HERV15, LTR12_.
2. RNA vs ATAC concordance per subfamily and group (which induced families are also opened?).
3. Locus level: IRescue `-l` and/or Stellarscope (with the Telescope HERV annotation), plus strand /
   host-gene checks; name the induced copies.
4. Use the diff-multiome group peak set (pending a method choice in the pipeline chat) for TE-family
   enrichment in cell-type peaks and TF motifs inside induced LTRs (chromVAR / HOMER / ChromBPNet).
5. Exclude or flag tiny groups (ductal 52, late_SC_delta 56, exocrine 131 cells).

## Gotchas already solved
- scTE uses file mode `'rU'` (fails on Python 3.11): patched to `'r'` (in the env build script).
- scTE needs BAMs filtered to whitelisted cells with a corrected UMI:
  `samtools view -D CB:<whitelist> -e 'exists([UB])'` (done inside `04_scte_rna.sb`).
- scTE workers each load the full index: run on `amem` (6 workers, 320 GB); `acpu` 60-80 GB got OOM.
- `scTE_build -o X` writes `X.<mode>.idx`; links `hg38.<mode>.idx` exist.
- scTEATAC reads the cell barcode from the read name, not the CB tag -> unusable for Cell Ranger;
  replaced by `05_atac_te_counts.sb` (fragment files).
- `amem` jobs need QoS `mem-normal` and a large memory request (small requests are rejected);
  `acompile` needs QoS `compile`.
