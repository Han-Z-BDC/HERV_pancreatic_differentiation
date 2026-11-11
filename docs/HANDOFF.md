# Handoff: HERV / TE analysis (state on 2026-10-03, updated after the scTE mode comparison)

Read this first when continuing the project in a new chat.

## Question
Are endogenous retroviruses (ERVs) and other viral-origin repeats induced at specific stages or in
specific cell types during SC-islet differentiation? Data: the diff-multiome 10x multiome time course
(D0-D45, 15 samples, 98,236 kept cells, RNA + ATAC from the same nuclei).

## Where things are
- Repo (Alpine): `/pl/active/ZhuLab_BDC/GitHub/HERV_pancreatic_differentiation` (local git, not pushed).
  Local mirror on the Mac: `~/Documents/Claude/HERV_pancreatic_differentiation` (scripts, docs, small results only).
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
| scTE RNA (exclusive + nointron), all 15 | `04_scte_rna.sb` (jobs 33356420, 33356418_10) | `$TE_OUT/scte_rna_{exclusive,nointron}/` |
| scTE merge, pseudobulk (same 119 units), DESeq2 | `10_merge_scte.py` / `.sb` | `scte_<mode>_merged.h5ad`, `$TE_OUT/deseq2/scte_{exclusive,nointron}/` |
| exclusive vs nointron read-through flags | `11_compare_scte_modes.py` | `$TE_OUT/qc/scte_mode_*` (local copy `results/scte_modes/`) |
| scTEATAC, fixed (CB-prefixed read names, name-sorted BAM) | `05b_scte_atac.sb` | `$TE_OUT/scte_atac/<sample>.h5ad` (all 15) |
| scTEATAC merge + DESeq2 + compare with cut sites | `12_merge_scte_atac.py`, `12_scte_atac_deseq2.sb`, `13_compare_atac_methods.py` | `$TE_OUT/deseq2/scte_atac/`, `$TE_OUT/qc/atac_method_*` |
| candidate tiers (RNA x read-through x ATAC; genome-derived int-LTR pairs) | `14_candidate_tiers.py` | `$TE_OUT/qc/candidate_tiers_{group,stage}.csv`, `erv_int_ltr_pairs.csv` |
| ZNF808 / MER11 check | `17_znf808_mer11.py` / `.sb` | `$TE_OUT/qc/znf808_mer11*` |
| progress report (HTML, all scripts embedded) | `docs/build_report.py` + `docs/report_template.html` | `docs/HERV_progress_report.html` |
| full result explorer (HTML, interactive) | `scripts/23_build_explorer.py` + `docs/explorer_template.html` | `docs/HERV_explorer.html` |

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

### scTE exclusive vs nointron (step 1, done 2026-10-03)
nointron keeps only TE reads outside gene bodies, so nointron / exclusive = intergenic share of a
subfamily's signal (overall ~10-20%; LTR class 21%). In exclusive, TE columns take 50-78% of all scTE
counts (intronic reads). IRescue vs scTE-exclusive group log2FC r = 0.82; exclusive vs nointron r = 0.45.
Each IRescue/exclusive call is labelled intergenic / weak / no_intergenic / intronic_only, plus
`frac_ratio` = intergenic share in the called group / in the other groups (low = genic copies drive it).
- Robust (also up in nointron): HERVH-int, LTR7/7Y, HERVS71-int, HERVK11-int, LTR6B in hPSC-DE-GT-FG-
  OTX2_endo; LTR12_ / LTR12E / HERV9(NC)-int in ENP1-3 and early endocrine (intergenic share rises to
  0.8-0.9 there); MER110-int and MLT1F2-int in ENP1-3; LTR60B in early beta / EC-beta; liver
  HERVIP10B3, HERVI, PABL_A/B, LTR10B1/B2/E.
- Read-through suspects (lost in nointron): MER92-int and MER110-int in late beta / EC-beta (and D22/D45;
  MER110's ENP signal is intergenic, its late-beta signal is not); LTR48B (ENP, D15); LTR19C and LTR25-int
  in PP2; HERVK-int in OTX2_endo; HERVE-int in hPSC; HERV4_I-int at D0 (stage level); LTR21A in liver/stroma.
- HERV15-int in PP1/PP2: mixed. Still up in nointron (log2FC 1.5 / 1.7, padj < 0.05) but far below IRescue
  (3.4 / 4.5); in PP2 only 9% of its reads are intergenic vs 34% elsewhere (frac_ratio 0.26), so most of
  the PP2 boost is genic. Needs locus level.
- HERVS71-int: robust (intergenic share ~0.5 everywhere).

## Next steps
### Since then (evening 2026-10-03)
- scTEATAC fixed: it took the barcode from the read name (Cell Ranger: sequencer ID -> 2 "cells") and
  needs a name-sorted BAM for bamToBed -bedpe. Now agrees with cut sites (group log2FC r = 0.95;
  one-method calls are near-misses or the 52-cell ductal group); multi-mappers only add L1HS/SVA/Alu.
  Cut sites stay the main ATAC measure; ATAC support = either method.
- Calling rule agreed with the user: no single main RNA method. candidate = up in IRescue OR scTE-excl;
  robust = other method log2FC > 0.5 and not intronic_only in nointron; else discordant + reason
  (read_through / scte_only (likely antisense) / irescue_only / absent). Group: 242 candidates, 170 robust
  (44 with ATAC support), 72 discordant.
- Literature (user pointed out): LTR6A/B / HERVS71 in definitive endoderm is KNOWN -- Wu et al. Cell Rep
  2022 (LTR6B DE enhancer, FOXA2/GATA4, TET1, ERBB4) and Babarinde et al. Nat Commun 2025 (LTR6A
  endoderm transcripts). Treat as positive control, not a finding. De Franco et al. Nat Genet 2023:
  ZNF808 represses MER11A/B/C; KO -> pancreatic agenesis with liver genes. Our check: liver-like cells
  have ZNF808 1.4-1.6x HIGHER than PP and no MER11 gain at subfamily level (per-copy check pending).
- Apparently new (keyword search only): HERV15/LTR15 in PP (partly genic), LTR12_/LTR12E/HERV9 in
  ENP/early endocrine (intergenic, RNA only), LTR22/LTR22B opened in endocrine (ATAC only), MLT1F2-int
  in ENP, HERVIP10B3/LTR10B and PABL in the liver-like cells. MER110-int ENP = scTE-only (antisense?).
- Locus level DONE (2026-10-04): IRescue -l on all 725,763 LTR-class copies (15), copies per candidate
  (16), segdup / neighbour genes / host-gene expression (18), per-copy ATAC from fragments (20), MER11
  per copy (19), chr5 HERV9 provirus detail (21).
  HEADLINE: LTR12_/HERV9 endocrine signal = ONE full-length HERV9 provirus chr5:30,197,738-30,205,134 (+)
  (LTR12B 5' LTR - HERV9NC-int - HERV9-int - LTR12_ 3' LTR), intergenic (CDH6 1 Mb, down), not segdup.
  RNA 0 hPSC..FG, ENP1 26, ENP2 160, early EC-beta 396, early delta 252, late 29-108 (per M gene UMIs);
  5' LTR ATAC closed to PP2, opens ENP1/2, early EC-beta 16 per M cut sites; all units.
  Others: HERV15 PP = copies in ACSM1/ACSM3 (ACSM3 up 9x in PP2; LTR15 there opened) + lncRNA next to
  KLRK1; ENP LTR12E/MER110/MLT1F2 = chr5q13 SMN segdup (>=99.7%; copy unresolvable, no ATAC reads);
  MER92/MER110 late beta = introns of lncRNA AC112770.1 (up 4x); LTR19C/LTR25 PP2 = ANXA4 intron;
  LTR48B = HIST1H4I exon; HERVK OTX2 = PCAT14 exons. MER11 per copy: only 2 liver-on copies -> no
  ZNF808-loss signature. Report v3 has all of it (section 5.10).

1. (done) scTE exclusive vs nointron, ATAC methods, tiers, locus level -- see above.
   NEXT on the chr5 HERV9 provirus: read coverage along it (where transcripts start/end, splicing), public
   islet / fetal pancreas data, TF motifs in LTR12B (NEUROG3/NEUROD1/RFX/NKX2-2, NF-Y), co-accessible genes.
   HERV15/ACSM3: does the PP ACSM3 isoform start in LTR15?
2. RNA vs ATAC concordance per subfamily and group (which induced families are also opened?).
3. Locus level: IRescue `-l` and/or Stellarscope (with the Telescope HERV annotation), plus strand /
   host-gene checks; name the induced copies. Priority: HERV15-int (PP), MER92-int / MER110-int (late
   beta; find the host genes), LTR12_ (intergenic, endocrine).
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
- scTEATAC: prefix read names with CB and name-sort (05b); -f 2 -F 0xD0C, all MAPQ kept.
- `config.sh` does not export: Python jobs need `set -a; source config.sh; set +a`.
- `amem` jobs need QoS `mem-normal` and a large memory request (small requests are rejected);
  `acompile` needs QoS `compile`.
