# Idea: the chr5 HERV9 provirus and the CDH6 switch at endocrine delamination

Status: hypothesis, 2026-10-04. Script drafted (`scripts/22_herv9_cdh6_trajectory.py` / `.sb`), not run yet.

## Starting point
- The LTR12_ / HERV9 endocrine signal is one full-length provirus, chr5:30,197,738-30,205,134 (+)
  (LTR12B 5' LTR - HERV9NC-int - HERV9-int - LTR12_ 3' LTR), intergenic, not segdup (script 21).
  RNA 0 from hPSC to PP, ENP1 26, ENP2 160, early EC-beta 396, early delta 252, late 29-108
  (per M gene UMIs); 5' LTR ATAC opens with it, in every endocrine unit.
- Nearest protein-coding gene: CDH6 (K-cadherin), ~1 Mb, nothing protein-coding in between
  (then DROSHA 1.2 Mb, C5orf22 1.3 Mb).
- HERV9/ERV9 is catarrhine-specific (entered ~30-40 Mya; most active in apes, Costas & Naveira 2000),
  so this copy is a candidate primate-specific addition to the endocrine program (precedent:
  ZNF808 / MER11 in pancreas fate, De Franco et al. Nat Genet 2023).

## CDH6 in our differentiation (existing tables)
| stage / group | CDH6 | source |
|---|---|---|
| hPSC | off (2% of cells, log2FC -5.4 vs rest) | diff-multiome annotation de_vs_rest |
| DE | off (0.3%, -7.8) | same |
| GT1 / GT2 | on (42%) / lower (14%) | same |
| OTX2_endo | ~139 per M | HERV locus_followup.csv |
| PP1/PP2 | peak: 81% of cells, log2FC +3.0 vs rest | annotation de_vs_rest |
| ENP1 -> ENP2 | falls; ENP1 3.6x ENP2 (log2FC +1.8) | de_ENP1_vs_ENP2.csv |
| early EC-beta / SC-EC | lowest (32 / 22 per M) | locus_followup.csv |
| early alpha / beta | 59 / 85 | same |
| early delta | partly kept (159) | same |
| late SC-beta | back up (log2FC +1.2 vs other endocrine) | deseq2_endocrine markers |

Mirror image of the provirus: CDH6 peaks in PP and drops at ENP1 -> ENP2, exactly where the provirus
switches on; the group with the most provirus (early EC-beta) has almost the least CDH6; delta keeps the
most CDH6 and has less provirus; late beta regains CDH6 as the provirus fades.

## Hypothesis
Not "the provirus activates CDH6" (they are anti-correlated), but:
1. The provirus marks, or takes part in, the epithelial-to-endocrine adhesion switch when endocrine
   progenitors delaminate from the PP epithelium (partial EMT, cadherin down-regulation downstream of
   NEUROG3 / SNAI2).
2. If it acts on CDH6 in cis, it would be repressive (or a competing / decoy element), or the two are
   co-regulated by the same endocrine TFs.
3. Speculative evolutionary angle: a catarrhine-specific element tuning progenitor adhesion could feed
   into primate islet features (mixed alpha/beta/delta architecture vs rodent core-mantle). This is a
   long step beyond the data; keep it as a framing question only.

## Tests
1. Order of events along PP -> ENP -> endocrine pseudotime: does CDH6 fall before, with, or after the
   provirus rises? Do CDH1 / EPCAM / CLDNs fall with CDH6 (general delamination) or is CDH6 special?
   -> script 22 (onset_summary.csv, trajectory.png).
2. Cis linkage on metacells: provirus RNA and 5' LTR accessibility vs every gene within +-2 Mb, against
   a random-gene null; 5' LTR vs CDH6 promoter accessibility. -> script 22 (cis_links.csv).
3. 3D contact: public islet / EndoC-betaH1 Hi-C or promoter-capture Hi-C - same TAD as CDH6?
4. Transcript structure: read coverage along the provirus (start in LTR12B? splicing? where it ends),
   to know whether it is a lncRNA that could act locally.
5. TF motifs in LTR12B (NEUROG3, NEUROD1, RFX, NKX2-2, INSM1, NF-Y) and whether they are bound /
   footprinted in ENP2.
6. Primary tissue: CDH6 and provirus in human fetal pancreas (NEUROG3+ vs trunk) and adult islets;
   orthology of the copy in macaque (if present, macaque islets are a test; mouse lacks HERV9).
7. Perturbation (later): CRISPRi of the LTR12B 5' LTR in the differentiation; read out CDH6, endocrine
   yield, and cell clustering / islet-like aggregate architecture.

## Caveats
- One copy, one differentiation protocol; stage and batch are confounded for most stages.
- Anti-correlation over a lineage is expected for many gene pairs; needs the cis / null comparison.
- CDH6 values above mix methods (Wilcoxon on log data, DESeq2, CPM unit means): directions only.
- 1 Mb is a long distance for an enhancer-promoter contact; no evidence yet that they share a TAD.

## Literature to start from
- Costas & Naveira 2000, Mol Biol Evol 17:320 (ERV9 evolutionary history).
- Belshaw et al. 2005, J Virol 79:6997 (ERV9 extinction).
- De Franco et al. 2023, Nat Genet (ZNF808 / MER11, pancreatic agenesis).
- Cabrera et al. 2006, PNAS (human islet cytoarchitecture).
- Gouzi et al. 2011, Dev Dyn; Rukstalis & Habener 2007, Gene Expr Patterns (NEUROG3 / Snail2, endocrine
  delamination). Check these before citing.

## Update 2026-10-04: the provirus is an annotated, spliced lncRNA (Ensembl)
- Ensembl (current release, source havana_tagene = GENCODE long-read pipeline) annotates
  ENSG00000297849 ("novel transcript", lncRNA, + strand) at chr5:30,197,989-30,205,018, i.e. the provirus.
  One transcript, ENST00000751379.1, 2 exons (829 nt):
    exon 1  30,197,989-30,198,214  = starts ~250 bp inside the LTR12B 5' LTR (30,197,738-30,198,183),
            ends ~30 bp into HERV9NC-int (leader / splice-donor region)
    exon 2  30,204,416-30,205,018  = end of HERV9-int into the LTR12_ 3' LTR (30,204,661-30,205,134)
  -> LTR-promoted, spliced, 5' LTR-to-3' LTR transcript (retroviral subgenomic-like splice).
  Our pipeline called it "intergenic" because the Cell Ranger 2020-A GTF (GENCODE v32) predates it.
  10x 3' reads should fall on exon 2 (3' LTR), matching where IRescue put most of the signal.
- Neighbourhood (Ensembl, chr5:29.2-31.2 Mb): no protein-coding genes until CDH6 (31,193,679);
  ~15 lncRNAs, mostly new long-read models (ENSG0000029x/30x), e.g. ENSG00000304776 (30.10 Mb, +),
  ENSG00000307515 (30.35 Mb, +), ENSG00000287940 (30.55-30.74 Mb, -), and ENSG00000289601
  (31.09-31.27 Mb, -), antisense across the CDH6 promoter / 5' end.
- PDZD2 (~1.5-2 Mb, beyond CDH6/DROSHA/C5orf22): Leung et al. 2008 Int J Biochem Cell Biol and
  2009 Stem Cells Dev report PDZD2 affects growth / differentiation of human fetal pancreatic progenitors
  and maturation of progenitor-derived islet-like clusters. Check whether it tracks the provirus.
To do: rerun locus/gene counting with a current GENCODE GTF (v47+) so the provirus and the CDH6
antisense lncRNA are counted as genes; check splice junction reads (exon1->exon2) in the GEX BAMs.

## Result 2026-10-05: script 22 run (job 33438462, scVI B latent `X_scvi_scviB_per_sample`)
- Order along diffusion pseudotime (half-max, 95% bootstrap CI): EPCAM fall 0.29, CDH1 fall 0.36,
  CDH6 fall 0.48 (0.47-0.76) together with NEUROG3 rise 0.48, NEUROD1 0.76, CHGA 0.86, CDH6 promoter
  ATAC fall 0.86, 5' LTR ATAC rise 0.89, provirus RNA rise 0.90 (0.89-0.90).
  -> CDH6 drops at NEUROG3 onset; the provirus switches on later, after CHGA (committed endocrine).
- Cis links (541 metacells): CDH6 rho -0.55 vs provirus RNA, empirical p 0.10 against 1,000 random
  genes; DROSHA -0.49 (p 0.17); 5' LTR vs CDH6 promoter accessibility -0.24. No cis evidence.
- Conclusion: the delamination / CDH6 hypothesis is not supported; the provirus is an endocrine-state
  marker. Caveat: ENP1 forms its own diffusion arm, so pseudotime between PP and endocrine is not linear.
  Report v4 section 5.11. Outputs: $TE_OUT/qc/herv9_cdh6/ (local copy results/herv9_cdh6/).
