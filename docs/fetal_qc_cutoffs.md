# Fetal pancreas atlas — QC cutoff matrix

Every threshold applied, or proposed, to the 10 fetal libraries. Rules follow the lab pipeline's
`scripts/1.1_rna_qc.py`:

- **counts / genes** — one-sided MAD: `median(log1p(x)) − n_mads · MAD(log1p(x))`, low tail only.
  No upper bound is applied anywhere; multiplets are step 1.3's job.
- **%MT** — `median + n_mads · MAD`, clamped to `[floor, cap]`; or a fixed value; or off.
- **%Hb** — fixed, applied *before* the MAD statistics so erythrocytes cannot inflate the MAD.
- **%ribo** — off.

---

## Olaniru 2023 — whole cells — APPLIED (approved 2026-10)

`n_mads: 3` · `mt_pct_n_mads: 5` · `mt_pct_floor: 3.0` · `mt_pct_cap: 15.0` · `hb_pct_max: 50`

| pool | called | erythro. removed | UMI cut | gene cut | %MT cut | no donor label | kept | % of called |
|---|---|---|---|---|---|---|---|---|
| pool1 | 8,228 | 106 | 806.0 | 537.5 | 11.88 | 773 | 6,732 | 81.8 |
| pool2 | 5,503 | 85 | 987.5 | 655.0 | 14.13 | 190 | 4,532 | 82.4 |
| pool3 | 6,775 | 27 | 1,044.5 | 711.4 | 13.69 | 323 | 5,560 | 82.1 |
| pool4 | 7,069 | **1,796** | 1,000.0 | 609.7 | 6.68 | 136 | 4,563 | 64.5 |
| | **27,575** | 2,014 | | | | 1,422 | **21,387** | **77.6** |

The %MT floor of 3.0 **never binds** here — every raw threshold (6.68–14.13) sits above it.
Erythrocyte removal is what makes a single `n_mads: 3` work across all four pools: pool4's
1,796 hemoglobin-dominated cells inflated its gene MAD to 0.609 vs 0.380 once removed.

After SoupX (1.2) and doublet removal (1.3): **19,830 cells**.

## Migliorini 2024 — nuclei — APPLIED (approved 2026-10-09)

`n_mads: 3` on counts/genes · **fixed %MT 3.0** · `hb_pct_max: 50` (no-op) · `ribo_max: None`

| sample | called | UMI cut (3 MAD) | gene cut (3 MAD) | lost to counts+genes | % |
|---|---|---|---|---|---|
| S1_14PCW | 14,072 | 1,852.9 | 1,320.1 | 562 | 3.99 |
| S2_14PCW | 12,092 | 1,294.2 | 1,020.4 | 517 | 4.28 |
| S3_16PCW | 10,676 | 2,745.0 | 1,772.9 | 490 | 4.59 |
| S4_16PCW | 13,748 | 1,896.7 | 1,391.4 | 964 | 7.01 |
| S5_18PCW | 12,030 | 1,784.5 | 1,343.6 | 838 | 6.97 |
| S6_18PCW | 13,831 | 1,797.3 | 1,310.3 | 803 | 5.81 |
| | **76,449** | | | 4,174 | 5.5 |

Erythrocyte removal is a **no-op** for nuclei: 0–2 nuclei per sample exceed 50% Hb, 7–25 exceed
10%. Mature red cells are anucleate and nucleated erythroid precursors are essentially absent.
The step stays in the script anyway, costing ~1 nucleus in total.

### The %MT decision for nuclei — settled at a fixed 3%

Median %MT is **0.02–0.10%** (vs 1.3–4.8% for whole cells) — mitochondria leave with the
discarded cytoplasm. Median %ribo is 0.50–0.81% vs 16.6–18.8%. Four candidate rules:

| rule | threshold | nuclei lost | verdict |
|---|---|---|---|
| MAD, `n_mads: 3`, no floor | 0.080 – 0.320% | 14.3 – 20.2% | rejected |
| MAD, `n_mads: 5`, no floor | 0.120 – 0.465% | 8.4 – 13.8% | rejected |
| fixed 1% | 1.0% | 0.43 – 6.14% | rejected |
| **fixed 3%** | **3.0%** | **0.05 – 2.03%** | **APPLIED** |
| fixed 5% | 5.0% | 0.03 – 1.16% | considered |
| off | — | 0 | considered |

Applied via `09_qc_filter.py --mt-fixed 3.0`, which bypasses the MAD rule outright rather than
relying on the `floor` clamp to coincidentally land on 3.0. Per-sample cost: S1 29 (0.21%),
S2 6 (0.05%), S3 134 (1.26%), S4 125 (0.91%), S5 207 (1.72%), S6 281 (2.03%).

Cost at each candidate fixed line, for the record:

| sample | PCW | median %MT | 1% | 2% | **3%** | 5% |
|---|---|---|---|---|---|---|
| S1_14PCW | 14 | 0.044% | 150 (1.07%) | 53 (0.38%) | **29 (0.21%)** | 17 (0.12%) |
| S2_14PCW | 14 | 0.020% | 52 (0.43%) | 11 (0.09%) | **6 (0.05%)** | 4 (0.03%) |
| S3_16PCW | 16 | 0.077% | 444 (4.16%) | 208 (1.95%) | **134 (1.26%)** | 70 (0.66%) |
| S4_16PCW | 16 | 0.060% | 447 (3.25%) | 211 (1.53%) | **125 (0.91%)** | 76 (0.55%) |
| S5_18PCW | 18 | 0.093% | 669 (5.56%) | 314 (2.61%) | **207 (1.72%)** | 96 (0.80%) |
| S6_18PCW | 18 | 0.103% | 849 (6.14%) | 440 (3.18%) | **281 (2.03%)** | 160 (1.16%) |

**Known limitation of any fixed %MT line here.** Median %MT is monotonic in gestational age
(0.020, 0.044 → 0.060, 0.077 → 0.093, 0.103) and so is the cost: at 3% it runs 0.05–0.21% at
14 PCW, 0.91–1.26% at 16 PCW, 1.72–2.03% at 18 PCW — about a tenfold spread with no overlap
between stages. The filter therefore thins the later stages hardest, along the axis the atlas
is built to measure. At 3% the absolute cost is small enough (max 2.03%) that this is a recorded
caveat rather than a confound, but it should be checked if a late-stage population looks thin.

The discarded nuclei are also enriched for the smallest cell types at every threshold tested
(`PTPRC`+ immune 1.45–2.29x, `PECAM1`+ endothelial 1.50–3.37x at the 1% line), because a small
nucleus receives the same absolute ambient dose but has less of its own RNA to dilute it.

**Why the unclamped MAD rules fail.** MAD is scale-free, so on a distribution centred at 0.05%
it collapses. At these depths it is also *quantization noise*: S2's median %MT of 0.0199% on
5,080 median UMI is exactly **one** mitochondrial read, and its MAD is also one read. The rule
reduces to "discard any nucleus with more than ~4–20 MT reads" — Poisson noise, not quality.

Two measurements confirm the cut would remove biology, not debris:

1. **The failing nuclei are shallower, not dirtier.** Median genes in the 5-MAD failing set is
   1,860–3,003 vs 2,554–3,385 passing. Those are good nuclei; they fail because one MT read is a
   larger fraction of a smaller denominator — and we already filter depth directly.
2. **The cut is lineage-correlated.** `PTPRC`⁺ immune and `PECAM1`⁺ endothelial nuclei are
   enriched **1.5–2.3× in the discarded set in all six samples**, while endocrine and acinar
   lineages are slightly depleted (0.89–0.99). The direction matters: lower depth normally
   *reduces* marker detection, so the enrichment is compositional, overcoming a depth handicap
   rather than caused by one. Immune and endothelial are the two rarest fetal pancreas
   populations and the ones an atlas most needs intact.

A fixed cutoff is immune to this failure mode — one stray MT read cannot reach 1% at 5,000+ UMI.
The open question for a fixed 1% is whether its cost is *even across samples*: the 99th
percentile of %MT ranges from 0.60% (S2) to 5.53% (S6), nearly tenfold, so one line will be
nearly free in S1/S2 and bite several percent in S5/S6. An unevenly applied filter becomes a
batch effect that integration then tries to correct.

---

## Not applied, deliberately

| | why |
|---|---|
| upper UMI / gene bound | one-sided by design; multiplets are 1.3's job, and an upper cut deletes large, active cells |
| %ribo cutoff | `ribo_max: None` — no evidence it separates anything here, and nuclei sit at <1% |
| paper cutoffs | Olaniru's (>500 UMI, >250 genes) removes 0–23 nuclei per sample (0.00–0.17%); inert at nuclei depth |
| cell-cycle regression | would remove proliferating progenitors, which are biology in a developmental atlas |
