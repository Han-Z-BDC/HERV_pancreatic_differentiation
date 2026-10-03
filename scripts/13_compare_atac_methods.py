#!/usr/bin/env python
"""
Compare ATAC TE-subfamily results from cut sites in fragments (05, MAPQ >= 30 only) and scTEATAC
(05b, fragment spans, all MAPQ incl. multi-mapping reads), on identical units and size factors.

  1. totals per subfamily: scTEATAC / cut-site ratio; a high ratio = signal recovered from
     multi-mapping reads (young, repetitive families)
  2. DESeq2 group / stage vs rest: log2FC correlation, call overlap (padj < 0.05, log2FC > 1)
  3. per call, agreement: both / cut_only / scte_only; the young families listed explicitly

Inputs : $TE_OUT/deseq2/{atac,scte_atac}/counts.csv and results/<tag>_vs_rest_all.csv
Outputs: $TE_OUT/qc/atac_method_totals.csv, atac_method_comparison_{group,stage}.csv, atac_method_scatter.png
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TE_OUT = Path(os.environ["TE_OUT"]); D = TE_OUT / "deseq2"; qc = TE_OUT / "qc"
YOUNG = ["LTR5_Hs", "HERVK-int", "LTR7", "LTR7Y", "HERVH-int", "LTR12C", "LTR12_", "L1HS", "L1PA2", "L1PA3",
         "SVA_D", "SVA_E", "SVA_F", "AluY", "LTR6A", "LTR6B", "LTR22", "LTR22B", "LTR15", "LTR10D"]
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)

cut = pd.read_csv(D / "atac/counts.csv", index_col=0).sum(1)
sct = pd.read_csv(D / "scte_atac/counts.csv", index_col=0).sum(1)
info = pd.read_csv(D / "atac/subfamily_info.csv", index_col=0)
tot = pd.concat([cut.rename("cut_sites"), sct.rename("scte_fragments")], axis=1).join(info)
tot["ratio"] = tot.scte_fragments / tot.cut_sites
tot.sort_values("ratio", ascending=False).to_csv(qc / "atac_method_totals.csv")
ok = tot[(tot.cut_sites > 1000) & (tot.scte_fragments > 0)]
print(f"per-subfamily totals: log10 r = {np.corrcoef(np.log10(ok.cut_sites), np.log10(ok.scte_fragments))[0, 1]:.3f}; "
      f"median scTEATAC/cut ratio {ok.ratio.median():.2f}")
print("highest ratios (gain from multi-mapping reads):\n",
      ok.sort_values("ratio", ascending=False).head(15)[["family", "class", "cut_sites", "scte_fragments", "ratio"]].round(2).to_string())
print("young / focus families:\n", tot.reindex([y for y in YOUNG if y in tot.index])[["family", "cut_sites", "scte_fragments", "ratio"]].round(2).to_string())


def compare(tag):
    r = {}
    for k, m in {"cut": "atac", "scte": "scte_atac"}.items():
        x = pd.read_csv(D / m / "results" / f"{tag}_vs_rest_all.csv").set_index(["subfamily", "level"])
        r[k] = x[["baseMean", "log2FoldChange", "padj"]].add_suffix(f"_{k}")
        if k == "cut":
            fam = x[["family", "class"]]
    t = r["cut"].join(r["scte"], how="outer").join(fam, how="left")
    for k in r:
        t[f"up_{k}"] = (t[f"padj_{k}"] < 0.05) & (t[f"log2FoldChange_{k}"] > 1)
    t["agreement"] = np.select([t.up_cut & t.up_scte, t.up_cut, t.up_scte], ["both", "cut_only", "scte_only"], "")
    t = t.reset_index().sort_values(["level", "log2FoldChange_cut"], ascending=[True, False])
    t.round(4).to_csv(qc / f"atac_method_comparison_{tag}.csv", index=False)
    b = t[(t.baseMean_cut > 20) & (t.baseMean_scte > 20)]
    print(f"\n===== {tag}: calls cut / scTEATAC / both = {t.up_cut.sum()} / {t.up_scte.sum()} / {(t.agreement == 'both').sum()}; "
          f"log2FC r = {b.log2FoldChange_cut.corr(b.log2FoldChange_scte):.2f} (baseMean > 20 in both)")
    c = t[t.agreement != ""]
    print("agreement by class:\n", pd.crosstab(c["class"], c.agreement).loc[lambda x: x.sum(1) >= 5].to_string())
    print("per level (LTR calls):\n", pd.crosstab(c[c["class"] == "LTR"].level, c[c["class"] == "LTR"].agreement).to_string())
    cols = ["subfamily", "level", "log2FoldChange_cut", "padj_cut", "log2FoldChange_scte", "padj_scte", "agreement"]
    one = c[(c["class"] == "LTR") & (c.agreement != "both") & (c[["baseMean_cut", "baseMean_scte"]].max(1) > 50)]
    print(f"LTR calls found by one method only (baseMean > 50): {len(one)}")
    print(one.sort_values(["agreement", "level"])[cols].round(2).head(60).to_string(index=False))
    y = t[t.subfamily.isin(YOUNG) & (t.up_cut | t.up_scte)]
    print("young / focus families, any call:\n", y[cols].round(2).to_string(index=False), flush=True)
    return t


g = compare("group")
compare("stage")

d = g[(g["class"] == "LTR") & (g.baseMean_cut > 20)]
fig, ax = plt.subplots(figsize=(6.5, 6))
ax.scatter(d.log2FoldChange_cut, d.log2FoldChange_scte, s=4, c="#cccccc", rasterized=True)
for a, col in {"both": "#1b7837", "cut_only": "#2166ac", "scte_only": "#b2182b"}.items():
    e = d[d.agreement == a]
    ax.scatter(e.log2FoldChange_cut, e.log2FoldChange_scte, s=9, c=col, label=f"{a} ({len(e)})")
lim = np.nanmax(np.abs(d[["log2FoldChange_cut", "log2FoldChange_scte"]].values)) * 1.05
ax.plot([-lim, lim], [-lim, lim], c="k", lw=0.5, ls="--"); ax.axhline(0, c="k", lw=0.5); ax.axvline(0, c="k", lw=0.5)
ax.set_xlabel("log2FC, cut sites (fragments, MAPQ >= 30)"); ax.set_ylabel("log2FC, scTEATAC (fragment spans, all MAPQ)")
ax.set_title("ATAC: ERV/LTR subfamily x group vs rest"); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(qc / "atac_method_scatter.png", dpi=110)
print("done", qc)
