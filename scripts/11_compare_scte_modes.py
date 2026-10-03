#!/usr/bin/env python
"""
Flag host-gene read-through among the TE calls: compare IRescue, scTE exclusive and scTE nointron.

For each subfamily x level (group or stage) called up vs rest (padj < 0.05, log2FC > 1) by IRescue or
scTE exclusive, ask whether the induction survives when only intergenic TE reads count (nointron):
  intergenic      nointron also up (padj < 0.05, log2FC > 1)
  weak            nointron log2FC > 0.5 but not called
  no_intergenic   nointron mean count < 10 (filtered) -- almost all signal is in genes
  intronic_only   nointron log2FC <= 0.5 -- induction carried by copies inside genes (read-through suspect)
Also: frac_level / frac_rest = nointron / exclusive counts in that level vs all other levels; a low
ratio means the induced level's TE reads come disproportionately from genic copies.

Inputs : $TE_OUT/deseq2/{rna,scte_exclusive,scte_nointron}/results/<tag>_vs_rest_all.csv, counts.csv, coldata.csv
Outputs: $TE_OUT/qc/scte_mode_comparison_{group,stage}.csv, scte_mode_lfc_scatter.png
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TE_OUT = Path(os.environ["TE_OUT"]); D = TE_OUT / "deseq2"; qc = TE_OUT / "qc"
MODES = {"ir": "rna", "ex": "scte_exclusive", "no": "scte_nointron"}
FOCUS = ["HERVS71-int", "HERV15-int", "LTR12_", "LTR12E", "MLT1F2", "LTR48B", "MER92-int", "MER110-int"]
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)


def level_counts(mode, tag):
    cnt = pd.read_csv(D / mode / "counts.csv", index_col=0)
    cd = pd.read_csv(D / mode / "coldata.csv", index_col=0).reindex(cnt.columns)
    return cnt.T.groupby(cd[tag].values).sum().T


def compare(tag):
    res = {}
    for k, m in MODES.items():
        r = pd.read_csv(D / m / "results" / f"{tag}_vs_rest_all.csv").set_index(["subfamily", "level"])
        res[k] = r[["baseMean", "log2FoldChange", "padj"]].add_suffix(f"_{k}")
        if k == "ex":
            fam = r[["family", "class"]]
    t = res["ex"].join(res["no"], how="left").join(res["ir"], how="outer").join(fam, how="left")
    for k in MODES:
        t[f"up_{k}"] = (t[f"padj_{k}"] < 0.05) & (t[f"log2FoldChange_{k}"] > 1)
    t = t.reset_index()
    ex, no = level_counts(MODES["ex"], tag), level_counts(MODES["no"], tag)
    tot_ex, tot_no = ex.sum(1), no.sum(1)
    t["frac_overall"] = (tot_no / tot_ex).reindex(t.subfamily).values
    lv_ex = ex.stack().reindex(list(zip(t.subfamily, t.level))).values
    lv_no = no.stack().reindex(list(zip(t.subfamily, t.level))).values
    with np.errstate(invalid="ignore", divide="ignore"):
        t["frac_level"] = lv_no / lv_ex
    t["frac_rest"] = (tot_no.reindex(t.subfamily).values - lv_no) / (tot_ex.reindex(t.subfamily).values - lv_ex)
    t["frac_ratio"] = t.frac_level / t.frac_rest
    called = t.up_ir | t.up_ex
    v = np.select([t.up_no, t.log2FoldChange_no > 0.5, t.log2FoldChange_no.isna()],
                  ["intergenic", "weak", "no_intergenic"], "intronic_only")
    t["verdict"] = np.where(called, v, "")
    t = t.sort_values(["level", "log2FoldChange_ex"], ascending=[True, False])
    t.round(4).to_csv(qc / f"scte_mode_comparison_{tag}.csv", index=False)

    print(f"\n===== {tag}: up vs rest in IRescue / exclusive / nointron: "
          f"{t.up_ir.sum()} / {t.up_ex.sum()} / {t.up_no.sum()}  (all classes)")
    c = t[t.verdict != ""]
    print("verdicts for calls (IRescue or exclusive), by class:\n",
          pd.crosstab(c["class"], c.verdict).loc[lambda x: x.sum(1) >= 3].to_string())
    print("agreement IRescue vs exclusive among calls:", pd.crosstab(c.up_ir, c.up_ex).to_string().replace("\n", " | "))
    ok = t.baseMean_ex.gt(20) & t.baseMean_ir.gt(20)
    print(f"log2FC correlation (baseMean > 20 in both): IRescue~exclusive r={t[ok].log2FoldChange_ir.corr(t[ok].log2FoldChange_ex):.2f}; "
          f"exclusive~nointron r={t[ok].log2FoldChange_ex.corr(t[ok].log2FoldChange_no):.2f}")
    cols = ["subfamily", "level", "family", "baseMean_ir", "log2FoldChange_ir", "log2FoldChange_ex", "log2FoldChange_no",
            "padj_no", "frac_overall", "frac_level", "frac_ratio", "verdict"]
    erv = c[(c["class"] == "LTR") & (c.baseMean_ex > 20)]
    print(f"ERV/LTR calls with baseMean_ex > 20 ({len(erv)}):")
    print(erv[cols].round(2).to_string(index=False))
    f = t[t.subfamily.isin(FOCUS) & ((t.verdict != "") | t.up_no)]
    print("focus subfamilies:\n", f[cols].round(2).to_string(index=False), flush=True)
    return t


g = compare("group")
s = compare("stage")

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
col = {"intergenic": "#1b7837", "weak": "#a6dba0", "no_intergenic": "#999999", "intronic_only": "#c51b7d"}
for ax, (x, y, ttl) in zip(axes, [("log2FoldChange_ir", "log2FoldChange_ex", "IRescue vs scTE exclusive"),
                                  ("log2FoldChange_ex", "log2FoldChange_no", "scTE exclusive vs nointron")]):
    d = g[(g["class"] == "LTR") & (g.baseMean_ex > 20)]
    ax.scatter(d[x], d[y], s=4, c="#dddddd", rasterized=True)
    for vv, cc in col.items():
        e = d[d.verdict == vv]
        ax.scatter(e[x], e[y], s=10, c=cc, label=f"{vv} ({len(e)})")
    ax.axhline(0, c="k", lw=0.5); ax.axvline(0, c="k", lw=0.5); ax.plot([-8, 8], [-8, 8], c="k", lw=0.5, ls="--")
    ax.set_xlabel(x); ax.set_ylabel(y); ax.set_title(f"ERV/LTR subfamily x group vs rest: {ttl}")
axes[1].legend(fontsize=8, frameon=False, title="verdict (calls)")
fig.tight_layout(); fig.savefig(qc / "scte_mode_lfc_scatter.png", dpi=100, bbox_inches="tight")
print("done", qc)
