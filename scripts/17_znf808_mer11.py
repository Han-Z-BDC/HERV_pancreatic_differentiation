#!/usr/bin/env python
"""
ZNF808 / MER11 check in the wild-type time course.

De Franco et al. 2023 (Nat Genet, ZNF808 pancreatic agenesis): ZNF808 binds MER11A/B/C; in ZNF808-KO
differentiations MER11 copies gain H3K27ac (some at the pancreatic progenitor stage) and progenitors
switch on liver genes (AFP, ALB, APOA1). Our differentiations make off-target liver-like cells, so ask:
across cell groups, is ZNF808 (and the KZFP co-repressors TRIM28 / SETDB1) lower where MER11 is more
transcribed or more open, in particular liver vs PP1/PP2/FG?

Gene counts: shared RNA object, layers/counts (read column-wise from the CSR without loading it all),
summed into the 119 RNA units, per million gene UMIs (lib_size of 07). MER11 signal: every TE count
table on the same units -- IRescue, scTE exclusive / nointron (RNA), cut sites / scTEATAC (ATAC).
Outputs: $TE_OUT/qc/znf808_mer11_by_group.csv, znf808_mer11.png
"""
import os
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TE_OUT = Path(os.environ["TE_OUT"]); SCP_OUT = Path(os.environ["SCP_OUT"]); ANNOT = Path(os.environ["ANNOTATION_CSV"])
D = TE_OUT / "deseq2"; qc = TE_OUT / "qc"
GENES = ["ZNF808", "TRIM28", "SETDB1", "AFP", "ALB", "APOA1", "PDX1", "NKX6-1", "INS"]
SUBS = ["MER11A", "MER11B", "MER11C", "HERVK11-int"]
ORDER = ["hPSC", "DE", "GT", "OTX2_endo", "FG", "liver", "PP1", "PP2", "ENP1", "ENP2", "ENP3", "early_SC_beta",
         "early_EC_beta", "early_SC_alpha", "early_SC_delta", "early_SC_EC", "late_SC_beta", "late_EC_beta",
         "late_SC_alpha", "late_SC_delta", "late_SC_EC", "ductal", "exocrine", "stroma"]
pd.set_option("display.width", 250)

units = pd.read_csv(D / "rna/coldata.csv", index_col=0)
ann = pd.read_csv(ANNOT, index_col=0, dtype=str)

# ---- gene counts for a few genes, from the CSR counts layer, chunked over rows
with h5py.File(SCP_OUT / "integration/h5ad/merged_rna_scvi_integrated.h5ad", "r") as f:
    cells = f["obs"]["_index"][:].astype(str)
    var = f["var"]["_index"][:].astype(str)
    cols = {g: int(np.where(var == g)[0][0]) for g in GENES}
    L = f["layers"]["counts"]; indptr = L["indptr"][:]
    gx = np.zeros((len(cells), len(GENES)), dtype=np.float32)
    colidx = np.array(list(cols.values()))
    step = 5000
    for a in range(0, len(cells), step):
        b = min(a + step, len(cells))
        idx = L["indices"][indptr[a]:indptr[b]]; dat = L["data"][indptr[a]:indptr[b]]
        row = np.repeat(np.arange(a, b), np.diff(indptr[a:b + 1]))
        for k, c in enumerate(colidx):
            m = idx == c
            np.add.at(gx[:, k], row[m], dat[m])
gx = pd.DataFrame(gx, index=cells, columns=GENES)
grp = ann["ATAC_WNN_annotation"].reindex(gx.index).astype(str)
unit = pd.Series([c.split("#")[0] for c in gx.index], index=gx.index) + "__" + grp
gu = gx.groupby(unit.values).sum().reindex(units.index)
gcpm = gu.div(units.lib_size, axis=0) * 1e6

# ---- MER11 counts on the same units from every method
te = {}
for name, m, lib in [("IRescue", "rna", "rna"), ("scTE_excl", "scte_exclusive", "rna"), ("scTE_nointron", "scte_nointron", "rna"),
                     ("ATAC_cut", "atac", "atac"), ("ATAC_scTE", "scte_atac", "atac")]:
    c = pd.read_csv(D / m / "counts.csv", index_col=0).reindex(SUBS).T
    cd = pd.read_csv(D / m / "coldata.csv", index_col=0)
    te[name] = c.div(cd.lib_size.reindex(c.index), axis=0) * 1e6
by_group = {}
G = units.group
for g, d in gcpm.groupby(G):
    row = {f"{x}_cpm": d[x].mean() for x in GENES}
    for name, t in te.items():
        tt = t.reindex(d.index)
        for s in SUBS:
            row[f"{s}|{name}"] = tt[s].mean()
    row["units"] = len(d); row["cells"] = units.loc[d.index, "n_cells"].sum()
    by_group[g] = row
bg = pd.DataFrame(by_group).T.reindex([o for o in ORDER if o in by_group])
bg.round(3).to_csv(qc / "znf808_mer11_by_group.csv")

print("genes per group (CPM, mean of units):\n", bg[[f"{x}_cpm" for x in GENES] + ["cells"]].round(1).to_string())
mer = bg[[c for c in bg.columns if c.startswith("MER11") or c.startswith("HERVK11")]]
print("\nMER11 / HERVK11 per group (per million gene UMIs or cut sites):\n", mer.round(1).to_string())
big = bg[bg.cells >= 200]
print("\nSpearman across groups (>= 200 cells): ZNF808 vs MER11 signal")
for c in mer.columns:
    print(f"  {c:28s} rho = {big['ZNF808_cpm'].corr(big[c], method='spearman'):+.2f}")
print("\nliver vs PP1/PP2/FG (ratio of group means):")
for c in ["ZNF808_cpm", "TRIM28_cpm"] + list(mer.columns):
    print(f"  {c:28s} liver/PP1 {bg.loc['liver', c] / bg.loc['PP1', c]:.2f}  liver/PP2 {bg.loc['liver', c] / bg.loc['PP2', c]:.2f}  "
          f"liver/FG {bg.loc['liver', c] / bg.loc['FG', c]:.2f}")

fig, axes = plt.subplots(1, 3, figsize=(17, 4.6), gridspec_kw={"width_ratios": [1.3, 1, 1]})
x = np.arange(len(bg))
axes[0].bar(x, bg.ZNF808_cpm, color=["#b2182b" if g == "liver" else "#2166ac" if g in ("PP1", "PP2") else "#999999" for g in bg.index])
axes[0].set_xticks(x, bg.index, rotation=70, ha="right", fontsize=8); axes[0].set_ylabel("ZNF808, UMIs per million gene UMIs")
axes[0].set_title("ZNF808 by group (liver red, PP blue)")
for ax, mod in zip(axes[1:], ["scTE_nointron", "ATAC_cut"]):
    for s, col in zip(SUBS, ["#1b9e77", "#d95f02", "#7570b3", "#444444"]):
        v = big[f"{s}|{mod}"]
        ax.scatter(big.ZNF808_cpm, v, s=18, c=col, label=s)
        for g in ("liver", "PP1", "PP2", "FG", "hPSC"):
            if g in big.index:
                ax.annotate(g, (big.loc[g, "ZNF808_cpm"], big.loc[g, f"{s}|{mod}"]), fontsize=6, alpha=.7)
    ax.set_xlabel("ZNF808 (CPM)"); ax.set_ylabel(f"MER11 signal, {mod}"); ax.set_yscale("log"); ax.set_title(f"ZNF808 vs MER11 ({mod})")
axes[1].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(qc / "znf808_mer11.png", dpi=110)
print("done", qc)
