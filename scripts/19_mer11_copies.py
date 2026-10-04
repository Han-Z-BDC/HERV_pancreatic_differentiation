#!/usr/bin/env python
"""
ZNF808 / MER11 at copy level (follows 17_znf808_mer11.py, which found no subfamily-level signal).

Do individual MER11A/B/C copies switch on in the liver-like off-target cells compared with the
pancreatic progenitors (PP1, PP2) and foregut (FG) they branch from? Uses the IRescue locus pseudobulk
($TE_OUT/irescue_locus/locus_units.npz from 16): UMIs per million gene UMIs, unit means per group.
A copy is "liver-on" if liver >= 2 CPM, log2((liver+0.5)/(max(PP1,PP2,FG)+0.5)) >= 2, and it is up
(> 2x the PP/FG maximum) in at least half of the liver units. Reported with host / nearest genes.
Comparison: the same rule for every other LTR-class subfamily, to see whether MER11 stands out.

Outputs: $TE_OUT/qc/mer11_liver_copies.csv, mer11_vs_other_subfamilies.csv
"""
import os
import subprocess
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

TE_OUT = Path(os.environ["TE_OUT"]); TE_ANNOT = Path(os.environ["TE_ANNOT_DIR"])
BEDTOOLS = str(Path(os.environ["TE_ENV"]) / "bin/bedtools")
qc = TE_OUT / "qc"; LOC = TE_OUT / "irescue_locus"
pd.set_option("display.width", 250)

M = sp.load_npz(LOC / "locus_units.npz").tocsr()
loci = pd.Index((LOC / "locus_units_loci.txt").read_text().split("\n"))
units = pd.read_csv(LOC / "locus_units_units.csv", index_col=0)
cpm = (sp.diags(1e6 / units.lib_size.values) @ M).tocsc()
grp = units.group.values


def gmean(g):
    return np.asarray(cpm[grp == g].mean(0)).ravel()


liver, pp1, pp2, fg = gmean("liver"), gmean("PP1"), gmean("PP2"), gmean("FG")
ref = np.maximum.reduce([pp1, pp2, fg])
lfc = np.log2((liver + 0.5) / (ref + 0.5))
lu = cpm[grp == "liver"].toarray()
units_up = (lu > 2 * ref[None, :]).sum(0)
n_liver = (grp == "liver").sum()
on = (liver >= 2) & (lfc >= 2) & (units_up >= n_liver / 2)
sub = pd.Series([l.split("|")[0] for l in loci], index=loci)

# per subfamily: fraction of expressed copies that are liver-on
expressed = (np.maximum(liver, ref) >= 2)
tab = pd.DataFrame({"sub": sub.values, "on": on, "expr": expressed}).groupby("sub").agg(
    copies=("on", "size"), expressed=("expr", "sum"), liver_on=("on", "sum"))
tab = tab[tab.expressed >= 20]
tab["frac_on"] = tab.liver_on / tab.expressed
tab["rank"] = tab.frac_on.rank(ascending=False)
tab.sort_values("frac_on", ascending=False).round(3).to_csv(qc / "mer11_vs_other_subfamilies.csv")
print(f"LTR subfamilies with >= 20 expressed copies: {len(tab)}; median fraction liver-on {tab.frac_on.median():.3f}")
print(tab.loc[[s for s in ["MER11A", "MER11B", "MER11C", "HERVK11-int"] if s in tab.index]].round(3).to_string())
print("top 15 subfamilies by fraction of liver-on copies:\n", tab.sort_values("frac_on", ascending=False).head(15).round(3).to_string())

mer = np.where(sub.str.match(r"^(MER11[ABC]|HERVK11-int)$").values & on)[0]
out = pd.DataFrame({"locus": loci[mer], "liver": liver[mer], "PP1": pp1[mer], "PP2": pp2[mer], "FG": fg[mer],
                    "log2fc_vs_max": lfc[mer], "liver_units_up": [f"{u}/{n_liver}" for u in units_up[mer]]})
ctx = pd.read_csv(qc / "ltr_loci_context.csv.gz", index_col=0)
out = out.join(ctx, on="locus")
if len(out):
    p = out.locus.str.split("|").str[1].str.rsplit("_", n=3, expand=True)
    bed = pd.DataFrame({"c": p[0], "s": p[1].astype(int), "e": p[2].astype(int), "id": out.locus}).sort_values(["c", "s"])
    bed.to_csv(qc / "tmp_mer11.bed", sep="\t", header=False, index=False)
    genes = []
    for line in open(TE_ANNOT / "scte/genes.gtf"):
        t = line.split("\t", 9)
        if len(t) > 8 and t[2] == "gene" and 'gene_type "protein_coding"' in t[8]:
            genes.append((t[0], int(t[3]) - 1, int(t[4]), t[8].split('gene_name "')[1].split('"')[0]))
    pd.DataFrame(genes).sort_values([0, 1]).to_csv(qc / "tmp_pc.bed", sep="\t", header=False, index=False)
    cl = pd.read_csv(StringIO(subprocess.run([BEDTOOLS, "closest", "-a", str(qc / "tmp_mer11.bed"), "-b", str(qc / "tmp_pc.bed"),
                                              "-k", "2", "-d", "-t", "all"], capture_output=True, text=True, check=True).stdout),
                     sep="\t", header=None)
    out["nearest_pc"] = out.locus.map(cl.groupby(3).apply(lambda d: "; ".join(f"{g} {int(x) // 1000}kb" for g, x in zip(d[7], d[8]))))
    (qc / "tmp_mer11.bed").unlink(); (qc / "tmp_pc.bed").unlink()
out.sort_values("log2fc_vs_max", ascending=False).round(2).to_csv(qc / "mer11_liver_copies.csv", index=False)
print(f"\nMER11A/B/C / HERVK11-int copies on in liver (vs max of PP1, PP2, FG): {len(out)}")
print(out.sort_values("liver", ascending=False).round(2).to_string(index=False))
print("done", qc)
