#!/usr/bin/env python
"""
The HERV9 provirus at chr5:30,197,738-30,205,134 (+): RNA and chromatin across cell groups.

Structure (RepeatMasker): LTR12B (5' LTR) - HERV9NC-int - HERV9-int - LTR12_ (3' LTR), ~7.4 kb,
1 Mb from the nearest protein-coding gene (CDH6). The per-copy analysis (16/18/20) showed this one copy
gives ~90-97% of the LTR12_ / HERV9 induction in ENP2 and early endocrine groups.

RNA : IRescue locus counts (sense strand) for the copies of this provirus, per million gene UMIs.
ATAC: Tn5 cut sites from the fragment files in 5 windows (5' LTR +-100 bp, internal 5' / 3' halves,
      3' LTR +-100 bp, 2 kb downstream), per million total cut sites; units -> group means.
Outputs: $TE_OUT/qc/herv9_chr5_provirus.csv, herv9_chr5_provirus.png
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[1]
TE_OUT = Path(os.environ["TE_OUT"]); ANNOT = Path(os.environ["ANNOTATION_CSV"]); qc = TE_OUT / "qc"
CHROM = "chr5"
WINDOWS = {"5' LTR (LTR12B)": (30197638, 30198283), "internal 5'": (30198283, 30201000),
           "internal 3'": (30201000, 30204661), "3' LTR (LTR12_)": (30204561, 30205234),
           "downstream 2 kb": (30205234, 30207234)}
COPIES = ["LTR12B|chr5_30197738_30198183_+", "HERV9NC-int|chr5_30198183_30203886_+",
          "HERV9-int|chr5_30203884_30204661_+", "LTR12_|chr5_30204661_30205134_+"]
ORDER = ["hPSC", "DE", "GT", "OTX2_endo", "FG", "liver", "PP1", "PP2", "ENP1", "ENP2", "ENP3", "early_SC_beta",
         "early_EC_beta", "early_SC_alpha", "early_SC_delta", "early_SC_EC", "late_SC_beta", "late_EC_beta",
         "late_SC_alpha", "late_SC_EC"]

# ---- RNA, per copy
M = sp.load_npz(TE_OUT / "irescue_locus/locus_units.npz").tocsc()
loci = pd.Index((TE_OUT / "irescue_locus/locus_units_loci.txt").read_text().split("\n"))
ru = pd.read_csv(TE_OUT / "irescue_locus/locus_units_units.csv", index_col=0)
idx = loci.get_indexer(COPIES)
rna = pd.DataFrame(M[:, idx[idx >= 0]].toarray(), index=ru.index, columns=[c for c, i in zip(COPIES, idx) if i >= 0])
rna = rna.div(ru.lib_size, axis=0) * 1e6
rna["provirus total"] = rna.sum(1)
rna_g = rna.groupby(ru.group).mean()

# ---- ATAC, cut sites per window
ann = pd.read_csv(ANNOT, index_col=0, dtype=str)
ann = ann[ann.include_downstream.str.lower() == "true"]
au = pd.read_csv(TE_OUT / "deseq2/atac/coldata.csv", index_col=0)
cnt = pd.DataFrame(0.0, index=au.index, columns=list(WINDOWS))
lo = min(w[0] for w in WINDOWS.values()) - 1000; hi = max(w[1] for w in WINDOWS.values()) + 1000
for s in pd.read_csv(REPO / "configs/samples.tsv", sep="\t").itertuples():
    keep = {c.split("#", 1)[1]: g for c, g in ann.ATAC_WNN_annotation.items() if c.startswith(s.sample + "#")}
    tb = pysam.TabixFile(str(Path(s.atac_bam).parent / "atac_fragments.tsv.gz"))
    for line in tb.fetch(CHROM, lo, hi):
        c, a, b, bc = line.split("\t")[:4]
        g = keep.get(bc)
        u = f"{s.sample}__{g}"
        if g is None or u not in cnt.index:
            continue
        for x in (int(a), int(b) - 1):
            for k, (w0, w1) in WINDOWS.items():
                if w0 <= x < w1:
                    cnt.loc[u, k] += 1
atac = cnt.div(au.lib_size, axis=0) * 1e6
atac_g = atac.groupby(au.group).mean()
n_open = (atac.groupby(au.group)["5' LTR (LTR12B)"].apply(lambda v: f"{(v > 2).sum()}/{len(v)}"))

tab = rna_g.add_prefix("RNA ").join(atac_g.add_prefix("ATAC ")).join(n_open.rename("units 5'LTR > 2")).reindex(ORDER)
tab.round(2).to_csv(qc / "herv9_chr5_provirus.csv")
pd.set_option("display.width", 250)
print(tab.round(2).to_string())

fig, ax1 = plt.subplots(figsize=(12, 4.2))
x = np.arange(len(ORDER))
ax1.bar(x - 0.2, tab["RNA provirus total"], width=0.4, color="#2166ac", label="RNA: provirus copies (sense UMIs per million gene UMIs)")
ax1.set_ylabel("RNA (per million gene UMIs)", color="#2166ac")
ax2 = ax1.twinx()
ax2.bar(x + 0.2, tab["ATAC 5' LTR (LTR12B)"], width=0.4, color="#b2182b", label="ATAC: 5' LTR cut sites (per million)")
ax2.set_ylabel("ATAC, 5' LTR (per million cut sites)", color="#b2182b")
ax1.set_xticks(x, ORDER, rotation=60, ha="right", fontsize=8)
ax1.set_title("HERV9 provirus chr5:30,197,738-30,205,134 (+): LTR12B - HERV9NC/HERV9-int - LTR12_")
h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, fontsize=8, frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(qc / "herv9_chr5_provirus.png", dpi=110)
print("done", qc)
