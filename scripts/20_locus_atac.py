#!/usr/bin/env python
"""
Chromatin at the copies behind each ERV candidate: Tn5 cut sites in each copy's region, per group.

Regions: the top copies from 18 (locus_followup.csv), each widened by 3 kb on both sides so the
flanking LTRs of a provirus are included. Cut sites come straight from the Cell Ranger fragment files
(tabix queries, kept cells only), are summed into the ATAC units and scaled by total cut sites
(per million); group = mean of its units. For each copy: accessibility in its candidate's group vs the
mean of the other groups (log2FC), and the units open (> 2x rest mean).
Also lists the LTR-class elements inside each window (the provirus structure).

Outputs: $TE_OUT/qc/locus_atac.csv
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

REPO = Path(__file__).resolve().parents[1]
TE_OUT = Path(os.environ["TE_OUT"]); ANNOT = Path(os.environ["ANNOTATION_CSV"]); TE_ANNOT = Path(os.environ["TE_ANNOT_DIR"])
qc = TE_OUT / "qc"; FLANK = 3000
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)

fu = pd.read_csv(qc / "locus_followup.csv")
p = fu.locus.str.rsplit("_", n=3, expand=True)
fu["chrom"], fu["start"], fu["end"] = p[0], p[1].astype(int), p[2].astype(int)
regions = fu[["locus", "chrom", "start", "end"]].drop_duplicates("locus").set_index("locus")
regions["w0"] = (regions.start - FLANK).clip(lower=0); regions["w1"] = regions.end + FLANK

samples = pd.read_csv(REPO / "configs/samples.tsv", sep="\t")
ann = pd.read_csv(ANNOT, index_col=0, dtype=str)
ann = ann[ann.include_downstream.str.lower() == "true"]
units = pd.read_csv(TE_OUT / "deseq2/atac/coldata.csv", index_col=0)

counts = pd.DataFrame(0.0, index=units.index, columns=regions.index)
for s in samples.itertuples():
    frag = Path(s.atac_bam).parent / "atac_fragments.tsv.gz"
    keep = {c.split("#", 1)[1]: g for c, g in ann.ATAC_WNN_annotation.items() if c.startswith(s.sample + "#")}
    tb = pysam.TabixFile(str(frag))
    for lid, r in regions.iterrows():
        n = {}
        for line in tb.fetch(r.chrom, r.w0, r.w1):
            c, a, b, bc = line.split("\t")[:4]
            g = keep.get(bc)
            if g is None:
                continue
            a, b = int(a), int(b) - 1
            hits = (r.w0 <= a < r.w1) + (r.w0 <= b < r.w1)
            u = f"{s.sample}__{g}"
            n[u] = n.get(u, 0) + hits
        for u, v in n.items():
            if u in counts.index:
                counts.loc[u, lid] += v
    print(f"{s.sample}: done", flush=True)
cpm = counts.div(units.lib_size, axis=0) * 1e6
gm = cpm.groupby(units.group).mean()

# LTR-class elements in each window
rm = pd.read_csv(TE_ANNOT / "rmsk_hg38.bed.gz", sep="\t", header=None, usecols=[0, 1, 2, 3, 5, 7],
                 names=["chrom", "start", "end", "sub", "strand", "class"])
rm = rm[rm["class"].str.startswith("LTR")]
struct = {}
for lid, r in regions.iterrows():
    e = rm[(rm.chrom == r.chrom) & (rm.end > r.w0) & (rm.start < r.w1)]
    struct[lid] = " ".join(f"{x.sub}({x.strand})" for x in e.itertuples())

rows = []
for x in fu.itertuples():
    g = x.level
    a = gm.loc[g, x.locus]; rest = gm.drop(g)[x.locus].mean()
    u = cpm.loc[units.group == g, x.locus]
    rows.append({"subfamily": x.subfamily, "level": g, "tier": x.tier, "rank": x.rank, "locus": x.locus,
                 "rna_log2fc": x.log2fc, "atac_cpm_group": a, "atac_cpm_rest": rest,
                 "atac_log2fc": np.log2((a + 0.5) / (rest + 0.5)), "atac_units_open": f"{(u > 2 * rest).sum()}/{len(u)}",
                 "atac_max_group": gm[x.locus].idxmax(), "window_ltrs": struct[x.locus], "context": x.context,
                 "host_gene": x.host_gene})
out = pd.DataFrame(rows)
out.round(3).to_csv(qc / "locus_atac.csv", index=False)
print(out.drop(columns=["window_ltrs"]).round(2).to_string(index=False))
print("done", qc)
