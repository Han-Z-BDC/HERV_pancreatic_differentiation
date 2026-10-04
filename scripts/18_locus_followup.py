#!/usr/bin/env python
"""
Follow-up on the copies behind each ERV candidate (16_locus_candidates.py): neighbourhood, duplication
and whether the host / nearby genes are induced in the same group.

For the top copies of every candidate (rank <= 5 and >= 5% of the increase):
  - segmental duplication: highest fractional identity of an overlapping UCSC genomicSuperDups pair
    ($TE_ANNOT_DIR/segdups_hg38.bed); >= 0.98 means reads can be shared with a paralog, so the copy
    assignment (not the regional signal) is uncertain
  - the 3 nearest protein-coding genes (Cell Ranger GTF; distance 0 = overlapping)
  - for the host gene(s) and the nearest genes: expression in the candidate's group vs the mean of the
    other groups (UMIs per million gene UMIs, unit means, from the shared counts layer), log2FC
    -> host gene up with the copy = read-through / ERV exon; nearby gene up = possible LTR enhancer or
       promoter; nothing up = autonomous ERV transcript
Also lists clusters: copies of different candidates within 2 Mb of each other.

Outputs: $TE_OUT/qc/locus_followup.csv, $TE_OUT/qc/locus_clusters.csv
"""
import os
import subprocess
from io import StringIO
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

TE_OUT = Path(os.environ["TE_OUT"]); SCP_OUT = Path(os.environ["SCP_OUT"]); ANNOT = Path(os.environ["ANNOTATION_CSV"])
TE_ANNOT = Path(os.environ["TE_ANNOT_DIR"]); BEDTOOLS = str(Path(os.environ["TE_ENV"]) / "bin/bedtools")
qc = TE_OUT / "qc"; tmp = qc / "tmp_followup"; tmp.mkdir(exist_ok=True)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400); pd.set_option("display.max_colwidth", 70)

loc = pd.read_csv(qc / "locus_candidates_group.csv")
loc = loc[(loc["rank"] <= 5) & (loc.share_of_increase >= 0.05)].copy()
parts = loc.locus.str.rsplit("_", n=3, expand=True)
loc["chrom"], loc["start"], loc["end"], loc["strand"] = parts[0], parts[1].astype(int), parts[2].astype(int), parts[3]
copies = loc[["chrom", "start", "end", "locus"]].drop_duplicates().sort_values(["chrom", "start"])
copies.to_csv(tmp / "copies.bed", sep="\t", header=False, index=False)


def run(*args):
    return subprocess.run([BEDTOOLS, *args], capture_output=True, text=True, check=True).stdout


# ---- segmental duplications
sd = pd.read_csv(StringIO(run("intersect", "-a", str(tmp / "copies.bed"), "-b", str(TE_ANNOT / "segdups_hg38.bed"), "-wa", "-wb")),
                 sep="\t", header=None)
segdup = sd.groupby(3)[8].max() if len(sd) else pd.Series(dtype=float)
loc["segdup_identity"] = loc.locus.map(segdup)

# ---- protein-coding genes, nearest 3
rows = []
with open(TE_ANNOT / "scte/genes.gtf") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        t = line.split("\t", 9)
        if t[2] != "gene":
            continue
        a = t[8]
        name = a.split('gene_name "')[1].split('"')[0]
        gtype = a.split('gene_type "')[1].split('"')[0] if 'gene_type "' in a else ""
        rows.append((t[0], int(t[3]) - 1, int(t[4]), name, gtype, t[6]))
genes = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene", "type", "strand"])
pc = genes[genes.type == "protein_coding"].sort_values(["chrom", "start"])
pc[["chrom", "start", "end", "gene", "type", "strand"]].to_csv(tmp / "pc.bed", sep="\t", header=False, index=False)
cl = pd.read_csv(StringIO(run("closest", "-a", str(tmp / "copies.bed"), "-b", str(tmp / "pc.bed"), "-k", "3", "-d", "-t", "all")),
                 sep="\t", header=None)
near = cl.groupby(3).apply(lambda d: list(zip(d[7], d[10])))
loc["nearest_pc"] = loc.locus.map(lambda l: "; ".join(f"{g} {int(dist) // 1000}kb" for g, dist in near.get(l, [])))

# ---- gene expression per group for host + nearest genes
want = set()
for h in loc.host_gene.dropna():
    want.update(h.split(","))
for l in near.index:
    want.update(g for g, _ in near[l])
units = pd.read_csv(TE_OUT / "deseq2/rna/coldata.csv", index_col=0)
ann = pd.read_csv(ANNOT, index_col=0, dtype=str)
with h5py.File(SCP_OUT / "integration/h5ad/merged_rna_scvi_integrated.h5ad", "r") as f:
    cells = f["obs"]["_index"][:].astype(str)
    var = pd.Index(f["var"]["_index"][:].astype(str))
    gl = [g for g in sorted(want) if g in var]
    col = var.get_indexer(gl)
    pos = {c: k for k, c in enumerate(col)}
    L = f["layers"]["counts"]; indptr = L["indptr"][:]
    gx = np.zeros((len(cells), len(gl)), dtype=np.float32)
    for a in range(0, len(cells), 5000):
        b = min(a + 5000, len(cells))
        idx = L["indices"][indptr[a]:indptr[b]]; dat = L["data"][indptr[a]:indptr[b]]
        row = np.repeat(np.arange(a, b), np.diff(indptr[a:b + 1]))
        m = np.isin(idx, col)
        np.add.at(gx, (row[m], np.array([pos[c] for c in idx[m]], dtype=int)), dat[m])
gx = pd.DataFrame(gx, index=cells, columns=gl)
unit = pd.Series([c.split("#")[0] for c in cells], index=cells) + "__" + ann["ATAC_WNN_annotation"].reindex(cells).astype(str).values
gu = gx.groupby(unit.values).sum().reindex(units.index)
cpm = gu.div(units.lib_size, axis=0) * 1e6
gmean = cpm.groupby(units.group).mean()


def gene_fc(genes_, level):
    out = []
    for g in genes_:
        if g not in gmean.columns:
            continue
        a = gmean.loc[level, g]; r = gmean.drop(level)[g].mean()
        out.append(f"{g} {a:.1f}/{r:.1f} ({np.log2((a + 0.5) / (r + 0.5)):+.1f})")
    return "; ".join(out)


loc["host_expr"] = [gene_fc(h.split(","), l) if isinstance(h, str) else "" for h, l in zip(loc.host_gene, loc.level)]
loc["nearest_expr"] = [gene_fc([g for g, _ in near.get(c, [])], l) for c, l in zip(loc.locus, loc.level)]
cols = ["subfamily", "level", "tier", "reasons", "rank", "locus", "cpm_group", "cpm_rest", "log2fc", "share_of_increase",
        "units_up", "context", "host_gene", "orientation", "segdup_identity", "host_expr", "nearest_pc", "nearest_expr"]
loc[cols].round(3).to_csv(qc / "locus_followup.csv", index=False)

# ---- clusters of copies across candidates (same chromosome, within 2 Mb)
c = loc.sort_values(["chrom", "start"])
clusters, cur = [], []
for r in c.itertuples():
    if cur and (r.chrom != cur[-1].chrom or r.start - cur[-1].end > 2_000_000):
        clusters.append(cur); cur = []
    cur.append(r)
clusters.append(cur)
cl_rows = []
for k in clusters:
    cands = sorted({f"{x.subfamily}@{x.level}" for x in k})
    if len(cands) >= 3:
        cl_rows.append({"region": f"{k[0].chrom}:{min(x.start for x in k):,}-{max(x.end for x in k):,}",
                        "n_copies": len({x.locus for x in k}), "candidates": ", ".join(cands),
                        "max_segdup": np.nanmax([x.segdup_identity for x in k] + [np.nan])})
pd.DataFrame(cl_rows).to_csv(qc / "locus_clusters.csv", index=False)

show = ["subfamily", "level", "tier", "rank", "locus", "share_of_increase", "units_up", "context", "segdup_identity", "host_expr", "nearest_expr"]
print(loc[show].round(2).to_string(index=False))
print("\nclusters (>= 3 candidate x group combinations within 2 Mb):\n", pd.DataFrame(cl_rows).to_string(index=False))
for f in tmp.iterdir():
    f.unlink()
tmp.rmdir()
print("done", qc)
