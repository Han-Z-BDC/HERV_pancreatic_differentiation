#!/usr/bin/env python
"""
Pseudobulk TE-subfamily counts for DESeq2, for RNA (IRescue) and ATAC (fragment cut sites).

RNA : $TE_OUT/irescue_subfamily/irescue_subfamily_merged.h5ad (from 06_merge_irescue.py)
ATAC: $TE_OUT/atac_te_subfamily/<sample>.counts.tsv.gz and .cell_totals.tsv (from 05_atac_te_counts.sb)
Cells: include_downstream == True. Units: sample x ATAC_WNN_annotation group with >= --min-cells cells.
Library size per unit (used as DESeq2 size factors): RNA = summed gene UMIs; ATAC = summed cut sites.

Outputs ($TE_OUT/deseq2/{rna,atac}/): counts.csv (subfamilies x units), coldata.csv
  (unit, sample, stage, diff_batch, group, n_cells, lib_size), subfamily_info.csv (family, class).
Also merges the ATAC per-cell counts into $TE_OUT/atac_te_subfamily/atac_te_subfamily_merged.h5ad.
"""
import argparse
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[1]
TE_OUT = Path(os.environ["TE_OUT"]); ANNOT = Path(os.environ["ANNOTATION_CSV"])
TE_ANNOT = Path(os.environ["TE_ANNOT_DIR"])
ap = argparse.ArgumentParser()
ap.add_argument("--min-cells", type=int, default=20)
args = ap.parse_args()

samples = pd.read_csv(REPO / "configs/samples.tsv", sep="\t")["sample"].tolist()
ann = pd.read_csv(ANNOT, index_col=0, dtype=str)
keep = ann.index[ann.include_downstream.str.lower() == "true"]
rm = pd.read_csv(TE_ANNOT / "rmsk_hg38.bed.gz", sep="\t", header=None, usecols=[3, 6, 7],
                 names=["subfamily", "family", "class"]).drop_duplicates("subfamily").set_index("subfamily")
rm["class"] = rm["class"].str.rstrip("?")

# ---- ATAC: merge per-sample long tables into cells x subfamilies
parts, tot = [], []
for s in samples:
    d = TE_OUT / "atac_te_subfamily"
    c = pd.read_csv(d / f"{s}.counts.tsv.gz", sep="\t", header=None, names=["barcode", "subfamily", "n"])
    c["cell"] = s + "#" + c.barcode
    parts.append(c[["cell", "subfamily", "n"]])
    t = pd.read_csv(d / f"{s}.cell_totals.tsv", sep="\t"); t.index = s + "#" + t.barcode
    tot.append(t[["total_cut_sites", "te_cut_sites"]])
long = pd.concat(parts, ignore_index=True)
tot = pd.concat(tot)
cells = pd.Index(tot.index.intersection(keep))
subs = pd.Index(sorted(long.subfamily.unique()))
long = long[long.cell.isin(cells)]
X = sp.csr_matrix((long.n.values.astype(np.float32), (cells.get_indexer(long.cell), subs.get_indexer(long.subfamily))),
                  shape=(len(cells), len(subs)))
atac = ad.AnnData(X=X, obs=ann.loc[cells, ["stage", "ATAC_WNN_annotation"]].join(tot), var=pd.DataFrame(index=subs))
atac.var = atac.var.join(rm)
atac.write_h5ad(TE_OUT / "atac_te_subfamily" / "atac_te_subfamily_merged.h5ad", compression="gzip")
print(f"ATAC merged: {atac.n_obs} cells x {atac.n_vars} subfamilies; "
      f"TE share of cut sites {100 * atac.obs.te_cut_sites.sum() / atac.obs.total_cut_sites.sum():.1f}%", flush=True)

# ---- RNA
rna = ad.read_h5ad(TE_OUT / "irescue_subfamily" / "irescue_subfamily_merged.h5ad")
rna = rna[rna.obs_names.isin(keep)].copy()


def pseudobulk(a, lib_col, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    o = a.obs.copy()
    o["sample"] = [n.split("#")[0] for n in a.obs_names]
    o["unit"] = o["sample"] + "__" + o["ATAC_WNN_annotation"].astype(str)
    n = o["unit"].value_counts()
    units = n.index[n >= args.min_cells].sort_values()
    codes = pd.Categorical(o["unit"], categories=units).codes
    m = codes >= 0
    ind = sp.csr_matrix((np.ones(m.sum()), (codes[m], np.where(m)[0])), shape=(len(units), a.n_obs))
    cnt = pd.DataFrame(np.rint((ind @ a.X).toarray()).astype(np.int64).T, index=a.var_names, columns=units)
    cd = pd.DataFrame(index=units)
    cd["sample"] = [u.split("__")[0] for u in units]
    cd["group"] = [u.split("__")[1] for u in units]
    cd["stage"] = [s.split("_")[1] for s in cd["sample"]]
    cd["diff_batch"] = [s.split("_")[0] for s in cd["sample"]]
    cd["n_cells"] = n.reindex(units).values
    cd["lib_size"] = o.groupby("unit")[lib_col].sum().reindex(units).values
    cd.index.name = "unit"
    cnt.to_csv(outdir / "counts.csv"); cd.to_csv(outdir / "coldata.csv")
    info = pd.DataFrame(index=a.var_names).join(rm); info.to_csv(outdir / "subfamily_info.csv")
    print(f"{outdir.name}: {cnt.shape[0]} subfamilies x {cnt.shape[1]} units "
          f"({cd.group.nunique()} groups; dropped {(n < args.min_cells).sum()} units under {args.min_cells} cells)", flush=True)
    return cd


cd_rna = pseudobulk(rna, "gene_umis", TE_OUT / "deseq2" / "rna")
cd_atac = pseudobulk(atac, "total_cut_sites", TE_OUT / "deseq2" / "atac")
print(cd_rna.groupby("group").size().to_string())
