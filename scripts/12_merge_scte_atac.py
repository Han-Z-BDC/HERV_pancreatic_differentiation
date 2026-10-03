#!/usr/bin/env python
"""
Merge the scTEATAC outputs (05b_scte_atac.sb) and pseudobulk them into the ATAC units of 07, so DESeq2
(08) runs on exactly the same units and size factors (total cut sites from 05) as the cut-site counts.

Inputs : $TE_OUT/scte_atac/<sample>.h5ad (cells x TE subfamilies; fragments overlapping each subfamily)
         $TE_OUT/deseq2/atac/coldata.csv (units, lib_size), shared annotation, cell totals from 05
Outputs: $TE_OUT/scte_atac/scte_atac_merged.h5ad
         $TE_OUT/deseq2/scte_atac/{counts,coldata,subfamily_info}.csv
"""
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[1]
TE_OUT = Path(os.environ["TE_OUT"]); ANNOT = Path(os.environ["ANNOTATION_CSV"])
TE_ANNOT = Path(os.environ["TE_ANNOT_DIR"])

samples = pd.read_csv(REPO / "configs/samples.tsv", sep="\t")["sample"].tolist()
ann = pd.read_csv(ANNOT, index_col=0, dtype=str)
ann = ann[ann.include_downstream.str.lower() == "true"]
rm = pd.read_csv(TE_ANNOT / "rmsk_hg38.bed.gz", sep="\t", header=None, usecols=[3, 6, 7],
                 names=["subfamily", "family", "class"]).drop_duplicates("subfamily").set_index("subfamily")
rm["class"] = rm["class"].str.rstrip("?")
units_ref = pd.read_csv(TE_OUT / "deseq2/atac/coldata.csv", index_col=0)

blocks, names, subs = [], [], None
for s in samples:
    a = ad.read_h5ad(TE_OUT / "scte_atac" / f"{s}.h5ad")
    a.obs_names = s + "#" + a.obs_names.astype(str)
    a = a[a.obs_names.isin(ann.index)]
    if subs is None:
        subs = a.var_names
    blocks.append(sp.csr_matrix(a.X)[:, a.var_names.get_indexer(subs)])
    names += list(a.obs_names)
    tot = pd.read_csv(TE_OUT / "atac_te_subfamily" / f"{s}.cell_totals.tsv", sep="\t")
    print(f"{s}: {a.n_obs}/{ann.index.str.startswith(s + '#').sum()} kept cells, "
          f"{blocks[-1].sum():,.0f} TE fragment hits vs {tot.te_cut_sites.sum():,} TE cut sites", flush=True)
te = ad.AnnData(X=sp.vstack(blocks).tocsr().astype(np.float32), var=pd.DataFrame(index=subs))
te.obs_names = names
te.obs = ann.loc[names, ["stage", "ATAC_WNN_annotation"]].copy()
te.obs["sample"] = [n.split("#")[0] for n in names]
te.var = te.var.join(rm)
te.write_h5ad(TE_OUT / "scte_atac" / "scte_atac_merged.h5ad", compression="gzip")

unit = (te.obs["sample"].astype(str) + "__" + te.obs["ATAC_WNN_annotation"].astype(str)).values
units = units_ref.index
codes = pd.Categorical(unit, categories=units).codes
m = codes >= 0
ind = sp.csr_matrix((np.ones(m.sum()), (codes[m], np.where(m)[0])), shape=(len(units), te.n_obs))
cnt = pd.DataFrame(np.rint((ind @ te.X).toarray()).astype(np.int64).T, index=te.var_names, columns=units)
cd = units_ref.copy()          # lib_size = total cut sites (05), as for the cut-site counts
n = pd.Series(unit).value_counts().reindex(units).fillna(0).astype(int)
assert (n.values == cd.n_cells.values).all(), "cells per unit differ from the cut-site units"
out = TE_OUT / "deseq2" / "scte_atac"; out.mkdir(parents=True, exist_ok=True)
cnt.to_csv(out / "counts.csv"); cd.to_csv(out / "coldata.csv")
te.var[["family", "class"]].to_csv(out / "subfamily_info.csv")
print(f"scte_atac: {cnt.shape[0]} subfamilies x {cnt.shape[1]} units", flush=True)
