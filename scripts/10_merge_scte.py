#!/usr/bin/env python
"""
Merge the scTE RNA outputs (exclusive and nointron indexes) and pseudobulk them like the IRescue counts.

scTE writes genes and TE subfamilies into one cells x features matrix; TE columns are the subfamily
names of the scTE TE BED ($TE_ANNOT_DIR/scte/hg38_TE.bed). The two indexes differ in where a read on
a TE is kept as a TE read:
  exclusive : TE reads on exons of genes go to the gene; intronic and intergenic TE reads count as TE
  nointron  : TE reads anywhere in a gene body go to the gene; only intergenic TE reads count as TE
so nointron / exclusive per subfamily = share of its TE signal from copies outside genes. A subfamily
whose induction rests on intronic copies (host-gene pre-mRNA, common in single-nucleus RNA) loses it
in nointron.

Inputs : $TE_OUT/scte_rna_{exclusive,nointron}/<sample>.h5ad (04_scte_rna.sb)
         shared annotation (include_downstream == True cells only) and per-cell gene UMIs
Outputs: $TE_OUT/scte_rna_<mode>/scte_<mode>_merged.h5ad       (kept cells x TE subfamilies)
         $TE_OUT/deseq2/scte_<mode>/{counts,coldata,subfamily_info}.csv
             same units as $TE_OUT/deseq2/rna (IRescue); lib_size = gene UMIs of the cells scTE kept
         $TE_OUT/qc/scte_mode_totals.csv          per subfamily: IRescue / exclusive / nointron totals
         $TE_OUT/qc/scte_intergenic_fraction_by_group.csv   nointron / exclusive per subfamily x group
"""
import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[1]
TE_OUT = Path(os.environ["TE_OUT"]); ANNOT = Path(os.environ["ANNOTATION_CSV"])
SCP_OUT = Path(os.environ["SCP_OUT"]); TE_ANNOT = Path(os.environ["TE_ANNOT_DIR"])
MODES = ["exclusive", "nointron"]
qc = TE_OUT / "qc"; qc.mkdir(parents=True, exist_ok=True)

samples = pd.read_csv(REPO / "configs/samples.tsv", sep="\t")["sample"].tolist()
ann = pd.read_csv(ANNOT, index_col=0, dtype=str)
ann = ann[ann.include_downstream.str.lower() == "true"]
te_names = set(pd.read_csv(TE_ANNOT / "scte/hg38_TE.bed", sep="\t", header=None, usecols=[3])[3].unique())
rm = pd.read_csv(TE_ANNOT / "rmsk_hg38.bed.gz", sep="\t", header=None, usecols=[3, 6, 7],
                 names=["subfamily", "family", "class"]).drop_duplicates("subfamily").set_index("subfamily")
rm["class"] = rm["class"].str.rstrip("?")
with h5py.File(SCP_OUT / "integration/h5ad/merged_rna_scvi_integrated.h5ad", "r") as h:
    gene_umis = pd.Series(h["obs"]["total_counts"][:], index=h["obs"]["_index"][:].astype(str))
units_ref = pd.read_csv(TE_OUT / "deseq2/rna/coldata.csv", index_col=0)   # IRescue units


def merge(mode):
    blocks, names, gcounts, subs = [], [], [], None
    for s in samples:
        a = ad.read_h5ad(TE_OUT / f"scte_rna_{mode}" / f"{s}.h5ad")
        a.obs_names = s + "#" + a.obs_names.astype(str)
        a = a[a.obs_names.isin(ann.index)]
        is_te = a.var_names.isin(te_names)
        if subs is None:
            subs = a.var_names[is_te]
        X = a.X.tocsr()
        gcounts.append(np.asarray(X[:, ~is_te].sum(1)).ravel())
        blocks.append(X[:, a.var_names.get_indexer(subs)])
        names += list(a.obs_names)
        n_ann = ann.index.str.startswith(s + "#").sum()
        print(f"{mode} {s}: {a.n_obs}/{n_ann} kept cells in scTE output, {is_te.sum()} TE columns, "
              f"TE share of counts {100 * blocks[-1].sum() / X.sum():.1f}%", flush=True)
    te = ad.AnnData(X=sp.vstack(blocks).tocsr().astype(np.float32), var=pd.DataFrame(index=subs))
    te.obs_names = names
    te.obs = ann.loc[names, ["stage", "ATAC_WNN_annotation"]].copy()
    te.obs["sample"] = [n.split("#")[0] for n in names]
    te.obs["scte_gene_counts"] = np.concatenate(gcounts)
    te.obs["gene_umis"] = gene_umis.reindex(names).values
    te.var = te.var.join(rm)
    te.write_h5ad(TE_OUT / f"scte_rna_{mode}" / f"scte_{mode}_merged.h5ad", compression="gzip")
    return te


def pseudobulk(te, mode):
    o = te.obs
    unit = (o["sample"].astype(str) + "__" + o["ATAC_WNN_annotation"].astype(str)).values
    units = units_ref.index
    codes = pd.Categorical(unit, categories=units).codes
    m = codes >= 0
    ind = sp.csr_matrix((np.ones(m.sum()), (codes[m], np.where(m)[0])), shape=(len(units), te.n_obs))
    cnt = pd.DataFrame(np.rint((ind @ te.X).toarray()).astype(np.int64).T, index=te.var_names, columns=units)
    cd = units_ref.copy()
    cd["n_cells"] = pd.Series(unit).value_counts().reindex(units).fillna(0).astype(int).values
    cd["lib_size"] = pd.Series(o["gene_umis"].values, index=unit).groupby(level=0).sum().reindex(units).values
    out = TE_OUT / "deseq2" / f"scte_{mode}"; out.mkdir(parents=True, exist_ok=True)
    cnt.to_csv(out / "counts.csv"); cd.to_csv(out / "coldata.csv")
    te.var[["family", "class"]].to_csv(out / "subfamily_info.csv")
    lost = (cd.n_cells / units_ref.n_cells)
    print(f"{mode}: {cnt.shape[0]} subfamilies x {cnt.shape[1]} units; cells kept per unit vs IRescue: "
          f"median {lost.median():.3f}, min {lost.min():.3f} ({lost.idxmin()})", flush=True)
    return cnt


cnts = {m: pseudobulk(merge(m), m) for m in MODES}
ir = pd.read_csv(TE_OUT / "deseq2/rna/counts.csv", index_col=0)
grp = units_ref["group"]

tot = pd.DataFrame({"irescue": ir.sum(1), "exclusive": cnts["exclusive"].sum(1), "nointron": cnts["nointron"].sum(1)})
tot = tot.join(rm)
tot["intergenic_frac"] = tot.nointron / tot.exclusive
tot["exclusive_vs_irescue"] = tot.exclusive / tot.irescue
tot.sort_values("exclusive", ascending=False).to_csv(qc / "scte_mode_totals.csv")
by_grp = {m: cnts[m].T.groupby(grp.reindex(cnts[m].columns).values).sum().T for m in MODES}
frac = (by_grp["nointron"] / by_grp["exclusive"].where(by_grp["exclusive"] >= 20)).round(3)
frac.to_csv(qc / "scte_intergenic_fraction_by_group.csv")

pd.set_option("display.width", 250)
s = tot.groupby("class")[["irescue", "exclusive", "nointron"]].sum()
s["intergenic_frac"] = s.nointron / s.exclusive
print("totals by class:\n", s.sort_values("exclusive", ascending=False).head(8).round(3).to_string())
focus = ["HERVH-int", "LTR7", "HERVS71-int", "HERV15-int", "LTR12_", "LTR12E", "MLT1F2", "LTR48B", "MER92-int",
         "MER110-int", "HERVK11-int", "HERVK-int", "HERVE-int", "HERV4_I-int", "HERVIP10B3-int", "HERVI-int", "PABL_A-int"]
print("focus subfamilies:\n", tot.reindex([f for f in focus if f in tot.index])
      [["family", "irescue", "exclusive", "nointron", "intergenic_frac", "exclusive_vs_irescue"]].round(3).to_string())
print("done", flush=True)
