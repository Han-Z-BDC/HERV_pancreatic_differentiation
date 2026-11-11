#!/usr/bin/env python
"""Annotate Olaniru 2023 on the scvi_library embedding, using the paper's own marker panels.

Clustering uses X_scvi_library (chosen in the within-study comparison: the only correction that
raised within-pool4 late-stage recovery, 0.582 -> 0.626, while lifting pool mixing 0.414 -> 0.499).
Marker scoring uses the FULL gene matrix, not the 2,000 HVGs -- several panel genes (PPY, LEPR,
CDH19, ACSL1 ...) are cell-type markers precisely because they are restricted, so they need not
survive HVG selection, and scoring on the HVG subset would silently drop them.

Panels are ported verbatim from 04_seurat_replicate.R, which took them from the paper, so these
labels stay comparable to Olaniru's own. sc.tl.score_genes is scanpy's equivalent of Seurat's
AddModuleScore: mean of the panel minus the mean of a reference set sampled from matched
expression bins.

The polyhormonal rule is the paper's: a sub-cluster is polyhormonal when >=2 of INS/GCG/SST/GHRL
each reach >=30% of that sub-cluster's top hormone and the top hormone exceeds 1.
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad

BASE = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas")
OUT = BASE / "olaniru2023" / "annotation"
RS = 1234

MAJOR = {
    "Endocrine":    ["INS", "GCG", "SST", "PAX6", "CHGA", "NEUROD1"],
    "Acinar":       ["CPA1", "PRSS1", "CLPS", "CTRB2"],
    "Ductal":       ["SLC4A4", "CFTR", "ANXA4", "SOX9", "KRT19"],
    "Mesenchymal":  ["COL3A1", "DCN", "VIM", "COL1A1"],
    "Immune":       ["RAC2", "LYZ", "TRAC", "PTPRC"],
    "Endothelial":  ["VWF", "ADGRL4", "ANGPT2", "PECAM1"],
    "Schwann":      ["CRYAB", "CDH19", "SOX10"],
    "Erythroblast": ["HBB", "HBG2", "HBA1"],
}
ENDO = {
    "EP_NEUROG3": ["NEUROG3", "NEUROD1", "INSM1", "HES6", "SOX4"],
    "Beta":       ["INS", "IAPP", "NKX6-1", "PDX1"],
    "Alpha":      ["GCG", "ARX", "IRX2", "TTR"],
    "Delta":      ["SST", "HHEX", "LEPR"],
    "Epsilon":    ["GHRL", "ACSL1"],
    "PP_cell":    ["PPY"],
}
HORMONES = ["INS", "GCG", "SST", "GHRL"]
HORMONE_MIN = 1.0    # log1p CP10K; a cell "has" a hormone above this
POLY_FRAC = 0.30     # a sub-cluster is polyhormonal when >=30% of its cells carry >=2 hormones
RES_SWEEP = [0.4, 0.6, 0.8, 1.0, 1.5]
MARGIN_MIN = 0.5     # below this, the best panel barely beats the runner-up


def load_full_with_latent():
    """Full gene matrix + the scvi_library latent, joined on the library|barcode key."""
    cmp_ = sc.read_h5ad(BASE / "olaniru2023" / "integration" / "olaniru_batch_compare.h5ad")
    lat = pd.DataFrame(cmp_.obsm["X_scvi_library"], index=cmp_.obs_names)
    meta = pd.read_csv(BASE / "qc_review" / "fetal_cell_metadata.csv.gz")
    meta = meta[(meta.in_atlas) & (meta.study == "olaniru2023")]
    parts = []
    h5d = BASE / "olaniru2023" / "qc_filtered" / "detect_doublets" / "h5ad"
    for p in sorted(h5d.glob("*.doublets_removed.h5ad")):
        if "autothresh" in p.name:
            continue
        lib = p.name.split(".")[0]
        a = sc.read_h5ad(p)
        if "soupx_counts" in a.layers:
            a.X = a.layers["soupx_counts"].copy()
        a.layers.clear(); a.obsm.clear(); a.varm.clear(); a.uns.clear()
        a.var_names_make_unique()
        m = meta[meta.library == lib].set_index("barcode")
        a = a[[b for b in a.obs_names if b in m.index]].copy()
        o = pd.DataFrame(index=a.obs_names)
        for c in ["library", "donor_label", "donor_id", "pcw", "stage", "stage_bin"]:
            o[c] = m.loc[a.obs_names, c].to_numpy()
        a.obs = o
        a.obs_names = [f"{lib}|{b}" for b in a.obs_names]
        parts.append(a)
    full = ad.concat(parts, join="inner", merge="unique")
    full = full[lat.index].copy()                      # same order as the latent
    full.obsm["X_scvi_library"] = lat.to_numpy()
    print(f"full matrix {full.shape}, latent {full.obsm['X_scvi_library'].shape}", flush=True)
    return full


def score_panels(adata, panels, prefix):
    missing = {}
    for name, genes in panels.items():
        have = [g for g in genes if g in adata.var_names]
        miss = [g for g in genes if g not in adata.var_names]
        if miss:
            missing[name] = miss
        if not have:
            print(f"  WARNING: no genes found for {name}; score set to 0", flush=True)
            adata.obs[f"{prefix}{name}"] = 0.0
            continue
        sc.tl.score_genes(adata, have, score_name=f"{prefix}{name}", random_state=RS)
    if missing:
        print(f"  panel genes absent from the matrix: {missing}", flush=True)
    return missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--endo-resolution", type=float, default=1.0)
    a_ = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    adata = load_full_with_latent()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.neighbors(adata, use_rep="X_scvi_library", random_state=RS)
    sc.tl.umap(adata, random_state=RS)
    for r in RES_SWEEP:
        sc.tl.leiden(adata, resolution=r, key_added=f"leiden_{r}", flavor="igraph",
                     n_iterations=2, directed=False, random_state=RS)
        print(f"  resolution {r}: {adata.obs[f'leiden_{r}'].nunique()} clusters", flush=True)
    adata.obs["cluster"] = adata.obs[f"leiden_{a_.resolution}"]

    print("\n=== major-type scoring (paper panels) ===", flush=True)
    score_panels(adata, MAJOR, "major_")
    cols = [f"major_{k}" for k in MAJOR]
    cl = adata.obs.groupby("cluster", observed=True)[cols].mean()
    cl.columns = list(MAJOR)
    cl["cell_type"] = cl.idxmax(axis=1)
    cl["margin"] = cl[list(MAJOR)].apply(
        lambda r: float(np.sort(r.to_numpy())[-1] - np.sort(r.to_numpy())[-2]), axis=1)
    cl["n_cells"] = adata.obs["cluster"].value_counts().reindex(cl.index).to_numpy()
    adata.obs["cell_type"] = cl["cell_type"].reindex(adata.obs["cluster"]).to_numpy()
    # a label is only as good as its margin over the runner-up; flag the weak ones and print
    # what they actually express, so an unsupported call is visible rather than silent
    cl["confident"] = cl["margin"] >= MARGIN_MIN
    cl.to_csv(OUT / "cluster_major_scores.csv")
    print(cl[["n_cells", "cell_type", "margin", "confident"]].to_string())
    weak = cl.index[~cl["confident"]].tolist()
    if weak:
        print(f"\n  {len(weak)} cluster(s) below margin {MARGIN_MIN} -- label not supported by the "
              f"panels: {weak}", flush=True)
    sc.tl.rank_genes_groups(adata, "cluster", method="wilcoxon", n_genes=15)
    top = pd.DataFrame({c: [adata.uns["rank_genes_groups"]["names"][i][c] for i in range(15)]
                        for c in cl.index})
    top.to_csv(OUT / "cluster_top_markers.csv", index=False)
    print("\n=== top markers, clusters whose label is NOT supported ===")
    for c in weak:
        print(f"  cluster {c:>2} (n={int(cl.loc[c,'n_cells']):>5}, called {cl.loc[c,'cell_type']:12} "
              f"margin {cl.loc[c,'margin']:.3f}): {', '.join(top[c][:12])}")

    # ---- endocrine sub-clustering
    en_mask = (adata.obs["cell_type"] == "Endocrine").to_numpy()
    print(f"\n=== endocrine sub-clustering ({en_mask.sum():,} cells) ===", flush=True)
    en = adata[en_mask].copy()
    sc.pp.neighbors(en, use_rep="X_scvi_library", random_state=RS)
    sc.tl.leiden(en, resolution=a_.endo_resolution, key_added="sub", flavor="igraph",
                 n_iterations=2, directed=False, random_state=RS)
    score_panels(en, ENDO, "endo_")
    ecols = [f"endo_{k}" for k in ENDO]
    es = en.obs.groupby("sub", observed=True)[ecols].mean()
    es.columns = list(ENDO)
    lab = es.idxmax(axis=1)

    # Polyhormonal, computed PER CELL rather than on cluster means.
    #
    # The paper's rule ("two hormones each >=30% of the top, top > 1") was written for Seurat's
    # AverageExpression, which back-transforms to count space; applied to mean log1p values it
    # fires far too readily -- in the first run it relabelled every Alpha sub-cluster (GCG 2.8
    # and 4.6) and a Delta one, leaving the dataset with no alpha cells at all.
    #
    # It is also wrong in principle on means: a cluster holding SEPARATE alpha and beta cells has
    # means that look co-expressing. Polyhormonality is a property of individual cells, so a cell
    # counts as polyhormonal when it carries >=2 hormones above a per-cell threshold, and a
    # sub-cluster is polyhormonal only when enough of its cells are.
    hv = [h for h in HORMONES if h in en.var_names]
    hx = pd.DataFrame(np.asarray(en[:, hv].X.todense()), index=en.obs_names, columns=hv)
    percell_pos = (hx >= HORMONE_MIN).sum(axis=1)
    en.obs["n_hormones"] = percell_pos.to_numpy()
    frac_poly = percell_pos.groupby(en.obs["sub"].to_numpy(), observed=True).apply(
        lambda v: float((v >= 2).mean()))
    hmean = hx.groupby(en.obs["sub"].to_numpy(), observed=True).mean()
    poly = frac_poly >= POLY_FRAC
    lab[poly.reindex(lab.index).fillna(False) & (lab != "EP_NEUROG3")] = "Polyhormonal"
    es_extra = pd.DataFrame({"frac_cells_2plus_hormones": frac_poly.round(3)})
    es["endo_type"] = lab
    es = es.join(hmean.add_prefix("mean_")).join(es_extra)
    es["n_cells"] = en.obs["sub"].value_counts().reindex(es.index).to_numpy()
    es.to_csv(OUT / "endocrine_subclusters.csv")
    print(es[["n_cells", "endo_type"]].to_string())

    adata.obs["endo_type"] = pd.Series(
        lab.reindex(en.obs["sub"]).to_numpy(), index=en.obs_names).reindex(adata.obs_names)
    adata.obs["endo_subcluster"] = pd.Series(
        en.obs["sub"].astype(str).to_numpy(), index=en.obs_names).reindex(adata.obs_names)
    adata.obs["cell_type_fine"] = adata.obs["endo_type"].fillna(adata.obs["cell_type"])

    um = adata.obsm["X_umap"]
    cells = adata.obs.copy()
    cells["UMAP_1"], cells["UMAP_2"] = um[:, 0], um[:, 1]
    cells.index.name = "cell"
    cells.to_csv(OUT / "cells.csv")
    adata.write_h5ad(OUT / "olaniru_annotated.h5ad")

    print("\n=== composition ===")
    print(adata.obs["cell_type"].value_counts().to_string())
    print("\n=== fine composition ===")
    print(adata.obs["cell_type_fine"].value_counts().to_string())
    print("\n=== cell type x stage ===")
    print(pd.crosstab(adata.obs["cell_type_fine"], adata.obs["stage"]).to_string())
    print(f"\nwrote {OUT}/cells.csv, cluster_major_scores.csv, endocrine_subclusters.csv, "
          f"olaniru_annotated.h5ad")


if __name__ == "__main__":
    main()
