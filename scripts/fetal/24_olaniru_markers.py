#!/usr/bin/env python
"""Material for annotating Olaniru 2023 by hand: clusters at several resolutions + marker panels.

No labels are assigned anywhere in this script. It produces the things needed to assign them:

  clusters_by_resolution.png   UMAP at leiden 0.2 / 0.4 / 0.6 / 0.8 / 1.0 / 1.5 / 2.0, each cluster
                               numbered in place, so a resolution can be chosen by eye
  markers_<panel>.png          per-gene expression on the UMAP, grouped by lineage
  dotplot_res<r>.png           marker x cluster dot plots at each resolution -- usually the fastest
                               route from clusters to labels
  qc_on_umap.png               UMI / genes / %MT / pool / stage / donor, to spot clusters that are
                               technical rather than biological
  cluster_top_markers_res<r>.csv   top 25 Wilcoxon markers per cluster
  cluster_composition_res<r>.csv   cluster x pool and cluster x stage counts

Embedding is X_scvi_library, chosen in the within-study comparison. Expression is log1p CP10K on
the FULL gene matrix, not the 2,000 HVGs -- restricted markers like PPY, GHRL, LEPR and CDH19 need
not survive HVG selection, and several of the panels below would be silently empty without this.
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas")
OUT = BASE / "olaniru2023" / "annotation"
RS = 1234
RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0]
DOTPLOT_RES = [0.4, 0.8, 1.0, 1.5]

PANELS = {
  "progenitor_trunk": ["SOX9","NKX6-1","PDX1","HNF1B","ONECUT1","SPP1","HES1","CPM","GATA4","FOXA2"],
  "ductal":           ["KRT19","CFTR","ANXA4","SLC4A4","KRT18","MUC1","TACSTD2","CLDN4"],
  "acinar_tip":       ["PTF1A","CPA1","PRSS1","CTRB2","CLPS","CELA3A","RBPJL","CEL","CPA2"],
  "endocrine_prog":   ["NEUROG3","FEV","NEUROD1","INSM1","HES6","SOX4","PAX4","NKX2-2","CHGA","CHGB"],
  "endocrine_hormone":["INS","GCG","SST","GHRL","PPY","IAPP","TTR","PCSK1","PCSK2","SCG2"],
  "endocrine_tf":     ["PAX6","ARX","IRX2","MAFB","HHEX","LEPR","ISL1","NEUROD1","RFX6","PDX1"],
  "mesenchyme":       ["COL1A1","COL3A1","DCN","LUM","PDGFRA","PDGFRB","VIM","ACTA2","RGS5","POSTN"],
  "endothelial":      ["PECAM1","VWF","CDH5","KDR","ADGRL4","ANGPT2","CLDN5","ESAM"],
  "immune_myeloid":   ["PTPRC","LYZ","CD68","C1QA","C1QB","AIF1","CSF1R","ITGAM","MRC1","S100A8"],
  "immune_lymphoid":  ["TRAC","CD3E","CD3D","IL7R","MS4A1","CD79A","NKG7","GNLY","RAC2","CCL5"],
  "erythroid":        ["HBB","HBG1","HBG2","HBA1","HBA2","GYPA","ALAS2","AHSP","SLC4A1","KLF1"],
  "neural_schwann":   ["SOX10","PLP1","S100B","CRYAB","CDH19","MPZ","ERBB3","PHOX2B","ELAVL4","STMN2"],
  "other":            ["MKI67","TOP2A","UPK3B","MSLN","KRT8","EPCAM","AFP","ALB","TTR","LGALS1"],
}


def load():
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
        for c in ["library", "donor_label", "donor_id", "pcw", "stage"]:
            o[c] = m.loc[a.obs_names, c].to_numpy()
        a.obs = o
        a.obs_names = [f"{lib}|{b}" for b in a.obs_names]
        parts.append(a)
    full = ad.concat(parts, join="inner", merge="unique")[lat.index].copy()
    full.obsm["X_scvi_library"] = lat.to_numpy()
    return full


def label_clusters(ax, U, labels):
    """Put each cluster's number at its own median position, so clusters can be read off the plot."""
    for c in pd.unique(labels):
        m = (labels == c)
        ax.text(np.median(U[m, 0]), np.median(U[m, 1]), str(c), fontsize=8, weight="bold",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.72))


def main():
    ap = argparse.ArgumentParser()
    ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    adata = load()
    print(f"loaded {adata.shape}", flush=True)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars=[], percent_top=None, log1p=False, inplace=True)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    sc.pp.neighbors(adata, use_rep="X_scvi_library", random_state=RS)
    sc.tl.umap(adata, random_state=RS)
    U = adata.obsm["X_umap"]
    for r in RESOLUTIONS:
        sc.tl.leiden(adata, resolution=r, key_added=f"res{r}", flavor="igraph", n_iterations=2,
                     directed=False, random_state=RS)
        print(f"  res {r}: {adata.obs[f'res{r}'].nunique()} clusters", flush=True)

    # ---- clusters at each resolution
    n = len(RESOLUTIONS); ncol = 4; nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.0*ncol, 4.6*nrow), squeeze=False)
    for i, r in enumerate(RESOLUTIONS):
        ax = axes[i//ncol][i%ncol]
        lab = adata.obs[f"res{r}"].astype(str).to_numpy()
        uq = sorted(pd.unique(lab), key=int)
        cmap = plt.get_cmap("tab20")
        for k, u in enumerate(uq):
            m = (lab == u)
            ax.scatter(U[m,0], U[m,1], s=1.2, color=cmap(k % 20), linewidths=0, rasterized=True)
        label_clusters(ax, U, lab)
        ax.set_title(f"leiden {r}  ({len(uq)} clusters)", fontsize=12)
        ax.set_xticks([]); ax.set_yticks([])
    for j in range(n, nrow*ncol): axes[j//ncol][j%ncol].axis("off")
    fig.suptitle("Olaniru 2023 — Leiden resolution sweep on X_scvi_library", fontsize=14, y=1.0)
    fig.savefig(OUT / "clusters_by_resolution.png", dpi=115, bbox_inches="tight"); plt.close(fig)

    # ---- QC and covariates on the same UMAP, to spot technical clusters
    qcs = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    cats = ["library", "stage", "donor_id"]
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.2))
    for j, q in enumerate(qcs):
        ax = axes[0][j]
        v = adata.obs[q].to_numpy()
        vv = np.clip(v, None, np.percentile(v, 99))
        sctr = ax.scatter(U[:,0], U[:,1], c=vv, s=1.0, cmap="viridis", linewidths=0, rasterized=True)
        ax.set_title(q, fontsize=11); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(sctr, ax=ax, fraction=0.045)
    for j, c in enumerate(cats):
        ax = axes[1][j]
        lab = adata.obs[c].astype(str).to_numpy(); uq = sorted(pd.unique(lab))
        cmap = plt.get_cmap("tab20" if len(uq) > 10 else "tab10")
        for k, u in enumerate(uq):
            m = (lab == u)
            ax.scatter(U[m,0], U[m,1], s=1.0, color=cmap(k % 20), label=u, linewidths=0, rasterized=True)
        ax.set_title(c, fontsize=11); ax.set_xticks([]); ax.set_yticks([])
        if len(uq) <= 12:
            ax.legend(fontsize=6, markerscale=5, frameon=False, loc="best")
    fig.suptitle("QC and covariates on the same embedding", fontsize=13, y=1.0)
    fig.savefig(OUT / "qc_on_umap.png", dpi=115, bbox_inches="tight"); plt.close(fig)

    # ---- marker panels
    missing_all = {}
    for panel, genes in PANELS.items():
        have = [g for g in dict.fromkeys(genes) if g in adata.var_names]
        miss = [g for g in genes if g not in adata.var_names]
        if miss: missing_all[panel] = miss
        if not have: continue
        ncol = 5; nrow = int(np.ceil(len(have)/ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.5*ncol, 3.2*nrow), squeeze=False)
        for i, g in enumerate(have):
            ax = axes[i//ncol][i%ncol]
            v = np.asarray(adata[:, g].X.todense()).ravel()
            o = np.argsort(v)                       # draw positives last so they are not buried
            s = ax.scatter(U[o,0], U[o,1], c=v[o], s=1.0, cmap="Reds", linewidths=0,
                           rasterized=True, vmin=0, vmax=max(np.percentile(v, 99.5), 1e-6))
            ax.set_title(g, fontsize=11); ax.set_xticks([]); ax.set_yticks([])
            fig.colorbar(s, ax=ax, fraction=0.045)
        for j in range(len(have), nrow*ncol): axes[j//ncol][j%ncol].axis("off")
        fig.suptitle(f"{panel}  —  log1p CP10K", fontsize=13, y=1.0)
        fig.savefig(OUT / f"markers_{panel}.png", dpi=110, bbox_inches="tight"); plt.close(fig)
        print(f"  panel {panel}: {len(have)} genes plotted", flush=True)
    if missing_all:
        print(f"\n  genes absent from the matrix: {missing_all}", flush=True)

    # ---- dot plots + per-cluster markers and composition
    flat = [g for p in PANELS.values() for g in p]
    flat = [g for g in dict.fromkeys(flat) if g in adata.var_names]
    for r in DOTPLOT_RES:
        key = f"res{r}"
        sc.pl.dotplot(adata, flat, groupby=key, standard_scale="var", show=False,
                      figsize=(0.26*len(flat)+3, 0.34*adata.obs[key].nunique()+2))
        plt.savefig(OUT / f"dotplot_res{r}.png", dpi=110, bbox_inches="tight"); plt.close()
        sc.tl.rank_genes_groups(adata, key, method="wilcoxon", n_genes=25)
        names = adata.uns["rank_genes_groups"]["names"]
        pd.DataFrame({c: [names[i][c] for i in range(25)] for c in names.dtype.names}).to_csv(
            OUT / f"cluster_top_markers_{key}.csv", index=False)
        comp = pd.concat([pd.crosstab(adata.obs[key], adata.obs["library"]),
                          pd.crosstab(adata.obs[key], adata.obs["stage"])], axis=1)
        comp["n"] = adata.obs[key].value_counts().reindex(comp.index)
        comp.to_csv(OUT / f"cluster_composition_{key}.csv")
        print(f"  dotplot + markers + composition for {key}", flush=True)

    adata.write_h5ad(OUT / "olaniru_for_annotation.h5ad")
    print(f"\nwrote figures and tables to {OUT}")
    print("resolutions available in .obs: " + ", ".join(f"res{r}" for r in RESOLUTIONS))


if __name__ == "__main__":
    main()
