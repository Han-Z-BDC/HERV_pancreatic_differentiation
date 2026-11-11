#!/usr/bin/env python
"""Within-study batch correction for Olaniru 2023, before annotating it.

Olaniru differs from Migliorini in the one way that decides this: pool IS crossed with stage.
PCW12 appears in pool1+pool2, PCW13 in pool1+pool3, PCW14 and PCW19 in pool2+pool3 -- 15,584 of
19,830 cells (79%). So a pool effect can be estimated separately from a developmental one, which
was impossible in Migliorini where every library held exactly one stage.

Two things are measured rather than assumed:

1. IS CORRECTION NEEDED AT ALL?  For annotation we only need cell-type clusters. If Leiden on the
   uncorrected embedding already mixes pools within clusters, correcting adds risk for nothing.
   Measured as NMI(cluster, pool) against NMI(cluster, stage) -- if clusters track pool more than
   stage, there is a pool effect worth removing.

2. DOES CORRECTION DAMAGE POOL4?  PCW15/18/20 exist only in pool4 (4,246 cells, 21%), so for the
   late stages pool and stage are perfectly confounded. Aligning pool4 to the others could flatten
   the end of the developmental trajectory. Scored two ways: stage recovery computed WITHIN pool4
   (which cannot be helped by knowing the pool), and whether pool4 stays separable at all.

donor_id is included only as a negative control: 0 of 12 Olaniru donors span more than one stage,
so it is expected to destroy stage signal exactly as it did in the joint atlas.
"""
import argparse
from pathlib import Path

import numpy as np, pandas as pd, scanpy as sc, anndata as ad

BASE = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas")
OUT = BASE / "olaniru2023" / "integration"
RS, N_HVG, N_PCS, N_LATENT, TARGET_SUM = 1234, 2000, 30, 30, 1e4
KEYS = ["library", "donor_id"]


def build():
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
        obs = pd.DataFrame(index=a.obs_names)
        for c in ["library", "donor_label", "donor_id", "pcw", "stage", "stage_bin"]:
            obs[c] = m.loc[a.obs_names, c].to_numpy()
        a.obs = obs
        a.obs_names = [f"{lib}|{b}" for b in a.obs_names]
        parts.append(a)
        print(f"  {lib}: {a.n_obs:,}", flush=True)
    x = ad.concat(parts, join="inner", merge="unique")
    print(f"Olaniru merged: {x.n_obs:,} x {x.n_vars:,}")
    return x


def nmi(a, b):
    from sklearn.metrics import normalized_mutual_info_score
    return float(normalized_mutual_info_score(a, b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(OUT))
    a_ = ap.parse_args()
    out = Path(a_.outdir); out.mkdir(parents=True, exist_ok=True)

    adata = build()
    adata.layers["counts"] = adata.X.copy()
    # HVGs re-selected WITHIN Olaniru -- the joint atlas gene set was chosen to span both studies
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=N_HVG,
                                batch_key="library", subset=True)
    print(f"HVGs: {adata.n_vars} (seurat_v3, batch_key=library, within Olaniru)")
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, n_comps=N_PCS, random_state=RS)

    from harmonypy import run_harmony, compute_lisi
    embs = {"uncorrected": adata.obsm["X_pca"]}
    for key in KEYS:
        ho = run_harmony(adata.obsm["X_pca"], adata.obs, [key], random_state=RS)
        z = np.asarray(ho.Z_corr)
        if z.shape != adata.obsm["X_pca"].shape:
            z = z.T
        embs[f"harmony_{key}"] = z
        print(f"  harmony {key} done", flush=True)

    import scvi
    for key in KEYS:
        b = ad.AnnData(adata.layers["counts"].copy(), obs=adata.obs.copy(), var=adata.var.copy())
        scvi.model.SCVI.setup_anndata(b, batch_key=key, layer=None)
        model = scvi.model.SCVI(b, n_layers=2, n_latent=N_LATENT, gene_likelihood="nb")
        model.train(accelerator="cpu", devices=1)
        embs[f"scvi_{key}"] = model.get_latent_representation()
        print(f"  scvi {key} done", flush=True)

    is_p4 = (adata.obs["library"] == "pool4").to_numpy()
    rows = []
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.model_selection import cross_val_score
    for name, E in embs.items():
        adata.obsm[f"X_{name}"] = E
        r = {"run": name}
        lis = pd.DataFrame(compute_lisi(np.ascontiguousarray(E),
                                        adata.obs[["library", "donor_id", "stage"]],
                                        ["library", "donor_id", "stage"]),
                           columns=["library", "donor_id", "stage"])
        r["mix_library"] = round((lis["library"].mean() - 1) / 3, 3)
        r["mix_donor"] = round((lis["donor_id"].mean() - 1) / 11, 3)
        r["cLISI_stage"] = round(float(lis["stage"].mean()), 3)
        y = adata.obs["stage"].to_numpy()
        r["stage_acc"] = round(float(cross_val_score(
            KNeighborsClassifier(30), E, y, cv=3).mean()), 3)
        # within pool4: knowing the pool cannot help, so this isolates real late-stage structure
        r["stage_acc_pool4"] = round(float(cross_val_score(
            KNeighborsClassifier(30), E[is_p4], y[is_p4], cv=3).mean()), 3)
        # does pool drive the clustering? compare against stage on the same clusters
        t = ad.AnnData(np.zeros((adata.n_obs, 1), dtype="float32"), obs=adata.obs.copy())
        t.obsm["X_e"] = np.ascontiguousarray(E)
        sc.pp.neighbors(t, use_rep="X_e", random_state=RS)
        sc.tl.leiden(t, resolution=1.0, key_added="cl", flavor="igraph", n_iterations=2,
                     directed=False)
        sc.tl.umap(t, random_state=RS)
        adata.obsm[f"X_umap_{name}"] = t.obsm["X_umap"]
        adata.obs[f"leiden_{name}"] = t.obs["cl"].to_numpy()
        r["n_clusters"] = int(t.obs["cl"].nunique())
        r["NMI_cluster_pool"] = round(nmi(t.obs["cl"], adata.obs["library"]), 3)
        r["NMI_cluster_stage"] = round(nmi(t.obs["cl"], adata.obs["stage"]), 3)
        rows.append(r); print(f"  scored {name}: {r}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(out / "olaniru_batch_metrics.csv", index=False)
    adata.write_h5ad(out / "olaniru_batch_compare.h5ad")

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = ["library", "stage", "donor_id"]
    fig, axes = plt.subplots(len(embs), len(cols), figsize=(4.2*len(cols), 3.8*len(embs)),
                             squeeze=False)
    for i, name in enumerate(embs):
        U = adata.obsm[f"X_umap_{name}"]
        for j, c in enumerate(cols):
            ax = axes[i][j]
            cats = adata.obs[c].astype(str); uq = sorted(cats.unique())
            cmap = plt.get_cmap("tab20" if len(uq) > 10 else "tab10")
            for k, u in enumerate(uq):
                s = (cats == u).to_numpy()
                ax.scatter(U[s, 0], U[s, 1], s=1.0, color=cmap(k % 20), label=u,
                           linewidths=0, rasterized=True)
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(c, fontsize=11)
                if len(uq) <= 12:
                    ax.legend(fontsize=6, markerscale=5, frameon=False, loc="best")
            if j == 0: ax.set_ylabel(name, fontsize=9)
    fig.suptitle(f"Olaniru 2023 within-study batch correction ({adata.n_obs:,} cells)",
                 fontsize=13, y=1.002)
    fig.savefig(out / "olaniru_batch_umaps.png", dpi=115, bbox_inches="tight")
    pd.set_option("display.width", 250)
    print("\n=== Olaniru within-study batch comparison ===")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
