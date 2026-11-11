#!/usr/bin/env python
"""Merge the 10 fetal libraries and run Harmony three ways, to see what each batch key costs.

Normalisation follows the lab pipeline's "log" path exactly (configs/normalize_rna.yaml +
5.3_normalize_merged.py): normalize_total -> log1p -> 2000 HVGs (flavor=seurat) -> 30 PCs,
no scaling and no covariate regression, random_state 1234. Harmony via harmonypy's
run_harmony, as in 5.4_harmony_integrate.py.

Three batch keys, plus the uncorrected baseline:
  study      2 levels   the only technical axis CROSSED with stage (PCW14 and PCW18 appear in
                        both studies), so it is the one key that can be estimated separately
                        from developmental biology.
  donor_id  18 levels   every donor sits at exactly ONE stage (0 of 18 span two), so donor is
                        perfectly NESTED in stage.
  library   10 levels   Olaniru pools span 2-3 stages each (safe), but every Migliorini library
                        is a single stage (nested).

Because of that nesting, "does it mix the batches" is not a sufficient question -- a key that
nests stage can score perfectly on batch mixing by deleting the biology. So every run is scored
on BOTH axes: batch mixing (iLISI, higher is better) and stage preservation (cLISI on stage,
lower is better; kNN stage purity; and a kNN classifier's accuracy at recovering a cell's stage
from the embedding, which is the most direct test of whether stage survived).
"""
import argparse, json
from pathlib import Path

import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from harmonypy import run_harmony, compute_lisi
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import silhouette_score

BASE = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas")
RS, N_HVG, N_PCS = 1234, 2000, 30
KEYS = ["study", "donor_id", "library"]


def merge():
    meta = pd.read_csv(BASE / "qc_review" / "fetal_cell_metadata.csv.gz")
    meta = meta[meta["in_atlas"]].copy()
    parts = []
    for ds in ("olaniru2023", "migliorini2024"):
        h5d = BASE / ds / "qc_filtered" / "detect_doublets" / "h5ad"
        for p in sorted(h5d.glob("*.doublets_removed.h5ad")):
            if "autothresh" in p.name:
                continue
            lib = p.name.split(".")[0]
            a = sc.read_h5ad(p)
            # the lab's 4.1/5.3 start from SoupX-corrected counts when present
            if "soupx_counts" in a.layers:
                a.X = a.layers["soupx_counts"].copy()
            a.var_names_make_unique()
            m = meta[meta["library"] == lib].set_index("barcode")
            keep = [b for b in a.obs_names if b in m.index]
            a = a[keep].copy()
            for c in ["study", "study_gse", "assay", "material", "multiplexed", "hto_available",
                      "donor_label", "donor_id", "pcw", "stage", "stage_bin", "library"]:
                a.obs[c] = m.loc[a.obs_names, c].to_numpy()
            a.obs_names = [f"{lib}|{b}" for b in a.obs_names]   # ~800 barcodes collide across libraries
            parts.append(a)
            print(f"  loaded {lib}: {a.n_obs:,} x {a.n_vars:,}", flush=True)
    merged = ad.concat(parts, join="inner", merge="unique")
    print(f"\nmerged: {merged.n_obs:,} cells x {merged.n_vars:,} genes (intersection)")
    assert merged.obs_names.is_unique
    return merged


def metrics(emb, obs, label="stage", k=30, sub=20000, rs=RS):
    rng = np.random.default_rng(rs)
    idx = rng.choice(len(obs), min(sub, len(obs)), replace=False)
    E, O = emb[idx], obs.iloc[idx]
    out = {}
    lis = compute_lisi(E, O[KEYS + [label]], KEYS + [label], perplexity=min(30, len(idx)//4))
    lis = pd.DataFrame(lis, columns=KEYS + [label])
    for key in KEYS:
        n = O[key].nunique()
        out[f"iLISI_{key}"] = round(float(lis[key].mean()), 3)
        out[f"iLISI_{key}_norm"] = round(float((lis[key].mean() - 1) / max(n - 1, 1)), 3)
    out[f"cLISI_{label}"] = round(float(lis[label].mean()), 3)
    y = O[label].to_numpy()
    knn = KNeighborsClassifier(n_neighbors=k)
    out[f"{label}_knn_acc"] = round(float(cross_val_score(knn, E, y, cv=3).mean()), 3)
    out[f"{label}_baseline_acc"] = round(float(pd.Series(y).value_counts(normalize=True).max()), 3)
    out[f"{label}_silhouette"] = round(float(silhouette_score(E, y, sample_size=8000,
                                                              random_state=rs)), 3)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(BASE / "integration"))
    ap.add_argument("--hvg-batch-key", default="study",
                    help="select HVGs WITHIN this batch then rank by how many batches a gene is "
                         "variable in. Default 'study'. Pass '' for the pooled selection the "
                         "lab's 5.3 uses -- pooled is wrong for this merge: pooled variance is "
                         "dominated by cells-vs-nuclei (MT-*, RPS/RPL, MALAT1 differ 20-300x "
                         "between materials), so it builds the study effect into the PCs Harmony "
                         "is then asked to remove.")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    mp = out / "fetal_merged.h5ad"
    adata = sc.read_h5ad(mp) if mp.exists() else merge()
    if not mp.exists():
        adata.write_h5ad(mp)

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    bk = a.hvg_batch_key or None
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat", batch_key=bk)
    hv = adata[:, adata.var["highly_variable"]].copy()
    sc.tl.pca(hv, n_comps=N_PCS, random_state=RS)
    adata.obsm["X_pca"] = hv.obsm["X_pca"]
    print(f"\nPCA done: {N_HVG} HVGs (batch_key={bk}) -> {N_PCS} PCs", flush=True)
    # how much of the HVG set is driven by the assay difference rather than biology
    for fam, pat in [("MT-", r"^MT-"), ("ribosomal", r"^RP[SL]"), ("hemoglobin", r"^HB[ABGDMZE]")]:
        n = int(hv.var_names.str.match(pat).sum())
        print(f"    {fam:11} genes among HVGs: {n}")
    if bk:
        nb = adata.var.loc[adata.var["highly_variable"], "highly_variable_nbatches"]
        print(f"    HVGs variable in both studies: {int((nb == 2).sum())} / {len(nb)}")

    embs = {"uncorrected": adata.obsm["X_pca"]}
    for key in KEYS:
        print(f"\nharmony on '{key}' ({adata.obs[key].nunique()} levels) ...", flush=True)
        ho = run_harmony(adata.obsm["X_pca"], adata.obs, [key], random_state=RS)
        z = np.asarray(ho.Z_corr)
        if z.shape != adata.obsm["X_pca"].shape:
            z = z.T
        assert z.shape == adata.obsm["X_pca"].shape, f"{z.shape}"
        embs[f"harmony_{key}"] = z
        adata.obsm[f"X_pca_harmony_{key}"] = z

    rows = []
    for name, E in embs.items():
        print(f"\nscoring {name} ...", flush=True)
        r = {"run": name}; r.update(metrics(E, adata.obs))
        rows.append(r); print("   ", r, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "integration_metrics.csv", index=False)

    # UMAPs: rows = runs, cols = colourings
    colorings = ["study", "stage", "donor_id", "library", "material"]
    fig, axes = plt.subplots(len(embs), len(colorings),
                             figsize=(4.0*len(colorings), 3.7*len(embs)), squeeze=False)
    for i, (name, E) in enumerate(embs.items()):
        tmp = ad.AnnData(np.zeros((adata.n_obs, 1)), obs=adata.obs.copy())
        tmp.obsm["X_pca"] = E
        sc.pp.neighbors(tmp, use_rep="X_pca", random_state=RS)
        sc.tl.umap(tmp, random_state=RS)
        adata.obsm[f"X_umap_{name}"] = tmp.obsm["X_umap"]
        for j, c in enumerate(colorings):
            axp = axes[i][j]
            cats = adata.obs[c].astype(str)
            uq = sorted(cats.unique())
            cmap = plt.get_cmap("tab20" if len(uq) > 10 else "tab10")
            for n_, u in enumerate(uq):
                s = (cats == u).to_numpy()
                axp.scatter(tmp.obsm["X_umap"][s, 0], tmp.obsm["X_umap"][s, 1], s=0.8,
                            color=cmap(n_ % 20), label=u, linewidths=0, rasterized=True)
            axp.set_xticks([]); axp.set_yticks([])
            if i == 0: axp.set_title(c, fontsize=11)
            if j == 0: axp.set_ylabel(name, fontsize=10)
            if len(uq) <= 10 and i == 0:
                axp.legend(fontsize=5.5, markerscale=5, frameon=False, loc="best")
        print(f"  umap {name} done", flush=True)
    fig.suptitle("Fetal pancreas atlas -- Harmony batch key comparison "
                 f"({adata.n_obs:,} cells/nuclei)", fontsize=13, y=1.002)
    fig.savefig(out / "integration_umaps.png", dpi=115, bbox_inches="tight")

    adata.write_h5ad(out / "fetal_integrated.h5ad")
    pd.set_option("display.width", 250)
    print("\n=== metrics ===")
    print(df.to_string(index=False))
    print(f"\nwrote {out}/integration_metrics.csv, integration_umaps.png, fetal_integrated.h5ad")


if __name__ == "__main__":
    main()
