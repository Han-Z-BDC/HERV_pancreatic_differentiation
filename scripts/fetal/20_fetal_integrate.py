#!/usr/bin/env python
"""Fetal pancreas atlas: Harmony AND scVI, each with three batch keys, scored the same way.

Run in three stages so the GPU job stays small and each piece is restartable:
  --stage prep      merge the 10 libraries, pick HVGs, normalise, PCA, run Harmony x3   (CPU)
  --stage scvi      scVI x3 on raw counts of the same HVG set                           (GPU)
  --stage score     metrics + UMAPs for all 7 embeddings                                (CPU)

Decisions and why
-----------------
target_sum=1e4 (NOT scanpy's default None). With target_sum=None, normalize_total gives every
cell "a total count equal to the median of total counts for observations before normalization"
-- a POOLED statistic. 76% of this atlas is nuclei, which are deeper, so the pooled median is
pulled toward them and the Olaniru cells get a scale factor they would not get on their own.
Fixing the target makes per-library and merged normalisation identical and makes the result
independent of which samples happen to be included. It also matches the 10000 the lab's own
remove_background.yaml uses.

ONE shared HVG set for both methods: seurat_v3 on raw counts, n_top_genes=2000,
batch_key="library", subset=True -- the same call 5.5_scvi_integrate.py makes (its
"sample_id" is our "library"), so the scVI arm is unchanged from the validated pipeline. Batch-aware because sc-best-practices is explicit that "it is crucial
that gene selection is performed in a batch-aware manner", and pooled variance here is
dominated by the assay: %MT differs ~80x between cells and nuclei, %ribo ~25x. Sharing one gene
set across both methods is a deliberate deviation from 5.3's log path (which uses flavor=seurat
on log data with no batch_key) -- otherwise method and gene set would be confounded and the
comparison would not isolate either. batch_key="study" rather than "library" because six of the
ten libraries are single-stage, so demanding variability across all ten biases against
stage-specific genes.

scVI gets RAW COUNTS (it models counts with its own library-size term), so the normalisation
choice above affects only Harmony's input. That makes scVI the cleaner of the two comparisons.

Batch keys, and the trap in them:
  study     2 levels   the only technical axis CROSSED with stage (PCW14, PCW18 in both studies)
  donor_id 18 levels   every donor sits at exactly ONE stage -- 0 of 18 span two
  library  10 levels   Olaniru pools span 2-3 stages; every Migliorini library spans one
Because donor and library nest stage, batch mixing alone cannot score these: a key that nests
stage can score perfectly by deleting the biology. Every embedding is therefore scored on BOTH
mixing (iLISI) and stage retention (cLISI, silhouette, and kNN stage-recovery accuracy against
the majority-class baseline). stage/pcw is NEVER used as a batch key or covariate.
"""
import argparse, json
from pathlib import Path

import numpy as np, pandas as pd, scanpy as sc, anndata as ad

BASE = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas")
OUT = BASE / "integration"
RS, N_HVG, N_LATENT, N_PCS, TARGET_SUM = 1234, 2000, 30, 30, 1e4
SUFFIX = ""
HVG_BATCH_KEY = "library"   # 5.5_scvi_integrate.py's "sample_id"
KEYS = ["study", "donor_id", "library"]
META_COLS = ["study", "study_gse", "assay", "material", "multiplexed", "hto_available",
             "donor_label", "donor_id", "pcw", "stage", "stage_bin", "library"]


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
            if "soupx_counts" in a.layers:
                a.X = a.layers["soupx_counts"].copy()
            a.layers.clear(); a.obsm.clear(); a.varm.clear(); a.uns.clear()
            a.var_names_make_unique()
            m = meta[meta["library"] == lib].set_index("barcode")
            a = a[[b for b in a.obs_names if b in m.index]].copy()
            keep_obs = pd.DataFrame(index=a.obs_names)
            for c in META_COLS:
                keep_obs[c] = m.loc[a.obs_names, c].to_numpy()
            a.obs = keep_obs
            a.obs_names = [f"{lib}|{b}" for b in a.obs_names]   # ~800 barcodes collide
            parts.append(a)
            print(f"  {lib}: {a.n_obs:,} x {a.n_vars:,}", flush=True)
    merged = ad.concat(parts, join="inner", merge="unique")
    assert merged.obs_names.is_unique
    print(f"\nmerged: {merged.n_obs:,} cells x {merged.n_vars:,} genes (intersection)")
    return merged


def prep():
    OUT.mkdir(parents=True, exist_ok=True)
    adata = merge()
    adata.layers["counts"] = adata.X.copy()

    # HVG selection matches 5.5_scvi_integrate.py exactly: seurat_v3 on raw counts with
    # batch_key at the library level (5.5's "sample_id"). ONE set shared by Harmony and scVI,
    # so method and gene set are not confounded -- this is the change that makes the Harmony
    # arm batch-aware, since 5.3's log path selects pooled with no batch_key.
    #
    # Both candidate keys are computed so the choice is visible rather than assumed. The
    # tradeoff: per-library selection ranks genes by how many of the 10 libraries they vary in,
    # and six libraries are single-stage, so a gene expressed only late is non-variable in the
    # libraries lacking that stage and gets deprioritised. Per-study (2 batches) is gentler on
    # stage-specific genes. Reported, not silently chosen.
    for cand in ("library", "study"):
        # share X by reference rather than adata.copy() -- seurat_v3 HVG only reads the matrix,
        # and a full copy would duplicate the counts layer too (~3.4 GB per iteration)
        tmp = ad.AnnData(adata.X, obs=adata.obs[[cand]].copy(),
                         var=pd.DataFrame(index=adata.var_names))
        sc.pp.highly_variable_genes(tmp, flavor="seurat_v3", n_top_genes=N_HVG,
                                    batch_key=cand, subset=False)
        globals()[f"_hvg_{cand}"] = set(tmp.var_names[tmp.var["highly_variable"]])
        del tmp
    a_, b_ = globals()["_hvg_library"], globals()["_hvg_study"]
    print(f"\nHVG set overlap: batch_key=library vs study -> {len(a_ & b_)} shared "
          f"of {N_HVG} ({100*len(a_ & b_)/N_HVG:.0f}%)")

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=N_HVG,
                                batch_key=HVG_BATCH_KEY, subset=True)
    print(f"using batch_key={HVG_BATCH_KEY} (matches 5.5): {adata.n_vars} HVGs")
    nb = adata.var.get("highly_variable_nbatches")
    if nb is not None:
        nlev = adata.obs[HVG_BATCH_KEY].nunique()
        print(f"  variable in all {nlev} batches: {int((nb == nlev).sum())} / {len(nb)}"
              f"   median nbatches {float(nb.median()):.1f}")
    for fam, pat in [("MT-", r"^MT-"), ("ribosomal", r"^RP[SL]"), ("hemoglobin", r"^HB[ABGDMZE]")]:
        print(f"  {fam:11} among HVGs: {int(adata.var_names.str.match(pat).sum())}")

    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, n_comps=N_PCS, random_state=RS)
    print(f"PCA: {N_PCS} PCs on target_sum={TARGET_SUM:g} log1p data", flush=True)

    from harmonypy import run_harmony
    for key in KEYS:
        n = adata.obs[key].nunique()
        print(f"\nharmony '{key}' ({n} levels) ...", flush=True)
        ho = run_harmony(adata.obsm["X_pca"], adata.obs, [key], random_state=RS)
        z = np.asarray(ho.Z_corr)
        if z.shape != adata.obsm["X_pca"].shape:
            z = z.T
        assert z.shape == adata.obsm["X_pca"].shape, f"harmonypy Z_corr {z.shape}"
        adata.obsm[f"X_harmony_{key}"] = z
    adata.write_h5ad(OUT / "fetal_hvg_harmony.h5ad")
    print(f"\nwrote {OUT}/fetal_hvg_harmony.h5ad")


def scvi_stage(max_epochs, accelerator="gpu"):
    import scvi, torch
    print(f"scvi-tools {scvi.__version__}  torch {torch.__version__}  "
          f"cuda={torch.cuda.is_available()}", flush=True)
    p = OUT / "fetal_hvg_harmony.h5ad"
    adata = sc.read_h5ad(p)
    for key in KEYS:
        print(f"\nscVI '{key}' ({adata.obs[key].nunique()} levels) ...", flush=True)
        a = ad.AnnData(adata.layers["counts"].copy(), obs=adata.obs.copy(),
                       var=adata.var.copy())
        scvi.model.SCVI.setup_anndata(a, batch_key=key, layer=None)
        model = scvi.model.SCVI(a, n_layers=2, n_latent=N_LATENT, gene_likelihood="nb")
        model.train(max_epochs=max_epochs, accelerator=accelerator, devices=1)
        adata.obsm[f"X_scvi_{key}"] = model.get_latent_representation()
        print(f"  latent {adata.obsm[f'X_scvi_{key}'].shape}", flush=True)
    dest = OUT / f"fetal_integrated{SUFFIX}.h5ad"
    adata.write_h5ad(dest)
    print(f"\nwrote {dest}")


def metrics(emb, obs, label="stage", k=30, sub=20000):
    from harmonypy import compute_lisi
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.model_selection import cross_val_score
    from sklearn.metrics import silhouette_score
    rng = np.random.default_rng(RS)
    idx = rng.choice(len(obs), min(sub, len(obs)), replace=False)
    E, O = np.ascontiguousarray(emb[idx]), obs.iloc[idx]
    out = {}
    lis = pd.DataFrame(compute_lisi(E, O[KEYS + [label]], KEYS + [label]),
                       columns=KEYS + [label])
    for key in KEYS:
        n = O[key].nunique()
        out[f"iLISI_{key}"] = round(float(lis[key].mean()), 3)
        out[f"mix_{key}"] = round(float((lis[key].mean() - 1) / max(n - 1, 1)), 3)
    out[f"cLISI_{label}"] = round(float(lis[label].mean()), 3)
    y = O[label].to_numpy()
    out[f"{label}_knn_acc"] = round(float(cross_val_score(
        KNeighborsClassifier(n_neighbors=k), E, y, cv=3).mean()), 3)
    out["stage_baseline"] = round(float(pd.Series(y).value_counts(normalize=True).max()), 3)
    out[f"{label}_sil"] = round(float(silhouette_score(E, y, sample_size=8000,
                                                       random_state=RS)), 3)
    return out


def score():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cands = [OUT / "fetal_integrated.h5ad", OUT / "fetal_integrated_cpu.h5ad"]
    p = next((c for c in cands if c.exists()), None)
    if p is None:
        p = OUT / "fetal_hvg_harmony.h5ad"
        print(f"NOTE: no scVI output yet -- scoring the Harmony arm only ({p.name}). "
              f"Re-run this stage once scVI lands to score all 7 embeddings.")
    else:
        print(f"scoring {p.name}")
    adata = sc.read_h5ad(p)
    embs = {"uncorrected": adata.obsm["X_pca"]}
    for m in ("harmony", "scvi"):
        for key in KEYS:
            k = f"X_{m}_{key}"
            if k in adata.obsm:
                embs[f"{m}_{key}"] = adata.obsm[k]

    rows = []
    for name, E in embs.items():
        r = {"run": name}; r.update(metrics(E, adata.obs)); rows.append(r)
        print(f"  {name}: {r}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"integration_metrics{SUFFIX}.csv", index=False)

    colorings = ["study", "stage", "donor_id", "library", "material"]
    fig, axes = plt.subplots(len(embs), len(colorings),
                             figsize=(3.9*len(colorings), 3.6*len(embs)), squeeze=False)
    for i, (name, E) in enumerate(embs.items()):
        t = ad.AnnData(np.zeros((adata.n_obs, 1), dtype="float32"), obs=adata.obs.copy())
        t.obsm["X_e"] = np.ascontiguousarray(E)
        sc.pp.neighbors(t, use_rep="X_e", random_state=RS)
        sc.tl.umap(t, random_state=RS)
        adata.obsm[f"X_umap_{name}"] = t.obsm["X_umap"]
        for j, c in enumerate(colorings):
            axp = axes[i][j]
            cats = adata.obs[c].astype(str)
            uq = sorted(cats.unique())
            cmap = plt.get_cmap("tab20" if len(uq) > 10 else "tab10")
            for n_, u in enumerate(uq):
                s = (cats == u).to_numpy()
                axp.scatter(t.obsm["X_umap"][s, 0], t.obsm["X_umap"][s, 1], s=0.7,
                            color=cmap(n_ % 20), label=u, linewidths=0, rasterized=True)
            axp.set_xticks([]); axp.set_yticks([])
            if i == 0:
                axp.set_title(c, fontsize=11)
                if len(uq) <= 10:
                    axp.legend(fontsize=5.5, markerscale=5, frameon=False, loc="best")
            if j == 0: axp.set_ylabel(name, fontsize=9)
        print(f"  umap {name} done", flush=True)
    fig.suptitle(f"Fetal pancreas atlas: Harmony vs scVI x 3 batch keys ({adata.n_obs:,} cells)",
                 fontsize=13, y=1.001)
    fig.savefig(OUT / f"integration_umaps{SUFFIX}.png", dpi=110, bbox_inches="tight")
    adata.write_h5ad(OUT / f"fetal_integrated_scored{SUFFIX}.h5ad")
    pd.set_option("display.width", 260)
    print("\n=== metrics (mix_* higher=better mixed; stage_knn_acc higher=stage preserved) ===")
    print(df.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["prep", "scvi", "score"])
    ap.add_argument("--max-epochs", type=int, default=None)
    ap.add_argument("--accelerator", default="gpu", choices=["gpu", "cpu"],
                    help="cpu is a real option here: the ah200 queue is ~10h deep while acpu "
                         "starts in minutes, and 81k cells x 2000 genes at n_latent=30 is small")
    ap.add_argument("--out-suffix", default="",
                    help="so a CPU and a GPU run can race without overwriting each other")
    a = ap.parse_args()
    SUFFIX = a.out_suffix
    {"prep": prep,
     "scvi": lambda: scvi_stage(a.max_epochs, a.accelerator),
     "score": score}[a.stage]()
