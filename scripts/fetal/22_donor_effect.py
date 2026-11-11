#!/usr/bin/env python
"""How large is the donor effect in Olaniru, really -- and how much of it is actually the pool?

Donor is a major variance component in single-cell data generally, so "donor is confounded with
stage here" is not an answer to "how big is it". This measures it directly, using the one pairing
in the design that isolates it.

  PURE DONOR   pool1 carries TWO embryos at the same gestational age, 13PCW (n=2,137) and
               13PCW_1 (n=2,213): same GEM well, same ambient pool, same doublet rate, same stage.
               Any separation between them is biology -- genotype, sex, tissue handling of that
               embryo -- with no technical contribution at all.

  DONOR+POOL   the same stage sampled in two different pools (PCW12 pool1/pool2, PCW13 pool1/pool3,
               PCW14 pool2/pool3, PCW19 pool2/pool3). Here donor and pool vary together.

Comparing the two tells you whether the apparent donor effect is the embryo or the batch, which is
what decides whether donor is worth correcting for in a study where it is nested in stage.

Metric: a 30-NN classifier's accuracy at naming which donor a cell came from, with 0.5 = the two
are indistinguishable and 1.0 = fully separated; plus iLISI mixing on the same pair. Each pair is
processed on its own HVGs/PCA so one pair's variation cannot leak into another's.
"""
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
from harmonypy import compute_lisi

BASE = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas")
OUT = BASE / "olaniru2023" / "integration"
RS, N_HVG, N_PCS = 1234, 2000, 30

PAIRS = [
    ("PURE DONOR  (same pool, same stage)", "PCW13", [("pool1", "13PCW"), ("pool1", "13PCW_1")]),
    ("donor+pool  (same stage, 2 pools)",   "PCW12", [("pool1", "12PCW"), ("pool2", "12PCW")]),
    ("donor+pool  (same stage, 2 pools)",   "PCW13", [("pool1", "13PCW"), ("pool3", "13PCW")]),
    ("donor+pool  (same stage, 2 pools)",   "PCW14", [("pool2", "14PCW"), ("pool3", "14PCW")]),
    ("donor+pool  (same stage, 2 pools)",   "PCW19", [("pool2", "19PCW"), ("pool3", "19PCW")]),
]


def load():
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
        for c in ["library", "donor_label", "donor_id", "stage"]:
            o[c] = m.loc[a.obs_names, c].to_numpy()
        a.obs = o
        a.obs_names = [f"{lib}|{b}" for b in a.obs_names]
        parts.append(a)
    return ad.concat(parts, join="inner", merge="unique")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    full = load()
    print(f"loaded {full.n_obs:,} Olaniru cells\n", flush=True)
    rows = []
    for kind, stage, members in PAIRS:
        sel = np.zeros(full.n_obs, bool)
        lab = np.empty(full.n_obs, dtype=object)
        for lib, dl in members:
            s = ((full.obs["library"] == lib) & (full.obs["donor_label"] == dl)).to_numpy()
            sel |= s
            lab[s] = f"{lib}/{dl}"
        sub = full[sel].copy()
        sub.obs["pair_label"] = lab[sel]
        # each pair on its own HVGs/PCA so other donors' variation cannot leak in
        sub.layers["counts"] = sub.X.copy()
        sc.pp.highly_variable_genes(sub, flavor="seurat_v3", n_top_genes=N_HVG, subset=True)
        sc.pp.normalize_total(sub, target_sum=1e4); sc.pp.log1p(sub)
        sc.tl.pca(sub, n_comps=N_PCS, random_state=RS)
        E = sub.obsm["X_pca"]; y = sub.obs["pair_label"].to_numpy()
        acc = float(cross_val_score(KNeighborsClassifier(30), E, y, cv=3).mean())
        maj = float(pd.Series(y).value_counts(normalize=True).max())
        lis = compute_lisi(np.ascontiguousarray(E), sub.obs[["pair_label"]], ["pair_label"])
        mix = float((lis.mean() - 1) / 1)
        rows.append({"comparison": kind, "stage": stage,
                     "donors": " vs ".join(f"{l}/{d}" for l, d in members),
                     "n_cells": int(sub.n_obs), "donor_acc": round(acc, 3),
                     "majority": round(maj, 3),
                     "separability": round((acc - maj) / (1 - maj), 3),
                     "mix_iLISI": round(mix, 3)})
        print(f"  {kind} {stage}: acc={acc:.3f} (majority {maj:.3f}) mix={mix:.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "donor_effect.csv", index=False)
    pd.set_option("display.width", 250)
    print("\n=== donor effect ===")
    print(df.to_string(index=False))
    pure = df[df.comparison.str.startswith("PURE")]["separability"].mean()
    both = df[df.comparison.str.startswith("donor+pool")]["separability"].mean()
    print(f"\n  pure donor separability      : {pure:.3f}")
    print(f"  donor+pool separability (mean): {both:.3f}")
    print(f"  attributable to pool          : {both - pure:+.3f}")
    print("\n  separability 0 = donors indistinguishable, 1 = fully separated")


if __name__ == "__main__":
    main()
