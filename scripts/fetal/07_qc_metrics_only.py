#!/usr/bin/env python
"""QC metrics for the fetal datasets -- COMPUTED AND PLOTTED ONLY, NOTHING IS FILTERED.

The point is to look at the distributions and choose cutoffs deliberately, before any cell is
removed. Metrics and candidate thresholds are computed exactly as the lab pipeline's
scripts/1.1_rna_qc.py would (same mt/ribo gene patterns, same one-sided MAD rule, same
median+n*MAD MT threshold clamped to [floor, cap]), so whatever is chosen here transfers over
by setting n_mads / mt_pct_* in configs/qc_RNA.yaml.

Candidate thresholds are drawn at n_mads = 3, 4, 5 so the cost of each is visible. The paper's
own fixed cutoffs are drawn too, for reference:
  Olaniru 2023   nCount > 500, nFeature > 250, mitoRatio < 0.25
  Migliorini 2024 reports 73,274 nuclei passing QC + doublet exclusion over 6 donors.

Usage: python 07_qc_metrics_only.py --dataset olaniru|migliorini|both --outdir <dir>
"""
import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

N_MADS_GRID = [3, 4, 5]
MT_FLOOR, MT_CAP, MT_N_MADS = 3.0, 15.0, 5
PAPER = {"min_counts": 500, "min_genes": 250, "max_mt_pct": 25.0}


def mad(v):
    return float(np.median(np.abs(v - np.median(v))))


def low_outlier_cut(v, n_mads):
    """Lower bound of the one-sided MAD rule used by 1.1_rna_qc.py:is_mad_outlier."""
    m, d = float(np.median(v)), mad(v)
    return None if d == 0 else m - n_mads * d


def mt_cut(v, n_mads=MT_N_MADS, floor=MT_FLOOR, cap=MT_CAP):
    t = float(np.median(v)) + n_mads * mad(v)
    if floor is not None: t = max(t, floor)
    if cap is not None: t = min(t, cap)
    return t


def load(sample_id, h5):
    a = sc.read_10x_h5(h5)
    a.var_names_make_unique()
    if "feature_types" in a.var:
        a = a[:, a.var["feature_types"] == "Gene Expression"].copy()
    a.var["mt"] = a.var_names.str.startswith("MT-")
    a.var["ribo"] = a.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt", "ribo"], percent_top=None, log1p=True, inplace=True)
    a.obs["sample"] = sample_id
    return a


def per_sample_figure(obs, sample_id, kind, png):
    """4 metric panels + a counts-vs-genes scatter, with candidate cutoffs drawn, nothing applied."""
    fig = plt.figure(figsize=(21, 4.2))
    gs = fig.add_gridspec(1, 5, wspace=0.28)
    panels = [
        ("log1p_total_counts", "log1p(total counts)", True),
        ("log1p_n_genes_by_counts", "log1p(n genes)", True),
        ("pct_counts_mt", "% MT counts", False),
        ("pct_counts_ribo", "% ribo counts", False),
    ]
    colors = {3: "#d1495b", 4: "#edae49", 5: "#00798c"}
    for i, (col, label, is_mad_metric) in enumerate(panels):
        ax = fig.add_subplot(gs[0, i])
        v = obs[col].to_numpy()
        ax.hist(v, bins=80, color="#9aa5b1", edgecolor="none")
        if is_mad_metric:
            for nm in N_MADS_GRID:
                c = low_outlier_cut(v, nm)
                if c is None: continue
                n_lost = int((v < c).sum())
                ax.axvline(c, color=colors[nm], lw=1.6, ls="--")
                ax.text(c, ax.get_ylim()[1] * (0.97 - 0.09 * N_MADS_GRID.index(nm)), f" {nm}MAD n={n_lost}",
                        color=colors[nm], fontsize=8, va="top", ha="left")
            # the paper's fixed cutoff, in the same log1p space
            pv = np.log1p(PAPER["min_counts"] if col.endswith("total_counts") else PAPER["min_genes"])
            ax.axvline(pv, color="black", lw=1.4, ls=":")
            ax.text(pv, ax.get_ylim()[1] * 0.55, " paper", fontsize=8, va="top", ha="left")
        elif col == "pct_counts_mt":
            t = mt_cut(v)
            ax.axvline(t, color="#00798c", lw=1.6, ls="--")
            ax.text(t, ax.get_ylim()[1] * 0.97, f" MAD rule {t:.1f}%  n={int((v > t).sum())}",
                    color="#00798c", fontsize=8, va="top", ha="left")
            ax.axvline(PAPER["max_mt_pct"], color="black", lw=1.4, ls=":")
            ax.text(PAPER["max_mt_pct"], ax.get_ylim()[1] * 0.6, " paper 25%", fontsize=8, va="top", ha="left")
            ax.text(0.98, 0.5, f"median {np.median(v):.2f}%", transform=ax.transAxes,
                    fontsize=8, ha="right", color="#444")
        ax.set_xlabel(label, fontsize=9); ax.set_ylabel("cells" if i == 0 else "", fontsize=9)
        ax.tick_params(labelsize=8)
        for s in ("top", "right"): ax.spines[s].set_visible(False)

    ax = fig.add_subplot(gs[0, 4])
    sc_ = ax.scatter(obs["log1p_total_counts"], obs["log1p_n_genes_by_counts"],
                     c=obs["pct_counts_mt"], s=2, cmap="viridis",
                     vmin=0, vmax=max(5, float(np.percentile(obs["pct_counts_mt"], 99))))
    plt.colorbar(sc_, ax=ax, label="% MT", fraction=0.046)
    for nm in N_MADS_GRID:
        c = low_outlier_cut(obs["log1p_total_counts"].to_numpy(), nm)
        if c is not None: ax.axvline(c, color=colors[nm], lw=1.1, ls="--")
        c = low_outlier_cut(obs["log1p_n_genes_by_counts"].to_numpy(), nm)
        if c is not None: ax.axhline(c, color=colors[nm], lw=1.1, ls="--")
    ax.set_xlabel("log1p(total counts)", fontsize=9); ax.set_ylabel("log1p(n genes)", fontsize=9)
    ax.tick_params(labelsize=8)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(handles=[Line2D([], [], color=colors[n], ls="--", label=f"{n} MAD") for n in N_MADS_GRID]
                      + [Line2D([], [], color="black", ls=":", label="paper cutoff")],
              fontsize=7, loc="lower right", frameon=False)

    fig.suptitle(f"{sample_id}   ({kind}, {obs.shape[0]:,} cells called by Cell Ranger)   "
                 f"NOTHING FILTERED -- candidate cutoffs only", fontsize=11, y=1.04)
    fig.savefig(png, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, help="TSV: sample_id<TAB>kind<TAB>h5_path")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    out = Path(a.outdir); (out / "figures").mkdir(parents=True, exist_ok=True)

    rows, obs_all = [], []
    for line in Path(a.samples).read_text().strip().splitlines():
        sid, kind, h5 = line.split("\t")
        if not Path(h5).exists():
            print(f"  skip {sid}: {h5} not ready", file=sys.stderr); continue
        print(f"  {sid} ...", flush=True)
        ad = load(sid, h5)
        o = ad.obs.copy(); o["kind"] = kind
        per_sample_figure(o, sid, kind, out / "figures" / f"qc_{sid}.png")
        obs_all.append(o)

        r = {"sample": sid, "kind": kind, "cells_called": int(ad.n_obs),
             "median_umi": float(np.median(o["total_counts"])),
             "median_genes": float(np.median(o["n_genes_by_counts"])),
             "median_pct_mt": float(np.median(o["pct_counts_mt"])),
             "median_pct_ribo": float(np.median(o["pct_counts_ribo"]))}
        for nm in N_MADS_GRID:
            cc = low_outlier_cut(o["log1p_total_counts"].to_numpy(), nm)
            gc = low_outlier_cut(o["log1p_n_genes_by_counts"].to_numpy(), nm)
            fail = np.zeros(ad.n_obs, bool)
            if cc is not None: fail |= o["log1p_total_counts"].to_numpy() < cc
            if gc is not None: fail |= o["log1p_n_genes_by_counts"].to_numpy() < gc
            r[f"umi_cut_{nm}mad"] = round(float(np.expm1(cc)), 1) if cc is not None else None
            r[f"gene_cut_{nm}mad"] = round(float(np.expm1(gc)), 1) if gc is not None else None
            r[f"lost_{nm}mad"] = int(fail.sum())
            r[f"lost_{nm}mad_pct"] = round(100 * fail.sum() / ad.n_obs, 2)
        t = mt_cut(o["pct_counts_mt"].to_numpy())
        r["mt_cut_mad_rule"] = round(t, 2)
        r["lost_mt"] = int((o["pct_counts_mt"].to_numpy() > t).sum())
        pf = ((o["total_counts"] <= PAPER["min_counts"]) | (o["n_genes_by_counts"] <= PAPER["min_genes"])
              | (o["pct_counts_mt"] >= PAPER["max_mt_pct"])).to_numpy()
        r["lost_paper_rule"] = int(pf.sum())
        r["lost_paper_rule_pct"] = round(100 * pf.sum() / ad.n_obs, 2)
        rows.append(r)

    if not rows:
        sys.exit("no samples ready")
    df = pd.DataFrame(rows)
    df.to_csv(out / "qc_metrics_summary.csv", index=False)
    pd.concat(obs_all)[["sample", "kind", "total_counts", "n_genes_by_counts",
                        "pct_counts_mt", "pct_counts_ribo"]].to_csv(out / "qc_per_cell.csv.gz", index=False)
    print("\n" + df.to_string(index=False))
    print(f"\nwrote {out}/qc_metrics_summary.csv and {len(rows)} figures under {out}/figures/")


if __name__ == "__main__":
    main()
