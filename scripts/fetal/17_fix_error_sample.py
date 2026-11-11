#!/usr/bin/env python
"""Repair ONE sample whose Scrublet auto-threshold misfired, by matching its doublet rate to
the average of its sibling samples.

Scrublet's automatic threshold is a KDE minimum in the simulated-doublet score histogram, and
it is unreliable when that histogram is not cleanly bimodal. On these libraries it failed in
both directions: it called 5-30 doublets on nine samples (so the pipeline's "union of two
tools" added 2 cells in 21,387 and the stated safety margin did not exist), and on S6_18PCW it
called 3,582, removing 1,983 nuclei on Scrublet's word alone -- 29.6% of the sample against
12.3-16.2% in its five siblings.

This fixes only the broken sample, and only its call COUNT. The target is the mean of the
sibling samples' own doublet rates. The cells to call are then taken from the top of a
combined ranking of both tools' scores: Scrublet's SCORE is informative here (AUROC 0.856-0.871
against scDblFinder's calls on the Olaniru pools), it was only the cut point that failed, so
discarding the score entirely would throw away a real second opinion.

Nothing else is touched: the other samples keep the calls step 1.3 made, the shared pipeline
under single-cell-pipeline/ is unmodified, and the original h5ad files are preserved as
*_autothresh.h5ad before anything is rewritten.
"""
import argparse, shutil
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

BASE = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample", required=True, help="the sample to repair")
    ap.add_argument("--match", choices=["scrublet", "union"], default="scrublet",
                    help="scrublet: normalise ONLY Scrublet's own call rate to the sibling mean, "
                         "leaving scDblFinder's calls untouched, then re-form the union (the "
                         "default -- it repairs the component that actually broke). "
                         "union: override the sample's total union rate instead.")
    ap.add_argument("--rank", choices=["mean", "scdbl", "scrublet"], default="scrublet",
                    help="which ranking picks the cells for the replaced call set")
    ap.add_argument("--rate", choices=["mean_of_rates", "pooled"], default="mean_of_rates")
    ap.add_argument("--apply", action="store_true", help="without this, report only")
    a = ap.parse_args()

    h5d = BASE / a.dataset / "qc_filtered" / "detect_doublets" / "h5ad"
    paths = sorted(h5d.glob("*.doublets_annotated.h5ad"))
    sibs, target_path = [], None
    for p in paths:
        sid = p.name.split(".")[0]
        if sid == a.sample:
            target_path = p
            continue
        adata = ad.read_h5ad(p, backed="r")
        un = adata.obs["doublet_filter_union"].to_numpy().astype(bool)
        sc_a = adata.obs["scrublet_predicted_doublet"].to_numpy().astype(bool)
        sd_a = adata.obs["scDblFinder_predicted_doublet"].to_numpy().astype(bool)
        sibs.append({"sample": sid, "cells": adata.n_obs,
                     "scDbl": int(sd_a.sum()), "scrublet": int(sc_a.sum()),
                     "union": int(un.sum()),
                     "scrublet_rate": float(sc_a.mean()), "rate": float(un.mean())})
    if target_path is None:
        raise SystemExit(f"{a.sample} not found under {h5d}")

    sdf = pd.DataFrame(sibs)
    col = "scrublet_rate" if a.match == "scrublet" else "rate"
    cnt = "scrublet" if a.match == "scrublet" else "union"
    rate_mean = float(sdf[col].mean())
    rate_pooled = float(sdf[cnt].sum() / sdf["cells"].sum())
    rate = rate_mean if a.rate == "mean_of_rates" else rate_pooled

    print("=== sibling samples (kept exactly as step 1.3 called them) ===")
    print(sdf.assign(scrublet_pct=(100*sdf["scrublet_rate"]).round(4),
                     union_pct=(100*sdf["rate"]).round(2)).to_string(index=False))
    print(f"\n  matching on             : {a.match}")
    print(f"  mean of per-sample rates: {100*rate_mean:.4f}%")
    print(f"  pooled rate             : {100*rate_pooled:.4f}%")
    print(f"  using                   : {100*rate:.4f}%  ({a.rate})")

    adata = ad.read_h5ad(target_path)
    o = adata.obs
    n = adata.n_obs
    n_target = int(round(rate * n))

    sd_call = o["scDblFinder_predicted_doublet"].to_numpy().astype(bool)
    sc_auto = o["scrublet_predicted_doublet"].to_numpy().astype(bool)
    sd_s = o["scDblFinder_doublet_score"].to_numpy()
    sc_s = o["scrublet_doublet_score"].to_numpy()

    if a.rank == "scdbl":
        score = pd.Series(sd_s).rank(pct=True).to_numpy()
    elif a.rank == "scrublet":
        score = pd.Series(sc_s).rank(pct=True).to_numpy()
    else:
        score = (pd.Series(sd_s).rank(pct=True).to_numpy()
                 + pd.Series(sc_s).rank(pct=True).to_numpy()) / 2

    # ties: keep everything at or above the n_target-th score rather than splitting a tie by index
    thr = np.partition(score, -n_target)[-n_target]
    picked = score >= thr

    old_union = sd_call | sc_auto
    if a.match == "scrublet":
        # replace ONLY Scrublet's call set, then re-form the union exactly as 1.3 does.
        # scDblFinder's calls are untouched -- on this sample its 14.2% already sits at the
        # sibling mean of 14.03%, so it was never the component that failed.
        new_scrub = picked
        new = sd_call | new_scrub
    else:
        new_scrub = sc_auto
        new = picked

    print(f"\n=== {a.sample} ===")
    print(f"  cells                          {n:,}")
    print(f"  scDblFinder (UNCHANGED)        {int(sd_call.sum()):,}  ({100*sd_call.mean():.2f}%)")
    print(f"  Scrublet, auto-threshold       {int(sc_auto.sum()):,}  ({100*sc_auto.mean():.2f}%)  <-- broken")
    print(f"  Scrublet, sibling-matched      {int(new_scrub.sum()):,}  ({100*new_scrub.mean():.3f}%)"
          f"   target {n_target:,}")
    print(f"  union, as 1.3 ran it           {int(old_union.sum()):,}  ({100*old_union.mean():.2f}%)")
    print(f"  union, repaired                {int(new.sum()):,}  ({100*new.mean():.2f}%)")
    print(f"  Scrublet-only calls remaining  {int((new_scrub & ~sd_call).sum()):,}"
          f"   (was {int((sc_auto & ~sd_call).sum()):,})")
    print(f"  nuclei kept: {int((~old_union).sum()):,} -> {int((~new).sum()):,} "
          f"(+{int(old_union.sum()) - int(new.sum()):,})")

    if not a.apply:
        print("\nreport only -- pass --apply to rewrite")
        return

    removed = h5d / f"{a.sample}.doublets_removed.h5ad"
    for src, dst in [(removed, h5d / f"{a.sample}.doublets_removed_autothresh.h5ad"),
                     (target_path, h5d / f"{a.sample}.doublets_annotated_autothresh.h5ad")]:
        if src.exists() and not dst.exists():
            shutil.copy2(str(src), str(dst))

    adata.obs["doublet_filter_union_autothresh"] = old_union
    adata.obs["scrublet_predicted_doublet_autothresh"] = sc_auto
    adata.obs["scrublet_predicted_doublet"] = new_scrub
    adata.obs["doublet_filter_union"] = new
    adata.obs["doublet_filter_consensus"] = sd_call & new_scrub
    adata.obs["doublet_rule"] = (f"Scrublet rate matched to sibling mean "
                                 f"({100*rate:.4f}%, {a.rank} rank); scDblFinder unchanged")
    adata.write_h5ad(target_path)
    adata[~new].copy().write_h5ad(removed)
    print(f"\nrewrote {removed.name}: {int((~new).sum()):,} nuclei kept")


if __name__ == "__main__":
    main()
