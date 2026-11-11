#!/usr/bin/env python
"""Re-derive step 1.3's doublet calls with Scrublet RATE-MATCHED instead of auto-thresholded.

Why this exists
---------------
Scrublet's automatic threshold comes from a KDE minimum in the simulated-doublet score
histogram. When that histogram is not cleanly bimodal it fails, and it failed in both
directions on this data:

  Olaniru pools   threshold landed AT the score maximum -> 5-30 calls, every one of them
                  already a scDblFinder call, so the "union of two tools" added 2 cells
                  in 21,387 and the stated safety margin did not exist.
  Migliorini S6   threshold landed far too low -> 3,582 calls against 6-26 in the five
                  sibling samples, removing 1,983 nuclei on Scrublet's word alone (29.6%
                  of the sample, vs 12.3-16.2% everywhere else).

Scrublet's SCORE is not the problem: its AUROC against scDblFinder's calls is 0.856-0.871
on the Olaniru pools, so it ranks doublets well. Only the cut point is unreliable. Rate
matching replaces that cut point with scDblFinder's own observed rate on the same sample --
the rate 1.3 already computes and feeds to Scrublet as expected_doublet_rate -- so Scrublet
contributes a genuine second ranking instead of either nothing or a catastrophe.

This does NOT rerun scDblFinder or Scrublet. Step 1.3 saves both continuous scores in
<sample>.doublets_annotated.h5ad precisely so the union/consensus decision can be revisited,
and that is all this script does. The shared pipeline under single-cell-pipeline/ is untouched.

Ties: scores are discretised (a fraction of simulated-doublet neighbours), so many cells can
share the Nth score. Taking a strict top-N by argsort would split a tie arbitrarily by index;
this keeps every cell at or above the Nth score, which can exceed N slightly and is reported.
"""
import argparse, shutil
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

BASE = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas")


def rate_matched(score, n_target):
    """Cells at or above the n_target-th highest score (ties included, so possibly > n_target)."""
    if n_target <= 0:
        return np.zeros(len(score), bool)
    if n_target >= len(score):
        return np.ones(len(score), bool)
    thr = np.partition(score, -n_target)[-n_target]
    return score >= thr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", required=True, help="dataset dirs under fetal_pancreas")
    ap.add_argument("--apply", nargs="*", default=[],
                    help="of those, the ones to actually rewrite; others are reported only")
    a = ap.parse_args()

    rows = []
    for ds in a.datasets:
        d = BASE / ds / "qc_filtered" / "detect_doublets"
        h5d = d / "h5ad"
        do_apply = ds in a.apply
        for p in sorted(h5d.glob("*.doublets_annotated.h5ad")):
            sid = p.name.split(".")[0]
            adata = ad.read_h5ad(p)
            o = adata.obs
            sd = o["scDblFinder_predicted_doublet"].to_numpy().astype(bool)
            sc_auto = o["scrublet_predicted_doublet"].to_numpy().astype(bool)
            sc_score = o["scrublet_doublet_score"].to_numpy()
            n = adata.n_obs
            N = int(sd.sum())

            sc_rm = rate_matched(sc_score, N)
            union_auto = sd | sc_auto
            union_rm = sd | sc_rm
            cons_rm = sd & sc_rm

            rows.append({
                "dataset": ds, "sample": sid, "cells": n,
                "scDbl": N, "scDbl_%": round(100*N/n, 1),
                "scrub_auto": int(sc_auto.sum()), "scrub_rm": int(sc_rm.sum()),
                "overlap_rm": int(cons_rm.sum()),
                "overlap_%_of_scDbl": round(100*cons_rm.sum()/N, 1) if N else None,
                "union_auto": int(union_auto.sum()), "union_auto_%": round(100*union_auto.mean(), 1),
                "union_rm": int(union_rm.sum()), "union_rm_%": round(100*union_rm.mean(), 1),
                "consensus_rm": int(cons_rm.sum()), "consensus_rm_%": round(100*cons_rm.mean(), 1),
                "kept_union_auto": int((~union_auto).sum()),
                "kept_union_rm": int((~union_rm).sum()),
                "kept_consensus_rm": int((~cons_rm).sum()),
                "applied": do_apply,
            })

            if not do_apply:
                continue

            # keep the auto-threshold result rather than destroying it
            removed = h5d / f"{sid}.doublets_removed.h5ad"
            backup = h5d / f"{sid}.doublets_removed_autothresh.h5ad"
            if removed.exists() and not backup.exists():
                shutil.move(str(removed), str(backup))
            ann_backup = h5d / f"{sid}.doublets_annotated_autothresh.h5ad"
            if not ann_backup.exists():
                shutil.copy2(str(p), str(ann_backup))

            adata.obs["scrublet_predicted_doublet_ratematched"] = sc_rm
            adata.obs["doublet_filter_union_ratematched"] = union_rm
            adata.obs["doublet_filter_consensus_ratematched"] = cons_rm
            adata.write_h5ad(p)
            adata[~union_rm].copy().write_h5ad(removed)
            print(f"  rewrote {sid}: {n:,} -> {int((~union_rm).sum()):,} "
                  f"({int(union_rm.sum()):,} doublets, union of scDblFinder + rate-matched Scrublet)",
                  flush=True)

    df = pd.DataFrame(rows)
    pd.set_option("display.width", 260)
    cols = ["dataset", "sample", "cells", "scDbl", "scDbl_%", "scrub_auto", "scrub_rm",
            "overlap_rm", "overlap_%_of_scDbl", "union_auto", "union_auto_%",
            "union_rm", "union_rm_%", "consensus_rm_%"]
    print("\n=== Scrublet auto-threshold vs rate-matched ===")
    print(df[cols].to_string(index=False))
    print("\n=== cells kept under each rule ===")
    for ds, g in df.groupby("dataset", sort=False):
        print(f"  {ds}: union_auto {g.kept_union_auto.sum():,} | "
              f"union_ratematched {g.kept_union_rm.sum():,} | "
              f"consensus_ratematched {g.kept_consensus_rm.sum():,} "
              f"(of {g.cells.sum():,} post-QC)")
    out = BASE / "doublet_rule_comparison.csv"
    df.to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
