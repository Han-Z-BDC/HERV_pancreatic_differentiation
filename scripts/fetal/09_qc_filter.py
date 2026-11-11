#!/usr/bin/env python
"""Apply the agreed QC to the Olaniru 2023 pools and write the cell whitelists the TE track uses.

Settings were chosen by looking at the distributions first (07_qc_metrics_only.py ->
08_qc_report.py -> docs/fetal_qc_review.html), not by taking a default:

  1. Erythrocyte removal FIRST.  Cells with pct_counts_hb > 50 are fetal erythrocytes -- in pool4
     they are 25.4% of all cells, with a median of 94.8% haemoglobin, 171 genes, 0.03% MT and 0.5%
     ribo.  They are real cells, not debris, but they wreck the MAD statistics: pool4's gene MAD is
     0.609 with them and 0.380 without, which drags its 3-MAD gene cutoff from 610 down to 230 and
     lets actual debris through.  Removing them first makes one global n_mads work for every pool.
  2. n_mads = 3 on log1p(total_counts) and log1p(n_genes_by_counts), one-sided (low only), exactly
     as single-cell-pipeline/scripts/1.1_rna_qc.py:is_mad_outlier().  The repo default of 5 does not
     transfer here: at 4-5 MAD the counts criterion removes ~0 cells in every pool, i.e. it is not
     filtering at all, while the plots show a clear low-quality left shoulder.
  3. %MT: median + 5*MAD clamped to [3, 15], as compute_mt_threshold().  On these pools it resolves
     to 6.7-14.1% and neither clamp engages.
  4. %ribo: computed and reported, never filtered (median 17-19% here; high ribosomal content is
     expected in translationally active progenitors).
  5. A GEO donor label is required for the whitelist -- a cell the authors could not assign to a
     gestational age cannot be placed in the atlas.  Reported separately from the QC loss.

Writes, under $FETAL_ROOT/qc_filtered/:
  whitelists/<pool>.txt      barcodes for IRescue (05_irescue_fetal.sb reads these)
  qc_filter_summary.csv      per-pool counts at every stage, and the cutoffs actually used
  per_cell_flags.csv.gz      one row per cell, every metric and flag, for traceability
  figures/filter_<pool>.png  what was removed and why
  rna_qc/h5ad/<pool>.qc_filtered.h5ad    kept cells, in the exact layout
  rna_qc/h5ad/<pool>.qc_annotated.h5ad   all cells + flags
     These two are the artifacts single-cell-pipeline/scripts/1.1_rna_qc.py would have written.
     This script stands in for that step (it adds the erythrocyte removal 1.1 has no notion of),
     so emitting the same files lets 1.2_remove_background.sb and 1.3_detect_doublets.sb run
     against the fetal data unmodified -- see scp_fetal/config.sh.

Usage: python 09_qc_filter.py --samples configs/fetal_olaniru_samples.tsv --outdir <dir>
"""
import argparse, os, re, sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N_MADS = 3.0
MT_N_MADS, MT_FLOOR, MT_CAP = 5.0, 3.0, 15.0
RIBO_MAX = None
HB_PCT_MAX = 50.0
HB_GENES = ["HBB", "HBA1", "HBA2", "HBG1", "HBG2", "HBD", "HBM", "HBZ", "HBE1",
            "ALAS2", "SLC4A1", "AHSP", "CA1"]


def mad(v):
    return float(np.median(np.abs(v - np.median(v))))


def low_cut(v, n_mads):
    d = mad(v)
    return None if d == 0 else float(np.median(v) - n_mads * d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cellranger", required=True, help="dir holding <pool>/outs/")
    ap.add_argument("--geo-suppl", default=None,
                    help="dir with the GEO *_demultiplexing_*.csv.gz; omit for a dataset that has "
                         "one donor per library and therefore no demultiplexing table")
    ap.add_argument("--mt-fixed", type=float, default=None,
                    help="use a FIXED %MT cutoff instead of the median+n*MAD rule. Required for "
                         "nuclei: their median %MT is 0.02-0.10%%, so MAD collapses to the "
                         "quantization of a single mitochondrial read and the MAD rule would "
                         "discard 8-20%% of nuclei, enriched 1.5-2.3x for immune/endothelial. "
                         "Pass 0 to disable the %MT criterion entirely.")
    a = ap.parse_args()
    out = Path(a.outdir)
    (out / "whitelists").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    samples = pd.read_csv(a.samples, sep="\t")
    rows, flags = [], []
    for _, srow in samples.iterrows():
        pool, demux = srow["pool"], srow["geo_demux"]
        h5 = Path(a.cellranger) / pool / "outs" / "filtered_feature_bc_matrix.h5"
        if not h5.exists():
            print(f"  skip {pool}: {h5} not found", file=sys.stderr); continue
        print(f"  {pool} ...", flush=True)
        ad = sc.read_10x_h5(h5); ad.var_names_make_unique()
        if "feature_types" in ad.var:
            ad = ad[:, ad.var["feature_types"] == "Gene Expression"].copy()
        ad.var["mt"] = ad.var_names.str.startswith("MT-")
        ad.var["ribo"] = ad.var_names.str.startswith(("RPS", "RPL"))
        ad.var["hb"] = ad.var_names.isin(HB_GENES)
        sc.pp.calculate_qc_metrics(ad, qc_vars=["mt", "ribo", "hb"], percent_top=None,
                                   log1p=True, inplace=True)
        o = ad.obs
        n0 = ad.n_obs

        # 1. erythrocytes
        ery = (o["pct_counts_hb"] > HB_PCT_MAX).to_numpy()

        # 2-3. MAD thresholds, computed on NON-erythrocytes only
        keep_for_stats = ~ery
        lc = o["log1p_total_counts"].to_numpy()
        lg = o["log1p_n_genes_by_counts"].to_numpy()
        mt = o["pct_counts_mt"].to_numpy()
        cc = low_cut(lc[keep_for_stats], N_MADS)
        cg = low_cut(lg[keep_for_stats], N_MADS)
        sub_mt = mt[keep_for_stats]
        cmt_raw = float(np.median(sub_mt) + MT_N_MADS * mad(sub_mt))
        if a.mt_fixed is None:
            cmt = min(max(cmt_raw, MT_FLOOR), MT_CAP)
        else:
            cmt = float(a.mt_fixed)

        fail_c = lc < cc if cc is not None else np.zeros(n0, bool)
        fail_g = lg < cg if cg is not None else np.zeros(n0, bool)
        fail_mt = np.zeros(n0, bool) if cmt == 0 else mt > cmt
        fail_rb = (o["pct_counts_ribo"].to_numpy() > RIBO_MAX) if RIBO_MAX else np.zeros(n0, bool)
        fail_qc = fail_c | fail_g | fail_mt | fail_rb

        # 5. GEO donor label (only for a multiplexed dataset that ships one)
        if a.geo_suppl is None or not isinstance(demux, str) or not demux or demux == "nan":
            # One donor per library, so the library IS the donor -- and it has to be labelled as
            # such, not left as None: an all-None object column cannot be written to h5ad, and
            # two libraries can share a gestational age while being different embryos
            # (S1_14PCW and S2_14PCW), so the age alone would merge two donors into one.
            lab = np.array([pool] * n0, dtype=object)
            has_donor = np.ones(n0, bool)
        else:
            geo = pd.read_csv(Path(a.geo_suppl) / demux)
            geo.columns = ["barcode", "geo_label"]
            geo["geo_label"] = geo["geo_label"].astype(str).str.strip()
            lab = pd.Series(o.index, index=o.index).map(
                dict(zip(geo["barcode"], geo["geo_label"]))).to_numpy()
        # A donor label is "<weeks>PCW", optionally with a suffix when one pool holds two embryos of
        # the SAME gestational age -- pool1 has both "13PCW" and "13PCW_1" (2,712 and 2,722 cells),
        # so testing endswith("PCW") would silently discard half that pool. The non-donor values are
        # "Doublet"/"Unassigned", whose capitalisation differs between pools ("doublet" in pool3).
            donor_re = re.compile(r"^\d+PCW(_\d+)?$")
            has_donor = np.array([isinstance(x, str) and bool(donor_re.match(x)) for x in lab])

        keep = (~ery) & (~fail_qc) & has_donor
        np.savetxt(out / "whitelists" / f"{pool}.txt", o.index[keep].to_numpy(), fmt="%s")

        ad.obs["donor"] = lab
        ad.obs["erythrocyte"] = ery
        ad.obs["fail_counts"] = fail_c
        ad.obs["fail_genes"] = fail_g
        ad.obs["fail_mt"] = fail_mt
        ad.obs["passed_qc_step1_2"] = keep   # name matches 1.1_rna_qc.py's own column
        h5d = out / "rna_qc" / "h5ad"; h5d.mkdir(parents=True, exist_ok=True)
        ad.write_h5ad(h5d / f"{pool}.qc_annotated.h5ad")
        ad[keep].copy().write_h5ad(h5d / f"{pool}.qc_filtered.h5ad")

        flags.append(pd.DataFrame({
            "pool": pool, "barcode": o.index, "geo_label": lab,
            "total_counts": o["total_counts"].to_numpy(),
            "n_genes": o["n_genes_by_counts"].to_numpy(),
            "pct_mt": mt, "pct_ribo": o["pct_counts_ribo"].to_numpy(),
            "pct_hb": o["pct_counts_hb"].to_numpy(),
            "erythrocyte": ery, "fail_counts": fail_c, "fail_genes": fail_g,
            "fail_mt": fail_mt, "no_donor": ~has_donor, "keep": keep}))

        rows.append({
            "pool": pool, "cells_called": n0,
            "erythrocytes": int(ery.sum()),
            "after_ery": int((~ery).sum()),
            "umi_cut": round(float(np.expm1(cc)), 1) if cc is not None else None,
            "gene_cut": round(float(np.expm1(cg)), 1) if cg is not None else None,
            "mt_cut_raw": round(cmt_raw, 2), "mt_cut_used": round(cmt, 2),
            "mt_rule": ("fixed" if a.mt_fixed is not None else
                        ("floor" if cmt > cmt_raw + 1e-9 else
                         ("cap" if cmt < cmt_raw - 1e-9 else "mad"))),
            "fail_counts": int((fail_c & ~ery).sum()), "fail_genes": int((fail_g & ~ery).sum()),
            "fail_mt": int((fail_mt & ~ery).sum()),
            "fail_qc_any": int((fail_qc & ~ery).sum()),
            "no_donor_label": int((~has_donor & ~ery & ~fail_qc).sum()),
            "kept": int(keep.sum()),
            "kept_pct_of_called": round(100 * keep.sum() / n0, 1),
        })

        # figure
        fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
        cat = np.where(ery, "erythrocyte", np.where(fail_qc, "fail QC",
                       np.where(~has_donor, "no donor label", "kept")))
        col = {"kept": "#2f6fb5", "fail QC": "#d1495b", "erythrocyte": "#e0a44a",
               "no donor label": "#9aa5b1"}
        for ax, (xv, yv, xl, yl) in zip(axes, [
                (lc, lg, "log1p(total counts)", "log1p(n genes)"),
                (lg, mt, "log1p(n genes)", "% MT"),
                (lg, o["pct_counts_hb"].to_numpy(), "log1p(n genes)", "% haemoglobin")]):
            for k, c in col.items():
                m = cat == k
                if m.sum(): ax.scatter(xv[m], yv[m], s=2, c=c, alpha=.5, label=f"{k} ({m.sum():,})",
                                       linewidths=0, rasterized=True)
            ax.set_xlabel(xl, fontsize=9); ax.set_ylabel(yl, fontsize=9); ax.tick_params(labelsize=8)
            for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        if cc is not None: axes[0].axvline(cc, ls="--", c="#d1495b", lw=1.2)
        if cg is not None:
            axes[0].axhline(cg, ls="--", c="#d1495b", lw=1.2)
            axes[1].axvline(cg, ls="--", c="#d1495b", lw=1.2); axes[2].axvline(cg, ls="--", c="#d1495b", lw=1.2)
        axes[1].axhline(cmt, ls="--", c="#d1495b", lw=1.2)
        axes[2].axhline(HB_PCT_MAX, ls="--", c="#e0a44a", lw=1.2)
        axes[0].legend(fontsize=7, markerscale=3, frameon=False, loc="lower right")
        fig.suptitle(f"{pool}: {int(keep.sum()):,} of {n0:,} kept "
                     f"(n_mads={N_MADS}, %MT<{cmt:.1f}, %HB<{HB_PCT_MAX:.0f})", fontsize=10, y=1.03)
        fig.savefig(out / "figures" / f"filter_{pool}.png", dpi=130, bbox_inches="tight")
        plt.close(fig)

    if not rows: sys.exit("no pools processed")
    df = pd.DataFrame(rows)
    df.to_csv(out / "qc_filter_summary.csv", index=False)
    pd.concat(flags).to_csv(out / "per_cell_flags.csv.gz", index=False)
    print("\n" + df.to_string(index=False))
    print(f"\nTOTAL called {df.cells_called.sum():,} -> kept {df.kept.sum():,} "
          f"({100*df.kept.sum()/df.cells_called.sum():.1f}%)")
    print(f"wrote whitelists to {out}/whitelists/")


if __name__ == "__main__":
    main()
