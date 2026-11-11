#!/usr/bin/env python
"""Write the FINAL cell whitelists -- after QC (09) + SoupX (1.2) + doublet removal (1.3).

09_qc_filter.py writes qc_filtered/whitelists/ (post-QC, 21,387 cells for Olaniru). Those are what
IRescue was first run against. Doublet detection then removed a further 1,557 cells, so the TE
matrices built from the earlier whitelists contain doublets.

This reads the doublet-removed h5ads from the lab pipeline's step 1.3 and writes
qc_filtered/whitelists_final/<sample>.txt, which 05_irescue_fetal.sb consumes.

(IRescue deduplicates UMIs and runs its EM per cell barcode, so a cell's TE counts do not depend on
which other barcodes are in the whitelist -- subsetting the existing matrices to these barcodes
would give the same numbers. Re-running is just the unambiguous version.)

Usage: python 10_final_whitelists.py --doublet-dir <...>/detect_doublets/h5ad --outdir <...>
"""
import argparse
from pathlib import Path

import anndata
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--doublet-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    total = 0
    for s in a.samples:
        h5 = Path(a.doublet_dir) / f"{s}.doublets_removed.h5ad"
        if not h5.exists():
            raise SystemExit(f"missing {h5} -- run the pipeline's 1.3_detect_doublets.sb first")
        ad = anndata.read_h5ad(h5, backed="r")
        bc = np.asarray(ad.obs_names)
        np.savetxt(out / f"{s}.txt", bc, fmt="%s")
        total += len(bc)
        print(f"  {s:10} {len(bc):>7,} cells -> {out / (s + '.txt')}")
    print(f"  {'TOTAL':10} {total:>7,} cells")


if __name__ == "__main__":
    main()
