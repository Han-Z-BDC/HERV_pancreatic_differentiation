"""Is Olaniru's rho = 0.01 a real estimate, or SoupX's contaminationRange floor?

All four pools reported exactly 0.01, which is the lower bound of autoEstCont's default
contaminationRange = c(0.01, 0.8). Rerunning the SAME estimation with the bound dropped to
0.0001 distinguishes the two: if rho stays ~0.01 it was a genuine estimate; if it falls to
e.g. 0.002 it was clamped, and the pipeline removed several times more counts than the data
supported. Nothing is written -- this only re-estimates and prints.
"""
import sys, numpy as np, scanpy as sc, anndata2ri, rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter

anndata2ri.activate()
r("suppressMessages(library(SoupX))")

BASE = "/scratch/alpine/hzhu10@xsede.org/fetal_pancreas"
CASES = [("olaniru2023", p) for p in ("pool1", "pool2", "pool3", "pool4")]

print(f"{'sample':10} {'rho default [0.01,0.8]':>24} {'rho with floor 1e-4':>22}  verdict")
for ds, sid in CASES:
    raw = sc.read_10x_h5(f"{BASE}/{ds}/cellranger/{sid}/outs/raw_feature_bc_matrix.h5")
    raw.var_names_make_unique()
    if "feature_types" in raw.var:
        raw = raw[:, raw.var["feature_types"] == "Gene Expression"].copy()
    flt = sc.read_h5ad(f"{BASE}/{ds}/qc_filtered/rna_qc/h5ad/{sid}.qc_filtered.h5ad")
    raw = raw[:, flt.var_names].copy()

    # the throwaway clustering scaffold SoupX needs to group cells, same recipe as 1.2
    a = flt.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=3000); a = a[:, a.var.highly_variable].copy()
    sc.pp.scale(a, max_value=10); sc.tl.pca(a, n_comps=50)
    sc.pp.neighbors(a, n_neighbors=30); sc.tl.leiden(a, resolution=1.0, key_added="cl")

    ro.globalenv["tod"] = raw.X.T.tocsc()
    ro.globalenv["toc"] = flt.X.T.tocsc()
    ro.globalenv["cl"] = ro.StrVector(a.obs["cl"].astype(str).tolist())
    # the matrices cross into R without dimnames, so colnames(toc) is NULL and setClusters
    # sees NAs -- the barcodes have to be attached explicitly on the R side
    ro.globalenv["bc"] = ro.StrVector(flt.obs_names.tolist())
    ro.globalenv["gn"] = ro.StrVector(flt.var_names.tolist())
    out = r("""
      colnames(toc) <- bc; rownames(toc) <- gn
      colnames(tod) <- paste0("d", seq_len(ncol(tod))); rownames(tod) <- gn
      sc <- SoupChannel(tod, toc)
      sc <- setClusters(sc, setNames(cl, bc))
      a <- tryCatch({ s<-autoEstCont(sc, doPlot=FALSE); s$metaData$rho[1] }, error=function(e) NA)
      b <- tryCatch({ s<-autoEstCont(sc, doPlot=FALSE, contaminationRange=c(1e-4,0.8));
                      s$metaData$rho[1] }, error=function(e) NA)
      c(a, b)
    """)
    d, f = float(out[0]), float(out[1])
    verdict = "CLAMPED -- over-corrected" if f < d * 0.8 else "genuine estimate"
    print(f"{sid:10} {d:24.5f} {f:22.5f}  {verdict}")
