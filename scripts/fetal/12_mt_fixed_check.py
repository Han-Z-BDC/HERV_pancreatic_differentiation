"""Cost of a FIXED %MT cutoff on the Migliorini nuclei, at 1 / 2 / 3 / 5 %.

A fixed cutoff avoids the MAD rule's quantization collapse. What still has to be checked is
whether it lands on a lineage: the 5-MAD cut removed immune/endothelial nuclei 1.5-2.3x over,
so the same enrichment test is run here at the proposed 1% line.
"""
import numpy as np, pandas as pd, scanpy as sc
from pathlib import Path

R = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas/migliorini2024/cellranger")
MARK = {"beta":"INS","alpha":"GCG","delta":"SST","acinar":"PRSS1","ductal":"KRT19",
        "mesen":"COL1A1","immune":"PTPRC","endo":"PECAM1","prolif":"MKI67"}
CUTS = [1.0, 2.0, 3.0, 5.0]

cost, enr = [], []
for s in sorted(p.name for p in R.iterdir() if p.is_dir()):
    h5 = R / s / "outs" / "filtered_feature_bc_matrix.h5"
    if not h5.exists():
        continue
    a = sc.read_10x_h5(h5); a.var_names_make_unique()
    if "feature_types" in a.var:
        a = a[:, a.var["feature_types"] == "Gene Expression"].copy()
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    v = a.obs["pct_counts_mt"].to_numpy()
    tot = a.obs["total_counts"].to_numpy()
    ng = a.obs["n_genes_by_counts"].to_numpy()

    c = {"sample": s, "n": a.n_obs}
    for t in CUTS:
        f = v > t
        c[f"lost_{t:g}%"] = int(f.sum()); c[f"pct_{t:g}%"] = round(100*f.mean(), 2)
    f1 = v > 1.0
    c["umi_pass"] = int(np.median(tot[~f1]))
    c["umi_fail"] = int(np.median(tot[f1])) if f1.any() else 0
    c["gene_fail"] = int(np.median(ng[f1])) if f1.any() else 0
    cost.append(c)

    sc.pp.normalize_total(a, target_sum=1e4)
    e = {"sample": s, "n_fail_1pct": int(f1.sum())}
    for lab, g in MARK.items():
        if g not in a.var_names or not f1.any():
            continue
        x = np.asarray(a[:, g].X.todense()).ravel()
        pf = (x[f1] > 0).mean(); pp = (x[~f1] > 0).mean()
        e[lab] = round(pf/pp, 2) if pp > 0 else np.nan
    enr.append(e)
    print(f"  {s} done", flush=True)

pd.set_option("display.width", 250)
print("\n=== nuclei lost at each fixed %MT cutoff ===")
print(pd.DataFrame(cost).to_string(index=False))
print("\n=== marker detection ratio, >1% MT vs <=1% MT  (>1.5 = the cut targets that lineage) ===")
print(pd.DataFrame(enr).to_string(index=False))
