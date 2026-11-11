"""Is the %MT-failing nuclei population a cell type or read noise?

If the nuclei a 5-MAD (unclamped) %MT cut would remove are enriched for a lineage marker,
the cut is deleting biology and the cost is not 12% of random cells. Also reports %Hb,
which 07 never computed, to settle whether nucleated erythroid precursors are present.
"""
import numpy as np, pandas as pd, scanpy as sc
from pathlib import Path

R = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas/migliorini2024/cellranger")
HB = ["HBB","HBA1","HBA2","HBG1","HBG2","HBD","HBM","HBZ","HBE1","ALAS2","SLC4A1","AHSP","CA1"]
MARK = {"beta":"INS","alpha":"GCG","delta":"SST","acinar":"PRSS1","acinar2":"CPA1",
        "ductal":"KRT19","mesen":"COL1A1","immune":"PTPRC","endo":"PECAM1","prolif":"MKI67"}

rows = []
for s in sorted(p.name for p in R.iterdir() if p.is_dir()):
    h5 = R / s / "outs" / "filtered_feature_bc_matrix.h5"
    if not h5.exists():
        continue
    a = sc.read_10x_h5(h5); a.var_names_make_unique()
    if "feature_types" in a.var:
        a = a[:, a.var["feature_types"] == "Gene Expression"].copy()
    a.var["mt"] = a.var_names.str.startswith("MT-")
    a.var["hb"] = a.var_names.isin(HB)
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt","hb"], percent_top=None, log1p=False, inplace=True)
    v = a.obs["pct_counts_mt"].to_numpy()
    med = np.median(v); m = np.median(np.abs(v - med))
    fail = v > med + 5*m                      # the cut under discussion: 5 MAD, no floor
    tot = a.obs["total_counts"].to_numpy()
    r = {"sample": s, "n": a.n_obs, "fail_mt5": int(fail.sum()),
         "fail_pct": round(100*fail.mean(), 1),
         "umi_pass": int(np.median(tot[~fail])), "umi_fail": int(np.median(tot[fail])),
         "gene_pass": int(np.median(a.obs["n_genes_by_counts"].to_numpy()[~fail])),
         "gene_fail": int(np.median(a.obs["n_genes_by_counts"].to_numpy()[fail])),
         "hb_gt50": int((a.obs["pct_counts_hb"].to_numpy() > 50).sum()),
         "hb_gt10": int((a.obs["pct_counts_hb"].to_numpy() > 10).sum()),
         "hb_p999": round(float(np.percentile(a.obs["pct_counts_hb"].to_numpy(), 99.9)), 2)}
    # marker detection rate in failing vs passing nuclei -- enrichment means a lineage, not noise
    sc.pp.normalize_total(a, target_sum=1e4)
    for lab, g in MARK.items():
        if g not in a.var_names:
            continue
        x = np.asarray(a[:, g].X.todense()).ravel()
        pf = 100*(x[fail] > 0).mean(); pp = 100*(x[~fail] > 0).mean()
        r[f"{lab}_fail%"] = round(pf, 1); r[f"{lab}_pass%"] = round(pp, 1)
        r[f"{lab}_OR"] = round(pf/pp, 2) if pp > 0 else np.nan
    rows.append(r)
    print(f"  {s} done", flush=True)

df = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print("\n=== depth of the %MT-failing nuclei (5 MAD, no floor) ===")
print(df[["sample","n","fail_mt5","fail_pct","umi_pass","umi_fail","gene_pass","gene_fail"]].to_string(index=False))
print("\n=== hemoglobin (nucleated erythroid precursors?) ===")
print(df[["sample","hb_gt50","hb_gt10","hb_p999"]].to_string(index=False))
print("\n=== marker detection rate, failing vs passing (OR>1.5 = the cut targets that lineage) ===")
cols = ["sample"] + [c for c in df.columns if c.endswith("_OR")]
print(df[cols].to_string(index=False))
