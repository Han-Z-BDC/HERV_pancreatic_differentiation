"""Is Scrublet's SCORE uninformative here, or just mis-thresholded?

The union rule in 1.3 adds ~0 cells because Scrublet calls 5-30 doublets where scDblFinder
calls 317-474. Two very different causes:
  (a) the score ranks doublets fine but the automatic threshold is set absurdly high  -> fixable
  (b) the score carries no signal on this data                                        -> drop it
Distinguished by AUROC of the scrublet score against scDblFinder's calls, and by how many
scDblFinder doublets fall in the top-N scrublet scores (N = scDblFinder's own count).
"""
import numpy as np, pandas as pd, anndata as ad
from pathlib import Path

D = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas/olaniru2023/qc_filtered/detect_doublets/h5ad")

def auroc(score, label):
    """Mann-Whitney U / (n_pos*n_neg) -- rank-based, ties handled."""
    r = pd.Series(score).rank().to_numpy()
    npos, nneg = label.sum(), (~label).sum()
    if npos == 0 or nneg == 0:
        return np.nan
    return (r[label].sum() - npos*(npos+1)/2) / (npos*nneg)

rows = []
for p in sorted(D.glob("*.doublets_annotated.h5ad")):
    pool = p.name.split(".")[0]
    a = ad.read_h5ad(p)
    o = a.obs
    sc_score = o["scrublet_doublet_score"].to_numpy()
    sd_score = o["scDblFinder_doublet_score"].to_numpy()
    sd_call = o["scDblFinder_predicted_doublet"].to_numpy().astype(bool)
    sc_call = o["scrublet_predicted_doublet"].to_numpy().astype(bool)
    n = len(o); N = int(sd_call.sum())

    # implied scrublet threshold: lowest score among its own calls
    thr = float(sc_score[sc_call].min()) if sc_call.any() else np.nan
    # top-N by scrublet score -- what scrublet WOULD call if rate-matched to scDblFinder
    topN = np.zeros(n, bool); topN[np.argsort(-sc_score)[:N]] = True

    rows.append({
        "pool": pool, "n": n, "scDbl": N, "scrub": int(sc_call.sum()),
        "scrub_thr": round(thr, 4),
        "scrub_score_med": round(float(np.median(sc_score)), 4),
        "scrub_score_p99": round(float(np.percentile(sc_score, 99)), 4),
        "scrub_score_max": round(float(sc_score.max()), 4),
        "AUROC_scrub_vs_scDbl": round(auroc(sc_score, sd_call), 3),
        "AUROC_scDbl_self": round(auroc(sd_score, sd_call), 3),
        "topN_recall_%": round(100*(topN & sd_call).sum()/N, 1),
        "topN_overlap": int((topN & sd_call).sum()),
    })
    print(f"  {pool} done", flush=True)

pd.set_option("display.width", 250)
df = pd.DataFrame(rows)
print("\n=== Scrublet score behaviour vs scDblFinder calls ===")
print(df.to_string(index=False))
print("""
Reading it:
  AUROC ~0.5           -> scrublet score is noise on this data (cause b, drop the tool)
  AUROC >0.8           -> score ranks doublets well, threshold is the bug (cause a, fixable)
  topN_recall_%        -> agreement if scrublet were rate-matched instead of auto-thresholded
  scrub_thr vs p99/max -> how far into the tail the automatic threshold was pushed""")
