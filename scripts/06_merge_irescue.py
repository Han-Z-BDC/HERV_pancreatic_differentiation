#!/usr/bin/env python
"""
Merge the 15 IRescue subfamily matrices into one AnnData and take a first look.

Inputs : $TE_OUT/irescue_subfamily/<sample>/counts/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
         shared annotation table (ATAC_WNN_annotation, stage, include_downstream)
         per-cell gene UMIs (total_counts) from the shared merged RNA object
         RepeatMasker subfamily -> family / class map
Outputs: $TE_OUT/irescue_subfamily/irescue_subfamily_merged.h5ad   (cells x TE subfamilies, raw counts)
         $TE_OUT/qc/irescue_*.csv / .png
   - TE UMIs per cell and per 10k gene UMIs, by stage
   - pseudobulk TE expression (UMIs per million gene UMIs) per stage and per annotation group
   - HERV-H / HERV-K / L1HS / SVA check across stages
   - heatmap of the most variable ERV (LTR-class) subfamilies across annotation groups

Note: TE reads can also be gene reads (intronic TEs in single-nucleus data), so "per gene UMIs" is a
normalisation, not a TE fraction of the transcriptome.
"""
import os
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad

REPO = Path(__file__).resolve().parents[1]
TE_OUT = Path(os.environ["TE_OUT"]); ANNOT = Path(os.environ["ANNOTATION_CSV"])
SCP_OUT = Path(os.environ["SCP_OUT"]); TE_ANNOT = Path(os.environ["TE_ANNOT_DIR"])
qc = TE_OUT / "qc"; qc.mkdir(parents=True, exist_ok=True)
STAGES = ["D0", "D4", "D7", "D9", "D12", "D15", "D22", "D45"]
GROUP_ORDER = ["hPSC", "DE", "GT", "OTX2_endo", "FG", "liver", "PP1", "PP2", "ductal", "exocrine", "ENP1", "ENP2",
               "ENP3", "early_SC_beta", "early_EC_beta", "early_SC_alpha", "early_SC_delta", "early_SC_EC",
               "late_SC_beta", "late_EC_beta", "late_SC_alpha", "late_SC_delta", "late_SC_EC", "stroma"]
FOCUS = ["HERVH-int", "LTR7", "LTR7Y", "LTR7B", "HERVK-int", "LTR5_Hs", "HERV9-int", "LTR12C", "MER41B",
         "HERVE-int", "HERVW-int", "LTR17", "L1HS", "L1PA2", "SVA_F", "SVA_D"]

samples = pd.read_csv(REPO / "configs/samples.tsv", sep="\t")
ann = pd.read_csv(ANNOT, index_col=0, dtype=str)

# ---- merge
mats, obs_names, feats = [], [], None
for s in samples["sample"]:
    d = TE_OUT / "irescue_subfamily" / s / "counts"
    m = scipy.io.mmread(d / "matrix.mtx.gz").tocsr().T          # features x cells -> cells x features
    f = pd.read_csv(d / "features.tsv.gz", sep="\t", header=None)[0].astype(str).values
    b = pd.read_csv(d / "barcodes.tsv.gz", header=None)[0].astype(str).values
    df = pd.DataFrame.sparse.from_spmatrix(m, columns=f)  # align feature order across samples
    if feats is None:
        feats = list(f)
    else:
        feats = list(dict.fromkeys(feats + list(f)))
    mats.append((m, f)); obs_names += [f"{s}#{x}" for x in b]
    print(f"{s}: {m.shape[0]} cells, {m.shape[1]} subfamilies, {int(m.sum()):,} TE UMIs", flush=True)
fidx = {f: i for i, f in enumerate(feats)}
blocks = []
for m, f in mats:
    cols = np.array([fidx[x] for x in f])
    m = m.tocoo()
    blocks.append(sp.csr_matrix((m.data, (m.row, cols[m.col])), shape=(m.shape[0], len(feats))))
X = sp.vstack(blocks).tocsr().astype(np.float32)

rm = pd.read_csv(TE_ANNOT / "rmsk_hg38.bed.gz", sep="\t", header=None, usecols=[3, 6, 7],
                 names=["subfamily", "family", "class"]).drop_duplicates("subfamily").set_index("subfamily")
var = pd.DataFrame(index=feats)
var["family"] = rm["family"].reindex(feats).fillna("unknown").values
var["class"] = rm["class"].reindex(feats).fillna("unknown").str.rstrip("?").values

obs = ann.reindex(obs_names)[["stage", "cluster", "single_stage_label", "ATAC_WNN_annotation", "include_downstream"]].copy()
obs["sample"] = [n.split("#")[0] for n in obs_names]
with h5py.File(SCP_OUT / "integration/h5ad/merged_rna_scvi_integrated.h5ad", "r") as h:
    idx = h["obs"]["_index"][:].astype(str)
    gene_umis = pd.Series(h["obs"]["total_counts"][:], index=idx)
obs["gene_umis"] = gene_umis.reindex(obs_names).values
te = ad.AnnData(X=X, obs=obs, var=var)
te.obs["te_umis"] = np.asarray(X.sum(1)).ravel()
te.obs["te_per_10k_gene"] = 1e4 * te.obs.te_umis / te.obs.gene_umis
print(f"merged: {te.n_obs} cells x {te.n_vars} subfamilies; missing annotation: {te.obs.stage.isna().sum()}", flush=True)
te.write_h5ad(TE_OUT / "irescue_subfamily" / "irescue_subfamily_merged.h5ad", compression="gzip")
te = te[te.obs.include_downstream.str.lower() == "true"].copy()

# ---- TE load by stage and class
load = te.obs.groupby("stage")[["te_umis", "gene_umis", "te_per_10k_gene"]].median().reindex(STAGES).round(1)
cls = pd.DataFrame(te.X.toarray(), columns=te.var_names).T.groupby(te.var["class"].values).sum().T
cls_share = cls.groupby(te.obs.stage.values).sum().div(cls.groupby(te.obs.stage.values).sum().sum(1), axis=0).mul(100).round(1).reindex(STAGES)
load.to_csv(qc / "irescue_te_load_by_stage.csv"); cls_share.to_csv(qc / "irescue_class_share_by_stage.csv")
pd.set_option("display.width", 250)
print("TE load by stage (medians per cell):\n", load.to_string())
print("share of TE UMIs by class (%):\n", cls_share[[c for c in ["SINE", "LINE", "LTR", "DNA", "Retroposon", "RC"] if c in cls_share]].to_string(), flush=True)


def pseudobulk_cpm(by):
    """TE UMIs per million gene UMIs, summed per group."""
    g = te.obs[by].values
    cnt = pd.DataFrame(te.X.toarray(), columns=te.var_names).groupby(g).sum()
    lib = te.obs.groupby(by)["gene_umis"].sum()
    return cnt.div(lib.reindex(cnt.index), axis=0) * 1e6


cpm_stage = pseudobulk_cpm("stage").reindex(STAGES)
cpm_group = pseudobulk_cpm("ATAC_WNN_annotation")
cpm_group = cpm_group.reindex([g for g in GROUP_ORDER if g in cpm_group.index])
cpm_stage.T.to_csv(qc / "irescue_cpm_by_stage.csv"); cpm_group.T.to_csv(qc / "irescue_cpm_by_group.csv")
focus = [f for f in FOCUS if f in cpm_stage.columns]
print("selected young / pluripotency-linked elements, UMIs per million gene UMIs, by stage:\n",
      cpm_stage[focus].T.round(1).to_string(), flush=True)

# ---- plots
fig, axes = plt.subplots(1, 2, figsize=(15, 4.5))
ax = axes[0]
data = [te.obs.loc[te.obs.stage == s, "te_per_10k_gene"].dropna().values for s in STAGES]
ax.boxplot(data, showfliers=False); ax.set_xticks(range(1, 9), STAGES); ax.set_ylabel("TE UMIs per 10k gene UMIs")
ax.set_title("TE load per cell by stage (IRescue, all TE classes)")
ax = axes[1]
for f in ["HERVH-int", "LTR7", "HERVK-int", "LTR5_Hs", "L1HS", "SVA_F"]:
    if f in cpm_stage:
        ax.plot(STAGES, cpm_stage[f].values, marker="o", label=f)
ax.set_yscale("log"); ax.set_ylabel("UMIs per million gene UMIs (pseudobulk)"); ax.legend(fontsize=9, frameon=False)
ax.set_title("Young elements across stages")
fig.tight_layout(); fig.savefig(qc / "irescue_stage_overview.png", dpi=100, bbox_inches="tight"); plt.close(fig)

# heatmap: most variable ERV (LTR-class) subfamilies across groups, among expressed ones
ltr = te.var_names[(te.var["class"] == "LTR").values]
sub = cpm_group[ltr]
sub = sub.loc[:, sub.max() >= 5]                       # at least 5 UMIs per million gene UMIs in some group
lg = np.log2(sub + 1)
top = lg.var().sort_values(ascending=False).index[:60]
z = (lg[top] - lg[top].mean()) / lg[top].std()
order = z.T.idxmax(axis=1).map({g: i for i, g in enumerate(z.index)}).sort_values().index  # sort by peak group
z = z[order]
fig, ax = plt.subplots(figsize=(14, 15))
im = ax.imshow(z.T.values, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
ax.set_xticks(range(len(z.index)), z.index, rotation=60, ha="right", fontsize=9)
ax.set_yticks(range(len(order)), [f"{o} ({te.var.loc[o, 'family']})" for o in order], fontsize=8)
fig.colorbar(im, ax=ax, fraction=0.025, label="z-score of log2(UMIs per million gene UMIs + 1)")
ax.set_title("60 most variable expressed ERV / LTR subfamilies across annotation groups (IRescue pseudobulk)")
fig.savefig(qc / "irescue_ERV_heatmap_by_group.png", dpi=80, bbox_inches="tight"); plt.close(fig)
z.T.round(3).to_csv(qc / "irescue_ERV_heatmap_values.csv")
peak = z.idxmax().rename("peak_group")
print("top variable ERV subfamilies and the group where each peaks:\n",
      pd.concat([peak, sub[order].max().round(1).rename("max_cpm")], axis=1).to_string(), flush=True)
print("done", qc, flush=True)
