#!/usr/bin/env python
"""
HERV9 chr5 provirus vs CDH6 along the PP -> ENP -> endocrine trajectory (idea note:
docs/IDEA_HERV9_CDH6_delamination.md).

Question: does the provirus switch on as CDH6 (and other epithelial adhesion genes) switch off at
ENP1 -> ENP2, i.e. during endocrine delamination? Which comes first, and is the 5' LTR (LTR12B)
accessibility linked to CDH6 or to any other gene in cis?

Steps
  1. Cells: include_downstream == True and ATAC_WNN_annotation in LINEAGE (PP1/PP2 -> ENP1-3 -> early
     endocrine; late groups optional with --with-late). At most --max-per-group cells per group.
  2. Pseudotime: diffusion pseudotime (scanpy) on the scVI latent, root = the PP1 cell at the extreme of
     the first diffusion component on the PP side. Also stores DC1/DC2 for plotting.
  3. Per-cell RNA: gene counts (shared counts layer of the scVI object) for the adhesion / EMT /
     endocrine panel; provirus = IRescue locus counts of the 4 copies (sense), per cell.
  4. Per-cell ATAC: Tn5 cut sites from the fragment files in the provirus 5' LTR, the CDH6 promoter,
     and every promoter (TSS +-1 kb) of protein-coding genes within +-2 Mb; normalised per bin by cut
     sites in a housekeeping-promoter background (no per-cell library size needed).
  5. Pseudotime bins (--n-bins, equal cell numbers): CPM-type means per bin; min-max scaled curves;
     half-max crossing point per feature (rising for provirus, falling for CDH6), with a bootstrap CI
     over cells.
  6. Cis links: metacells (k-means on the latent within each sample, ~--cells-per-metacell cells);
     Spearman of provirus RNA and of 5' LTR accessibility against every expressed gene within +-2 Mb,
     plus a null from --n-null random expressed genes genome-wide (empirical p).

Outputs ($TE_OUT/qc/herv9_cdh6/): cells.csv, trajectory_bins.csv, onset_summary.csv,
cis_links.csv, trajectory.png
Untested draft (2026-10-04): check the CHECK markers below (object paths / obsm key / IRescue
locus file layout) on the first run.
"""
import argparse
import glob
import gzip
import os
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import scipy.io
import scipy.sparse as sp
from scipy.stats import spearmanr

ap = argparse.ArgumentParser()
ap.add_argument("--with-late", action="store_true", help="also include late endocrine groups")
ap.add_argument("--max-per-group", type=int, default=4000)
ap.add_argument("--n-bins", type=int, default=40)
ap.add_argument("--cells-per-metacell", type=int, default=50)
ap.add_argument("--n-null", type=int, default=1000)
ap.add_argument("--n-boot", type=int, default=200)
ap.add_argument("--latent-key", default=os.environ.get("SCVI_LATENT_KEY", ""))
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()
rng = np.random.default_rng(args.seed)

REPO = Path(__file__).resolve().parents[1]
TE_OUT = Path(os.environ["TE_OUT"]); SCP_OUT = Path(os.environ["SCP_OUT"])
ANNOT = Path(os.environ["ANNOTATION_CSV"]); TE_ANNOT = Path(os.environ["TE_ANNOT_DIR"])
out = TE_OUT / "qc/herv9_cdh6"; out.mkdir(parents=True, exist_ok=True)
pd.set_option("display.width", 250)

LINEAGE = ["PP1", "PP2", "ENP1", "ENP2", "ENP3", "early_SC_beta", "early_EC_beta", "early_SC_alpha",
           "early_SC_delta", "early_SC_EC"]
if args.with_late:
    LINEAGE += ["late_SC_beta", "late_EC_beta", "late_SC_alpha", "late_SC_EC"]
PANEL = ["CDH6", "CDH1", "CDH2", "EPCAM", "CLDN4", "CLDN6", "CLDN7", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "VIM",
         "SOX9", "PDX1", "NKX6-1", "ONECUT1", "NEUROG3", "NEUROD1", "INSM1", "NKX2-2", "CHGA", "INS", "GCG", "SST",
         "TPH1", "DROSHA", "C5orf22"]
CHROM, PROV_START, PROV_END = "chr5", 30197738, 30205134
LTR5 = (30197638, 30198283)                                     # 5' LTR (LTR12B) +-100 bp
COPY_KEYS = ["chr5_30197738_30198183", "chr5_30198183_30203886", "chr5_30203884_30204661",
             "chr5_30204661_30205134"]
CIS = 2_000_000
HOUSEKEEPING = ["GAPDH", "ACTB", "TBP", "PGK1", "RPLP0", "HPRT1", "B2M", "PPIA", "YWHAZ", "SDHA", "UBC",
                "GUSB", "HMBS", "POLR2A", "EEF1A1", "TUBB", "RPL13A", "PSMB2", "VCP", "EIF4A2"]

# ---------------------------------------------------------------- 1. cells
ann = pd.read_csv(ANNOT, index_col=0, dtype=str)
ann = ann[(ann.include_downstream.str.lower() == "true") & ann.ATAC_WNN_annotation.isin(LINEAGE)]
ann = ann.groupby("ATAC_WNN_annotation", group_keys=False).apply(
    lambda d: d.sample(min(len(d), args.max_per_group), random_state=args.seed))
print("cells per group:\n", ann.ATAC_WNN_annotation.value_counts().reindex(LINEAGE).to_string(), flush=True)

# ---------------------------------------------------------------- 2. latent + pseudotime
SCVI_H5AD = SCP_OUT / "integration/h5ad/merged_rna_scvi_integrated.h5ad"     # CHECK: or scvi_filtered/


def find_latent(f):
    keys = list(f["obsm"].keys())
    if args.latent_key:
        return args.latent_key
    pref = [k for k in keys if "scvi" in k.lower() and "umap" not in k.lower()]
    if not pref:
        raise SystemExit(f"no scVI latent in obsm: {keys}; pass --latent-key")
    print("obsm keys:", keys, "-> using", pref[0], flush=True)
    return pref[0]


with h5py.File(SCVI_H5AD, "r") as f:
    cells_all = pd.Index(f["obs"]["_index"][:].astype(str))
    key = find_latent(f)
    sel = cells_all.get_indexer(ann.index)
    keep = sel >= 0
    ann = ann[keep]; sel = sel[keep]
    order = np.argsort(sel)
    lat = f["obsm"][key][np.sort(sel)]
    ann = ann.iloc[order]
    var = pd.Index(f["var"]["_index"][:].astype(str))
    # per-cell gene counts (raw counts layer, CSR) for the panel + all genes we may need later
    L = f["layers"]["counts"]; indptr = L["indptr"][:]
    rows = np.sort(sel)
    lib = np.zeros(len(rows), dtype=np.float64)
    data_blocks, ind_blocks, ptr = [], [], [0]
    for k, r in enumerate(rows):
        a, b = indptr[r], indptr[r + 1]
        d = L["data"][a:b]; i = L["indices"][a:b]
        lib[k] = d.sum(); data_blocks.append(d); ind_blocks.append(i); ptr.append(ptr[-1] + len(d))
    X = sp.csr_matrix((np.concatenate(data_blocks), np.concatenate(ind_blocks), np.array(ptr)),
                      shape=(len(rows), len(var)))
cells = ann.index
group = ann.ATAC_WNN_annotation.values
sample = np.array([c.split("#")[0] for c in cells])
print(f"{len(cells)} cells, latent {key} {lat.shape}", flush=True)

import anndata as ad
import scanpy as sc
A = ad.AnnData(X=np.zeros((len(cells), 1), dtype=np.float32), obs=pd.DataFrame({"group": group}, index=cells))
A.obsm["X_lat"] = np.asarray(lat, dtype=np.float32)
sc.pp.neighbors(A, use_rep="X_lat", n_neighbors=30)
sc.tl.diffmap(A, n_comps=10)
dc1 = A.obsm["X_diffmap"][:, 1]
pp_side = np.sign(np.median(dc1[group == "PP1"]) - np.median(dc1[np.isin(group, ["early_EC_beta", "early_SC_beta"])]))
cand = np.where(group == "PP1")[0]
A.uns["iroot"] = int(cand[np.argmax(pp_side * dc1[cand])])
sc.tl.dpt(A)
pt = A.obs["dpt_pseudotime"].values
print("median pseudotime per group:\n", pd.Series(pt).groupby(group).median().reindex(LINEAGE).round(3).to_string(), flush=True)

# ---------------------------------------------------------------- 3. provirus RNA per cell
# CHECK: assumes IRescue locus output per sample at $TE_OUT/irescue_locus/<sample>/counts/
#        {barcodes,features}.tsv.gz + matrix.mtx.gz (features x cells, 10x style)
prov = pd.Series(0.0, index=cells)
for s in np.unique(sample):
    d = TE_OUT / "irescue_locus" / s / "counts"
    if not d.exists():
        print("missing", d, flush=True); continue
    bc = pd.read_csv(d / "barcodes.tsv.gz", header=None)[0].astype(str).values
    feat = pd.read_csv(d / "features.tsv.gz", header=None, sep="\t")[0].astype(str).values
    fi = [k for k, x in enumerate(feat) if any(c in x for c in COPY_KEYS)]
    if not fi:
        print("provirus copies not found in", d, flush=True); continue
    m = scipy.io.mmread(d / "matrix.mtx.gz").tocsr()
    if m.shape[0] != len(feat):
        m = m.T.tocsr()
    v = np.asarray(m[fi].sum(0)).ravel()
    s_cells = pd.Series(v, index=[f"{s}#{b}" for b in bc])
    prov.update(s_cells.reindex(cells).dropna())

# ---------------------------------------------------------------- 4. genes: GTF, cis window, TSSs
rows_g = []
with open(TE_ANNOT / "scte/genes.gtf") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        t = line.split("\t", 9)
        if t[2] != "gene":
            continue
        a = t[8]
        name = a.split('gene_name "')[1].split('"')[0]
        gtype = a.split('gene_type "')[1].split('"')[0] if 'gene_type "' in a else ""
        tss = int(t[3]) - 1 if t[6] == "+" else int(t[4]) - 1
        rows_g.append((t[0], int(t[3]) - 1, int(t[4]), name, gtype, t[6], tss))
genes = pd.DataFrame(rows_g, columns=["chrom", "start", "end", "gene", "type", "strand", "tss"]).drop_duplicates("gene")
cis = genes[(genes.chrom == CHROM) & (genes.end > PROV_START - CIS) & (genes.start < PROV_END + CIS)]
print("genes within 2 Mb:", ", ".join(cis.gene), flush=True)

# ---------------------------------------------------------------- 5. ATAC cut sites per cell
windows = {"prov_5LTR": (CHROM, *LTR5)}
for g in cis.itertuples():
    windows[f"tss_{g.gene}"] = (g.chrom, g.tss - 1000, g.tss + 1000)
hk = genes[genes.gene.isin(HOUSEKEEPING)]
for g in hk.itertuples():
    windows[f"hk_{g.gene}"] = (g.chrom, g.tss - 1000, g.tss + 1000)
wnames = list(windows)
acc = np.zeros((len(cells), len(wnames)), dtype=np.float32)
cpos = {c: k for k, c in enumerate(cells)}
st = pd.read_csv(REPO / "configs/samples.tsv", sep="\t")
for s in st.itertuples():
    if s.sample not in set(sample):
        continue
    tb = pysam.TabixFile(str(Path(s.atac_bam).parent / "atac_fragments.tsv.gz"))     # as in script 21
    for j, w in enumerate(wnames):
        ch, a0, a1 = windows[w]
        try:
            it = tb.fetch(ch, max(0, a0 - 1000), a1 + 1000)
        except ValueError:
            continue
        for line in it:
            c, a, b, bc = line.split("\t")[:4]
            k = cpos.get(f"{s.sample}#{bc}")
            if k is None:
                continue
            for x in (int(a), int(b) - 1):
                if a0 <= x < a1:
                    acc[k, j] += 1
acc = pd.DataFrame(acc, index=cells, columns=wnames)
hk_cols = [w for w in wnames if w.startswith("hk_")]

# ---------------------------------------------------------------- per-cell table
panel = [g for g in PANEL if g in var]
gi = var.get_indexer(panel)
expr = pd.DataFrame(X[:, gi].toarray(), index=cells, columns=panel)
cell_tab = pd.DataFrame({"group": group, "sample": sample, "pseudotime": pt, "DC1": dc1,
                         "DC2": A.obsm["X_diffmap"][:, 2], "gene_umis": lib, "provirus": prov.values}, index=cells)
cell_tab = cell_tab.join(expr).join(acc[["prov_5LTR"] + [w for w in wnames if w == "tss_CDH6"]])
cell_tab["hk_cutsites"] = acc[hk_cols].sum(1).values
cell_tab.to_csv(out / "cells.csv")

# ---------------------------------------------------------------- 6. pseudotime bins + onset
cell_tab["bin"] = pd.qcut(cell_tab.pseudotime.rank(method="first"), args.n_bins, labels=False)


def bin_means(t):
    g = t.groupby("bin")
    res = pd.DataFrame({"pseudotime": g.pseudotime.median(), "n": g.size()})
    lib_b = g.gene_umis.sum()
    for c in ["provirus"] + panel:
        res[c] = g[c].sum() / lib_b * 1e6
    hk_b = g.hk_cutsites.sum().clip(lower=1)
    res["ATAC_prov_5LTR"] = g.prov_5LTR.sum() / hk_b * 1e3
    if "tss_CDH6" in t:
        res["ATAC_CDH6_promoter"] = g.tss_CDH6.sum() / hk_b * 1e3
    comp = pd.crosstab(t.bin, t.group, normalize="index")
    return res, comp


bins, comp = bin_means(cell_tab)
bins.join(comp.add_prefix("frac_")).round(4).to_csv(out / "trajectory_bins.csv")

FEATURES = {"provirus": "rise", "ATAC_prov_5LTR": "rise", "CDH6": "fall", "CDH1": "fall", "EPCAM": "fall",
            "SOX9": "fall", "NEUROG3": "rise", "NEUROD1": "rise", "CHGA": "rise", "SNAI2": "rise"}
if "ATAC_CDH6_promoter" in bins:
    FEATURES["ATAC_CDH6_promoter"] = "fall"


def halfmax(b, col, direction):
    y = b[col].rolling(3, center=True, min_periods=1).mean().values
    x = b.pseudotime.values
    lo_, hi_ = np.nanmin(y), np.nanmax(y)
    if hi_ - lo_ <= 0:
        return np.nan
    z = (y - lo_) / (hi_ - lo_)
    if direction == "rise":
        k = np.argmax(z >= 0.5)
    else:
        start = np.argmax(z)                        # fall after the peak
        after = np.where(z[start:] <= 0.5)[0]
        if not len(after):
            return np.nan
        k = start + after[0]
    return x[k]


feats = {c: d for c, d in FEATURES.items() if c in bins}
boot = {c: [] for c in feats}
for _ in range(args.n_boot):                      # resample cells, keep each cell's original bin
    t = cell_tab.sample(len(cell_tab), replace=True, random_state=int(rng.integers(1e9)))
    b, _ = bin_means(t)
    for c, d in feats.items():
        boot[c].append(halfmax(b, c, d))
onset = [{"feature": c, "direction": d, "halfmax_pt": halfmax(bins, c, d),
          "ci_lo": np.nanpercentile(boot[c], 2.5), "ci_hi": np.nanpercentile(boot[c], 97.5)}
         for c, d in feats.items()]
onset = pd.DataFrame(onset).sort_values("halfmax_pt")
onset.round(4).to_csv(out / "onset_summary.csv", index=False)
print("half-max crossing along pseudotime (rise = switch on, fall = switch off):\n", onset.round(3).to_string(index=False), flush=True)

# ---------------------------------------------------------------- 7. metacells + cis links
from sklearn.cluster import KMeans
mc = pd.Series("", index=cells)
for s in np.unique(sample):
    ix = np.where(sample == s)[0]
    k = max(1, len(ix) // args.cells_per_metacell)
    lab = KMeans(n_clusters=k, n_init=4, random_state=args.seed).fit_predict(np.asarray(lat)[ix]) if k > 1 else np.zeros(len(ix), int)
    mc.iloc[ix] = [f"{s}_m{x}" for x in lab]
M = sp.csr_matrix((np.ones(len(cells)), (pd.factorize(mc)[0], np.arange(len(cells)))))
mc_names = pd.factorize(mc)[1]
G = M @ X                                                     # metacell x gene counts
mlib = np.asarray(G.sum(1)).ravel()
expressed = np.where(np.asarray((G > 0).sum(0)).ravel() >= 0.2 * G.shape[0])[0]
mprov = M @ prov.values / mlib * 1e6
macc = (M @ acc.values)
mhk = macc[:, [wnames.index(w) for w in hk_cols]].sum(1).clip(min=1)
mltr = macc[:, wnames.index("prov_5LTR")] / mhk * 1e3


def gene_cpm(j):
    return np.asarray(G[:, j].todense()).ravel() / mlib * 1e6


null_j = rng.choice(expressed, size=min(args.n_null, len(expressed)), replace=False)
null_rna = np.array([spearmanr(mprov, gene_cpm(j))[0] for j in null_j])
null_acc = np.array([spearmanr(mltr, gene_cpm(j))[0] for j in null_j])
links = []
for g in cis.itertuples():
    if g.gene not in var:
        continue
    j = var.get_loc(g.gene)
    if j not in set(expressed):
        continue
    y = gene_cpm(j)
    r1 = spearmanr(mprov, y)[0]; r2 = spearmanr(mltr, y)[0]
    tsscol = f"tss_{g.gene}"
    r3 = spearmanr(mltr, macc[:, wnames.index(tsscol)] / mhk)[0] if tsscol in wnames else np.nan
    links.append({"gene": g.gene, "type": g.type, "dist_kb": round((g.tss - PROV_START) / 1000),
                  "rho_provRNA_gene": r1, "p_emp_rna": (np.sum(np.abs(null_rna) >= abs(r1)) + 1) / (len(null_rna) + 1),
                  "rho_5LTRacc_gene": r2, "p_emp_acc": (np.sum(np.abs(null_acc) >= abs(r2)) + 1) / (len(null_acc) + 1),
                  "rho_5LTRacc_promoterAcc": r3})
links = pd.DataFrame(links).sort_values("rho_provRNA_gene")
links.round(4).to_csv(out / "cis_links.csv", index=False)
print(f"\n{len(mc_names)} metacells; cis links (+-2 Mb), empirical p vs {len(null_j)} random genes:\n",
      links.round(3).to_string(index=False), flush=True)

# ---------------------------------------------------------------- 8. figure
fig, axs = plt.subplots(1, 3, figsize=(17, 4.6), gridspec_kw={"width_ratios": [1.2, 1.6, 1]})
ax = axs[0]
cmap = plt.get_cmap("tab20")
for k, g in enumerate(LINEAGE):
    m = group == g
    if m.any():
        ax.scatter(dc1[m], cell_tab.DC2.values[m], s=2, color=cmap(k % 20), label=g, rasterized=True)
ax.set_xlabel("DC1"); ax.set_ylabel("DC2"); ax.set_title("diffusion map (scVI latent)")
ax.legend(fontsize=6, markerscale=4, frameon=False, ncol=2)
ax = axs[1]
for col, c in [("provirus", "#b2182b"), ("ATAC_prov_5LTR", "#ef8a62"), ("CDH6", "#2166ac"),
               ("ATAC_CDH6_promoter", "#67a9cf"), ("CDH1", "#4d9221"), ("NEUROG3", "#762a83"), ("CHGA", "#999999")]:
    if col in bins:
        y = bins[col].rolling(3, center=True, min_periods=1).mean()
        ax.plot(bins.pseudotime, (y - y.min()) / (y.max() - y.min() + 1e-12), color=c, label=col,
                lw=2.2 if col in ("provirus", "CDH6") else 1.2)
ax.set_xlabel("diffusion pseudotime (PP1 root)"); ax.set_ylabel("scaled (bin means, min-max)")
ax.set_title("provirus vs CDH6 along PP -> ENP -> endocrine"); ax.legend(fontsize=7, frameon=False)
ax = axs[2]
if len(links):
    ax.barh(links.gene, links.rho_provRNA_gene, color=["#2166ac" if g == "CDH6" else "#bbbbbb" for g in links.gene])
    lo_, hi_ = np.percentile(null_rna, [2.5, 97.5])
    ax.axvspan(lo_, hi_, color="#fddbc7", alpha=0.6, label="95% random genes")
    ax.set_xlabel("Spearman, provirus RNA vs gene (metacells)"); ax.set_title("cis genes +-2 Mb")
    ax.legend(fontsize=7, frameon=False); ax.tick_params(axis="y", labelsize=7)
fig.tight_layout(); fig.savefig(out / "trajectory.png", dpi=120)
print("done", out)
