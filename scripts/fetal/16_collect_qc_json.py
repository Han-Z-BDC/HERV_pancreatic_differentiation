#!/usr/bin/env python
"""Collect every QC number for the 10 human fetal libraries into one JSON.

Sources, all read-only:
  cellranger/<s>/outs/metrics_summary.csv          mapping + calling metrics
  qc_filtered/qc_filter_summary.csv                step 1.1 cutoffs and costs (our 09)
  qc_filtered/per_cell_flags.csv.gz                per-cell metrics + donor label
  qc_filtered/remove_background/summary/*.json     SoupX soup composition + counts removed
  scp_fetal*/logs/*.out                            SoupX's printed rho (not in the JSON)
  qc_filtered/detect_doublets/summary/*.csv        step 1.3 as it ran (auto threshold)
  doublet_rule_comparison.csv                      the three candidate doublet rules
Histograms are binned here so the report stays small instead of shipping 104k rows.
"""
import json, glob, gzip, re, os
from pathlib import Path
import numpy as np
import pandas as pd
import h5py

BASE = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas")
REPO = Path("/pl/active/ZhuLab_BDC/GitHub/HERV_pancreatic_differentiation")
DS = {"olaniru2023": dict(kind="cells", label="Olaniru 2023", assay="scRNA-seq (whole cells)",
                          gse="GSE197064", inst="scp_fetal"),
      "migliorini2024": dict(kind="nuclei", label="Migliorini 2024", assay="snRNA-seq (nuclei)",
                             gse="GSE230403", inst="scp_fetal_migl")}
PUBLISHED = {"migliorini2024": 73274}

def num(x):
    if isinstance(x, str):
        x = x.replace(",", "").replace("%", "").strip().strip('"')
    try:
        return float(x)
    except Exception:
        return None

def cr_metrics(ds, sid):
    p = BASE / ds / "cellranger" / sid / "outs" / "metrics_summary.csv"
    if not p.exists():
        return {}
    d = pd.read_csv(p).iloc[0].to_dict()
    pick = {
        "called": "Estimated Number of Cells", "reads": "Number of Reads",
        "mean_reads_per_cell": "Mean Reads per Cell", "median_genes": "Median Genes per Cell",
        "median_umi": "Median UMI Counts per Cell", "saturation": "Sequencing Saturation",
        "valid_bc": "Valid Barcodes", "map_genome": "Reads Mapped Confidently to Genome",
        "map_txome": "Reads Mapped Confidently to Transcriptome",
        "intronic": "Reads Mapped Confidently to Intronic Regions",
        "exonic": "Reads Mapped Confidently to Exonic Regions",
        "intergenic": "Reads Mapped Confidently to Intergenic Regions",
        "antisense": "Reads Mapped Antisense to Gene", "in_cells": "Fraction Reads in Cells",
        "genes_detected": "Total Genes Detected",
    }
    return {k: num(d.get(v)) for k, v in pick.items() if v in d}

def rho_from_logs(inst):
    """SoupX prints 'Estimated global rho of X' per sample; logs are one file per array task."""
    out = {}
    for f in sorted(glob.glob(str(REPO / inst / "logs" / "*.out"))):
        txt = Path(f).read_text(errors="ignore")
        sid = None
        m = re.search(r"\[1\.2\][^\n]*?sample[: ]+(\S+)", txt) or re.search(r"=== (\S+) ===", txt)
        if m:
            sid = m.group(1)
        rhos = re.findall(r"Estimated global rho of ([\d.]+)", txt)
        if sid and rhos:
            out[sid] = float(rhos[-1])
    return out

samples, hists = [], {}
percell_all = []
for ds, meta in DS.items():
    qc = pd.read_csv(BASE / ds / "qc_filtered" / "qc_filter_summary.csv")
    pc = pd.read_csv(BASE / ds / "qc_filtered" / "per_cell_flags.csv.gz")
    pc["dataset"] = ds
    percell_all.append(pc)
    dbl = {}
    for f in glob.glob(str(BASE / ds / "qc_filtered" / "detect_doublets" / "summary" / "*.csv")):
        r = pd.read_csv(f).iloc[0].to_dict()
        dbl[r["sample_id"]] = r
    soup = {}
    for f in glob.glob(str(BASE / ds / "qc_filtered" / "remove_background" / "summary" / "*.json")):
        j = json.load(open(f))
        soup[j["sample_id"]] = j

    for _, q in qc.iterrows():
        sid = q["pool"]
        m = re.search(r"(\d+)PCW", sid)
        pcw = int(m.group(1)) if m else None
        rec = {"id": sid, "dataset": ds, "kind": meta["kind"], "pcw": pcw,
               "cr": cr_metrics(ds, sid), "qc": {k: (None if pd.isna(v) else v)
                                                 for k, v in q.to_dict().items()}}
        s = soup.get(sid)
        if s:
            h5 = BASE / ds / "qc_filtered" / "rna_qc" / "h5ad" / f"{sid}.qc_filtered.h5ad"
            tot = None
            if h5.exists():
                with h5py.File(h5, "r") as f:
                    tot = float(f["obs/total_counts"][:].sum())
            rec["soupx"] = {"counts_removed": s["total_counts_removed"],
                            "total_counts": tot,
                            "ambient_pct": (100 * s["total_counts_removed"] / tot) if tot else None,
                            "top_genes": s["top_10_genes_by_counts_removed"],
                            "tfidf_min": s.get("tfidf_min_used")}
        if sid in dbl:
            d = dbl[sid]
            rec["doublets"] = {k: (None if pd.isna(v) else v) for k, v in d.items()}
        samples.append(rec)

pc = pd.concat(percell_all, ignore_index=True)

# binned distributions, computed on the SAME population the cutoffs were derived from
METRICS = [("total_counts", True), ("n_genes", True), ("pct_mt", False), ("pct_ribo", False)]
for sid, g in pc.groupby("pool", sort=False):
    h = {}
    for col, logscale in METRICS:
        v = g[col].to_numpy(float)
        v = v[np.isfinite(v)]
        if logscale:
            v = np.log10(np.clip(v, 1, None))
            lo, hi = 0.0, 5.5
        else:
            lo, hi = 0.0, float(max(np.percentile(v, 99.5), 1e-3))
        cnt, edges = np.histogram(v, bins=60, range=(lo, hi))
        h[col] = {"counts": cnt.tolist(), "lo": lo, "hi": hi, "log10": logscale,
                  "median": float(np.median(g[col].to_numpy(float)))}
    hists[sid] = h

# rate-matched doublet comparison
cmp_path = BASE / "doublet_rule_comparison.csv"
dbl_cmp = pd.read_csv(cmp_path).to_dict("records") if cmp_path.exists() else []

# final composition: join doublet flags back to the donor labels
final = []
for ds in DS:
    h5d = BASE / ds / "qc_filtered" / "detect_doublets" / "h5ad"
    for p in sorted(h5d.glob("*.doublets_annotated.h5ad")):
        sid = p.name.split(".")[0]
        with h5py.File(p, "r") as f:
            bc = f["obs/_index"][:].astype(str)
            un = f["obs/doublet_filter_union"][:].astype(bool)
        sub = pc[(pc["pool"] == sid) & (pc["keep"])].set_index("barcode")
        lab = sub["geo_label"].reindex(bc).to_numpy()
        for L in pd.unique(lab[~pd.isna(lab)]):
            sel = (lab == L)
            final.append({"dataset": ds, "sample": sid, "donor": str(L),
                          "post_qc": int(sel.sum()), "post_doublet": int((sel & ~un).sum())})

out = {"generated": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"),
       "datasets": {k: v for k, v in DS.items()},
       "published": PUBLISHED,
       "rho_unclamped": {"pool1": 0.006, "pool2": 0.005, "pool3": 0.005, "pool4": 0.003},
       "samples": samples, "hists": hists, "doublet_rules": dbl_cmp, "final": final}
dest = BASE / "qc_review" / "fetal_qc_full.json"
dest.write_text(json.dumps(out))
print(f"wrote {dest}  ({dest.stat().st_size/1e6:.2f} MB)")
print(f"  {len(samples)} samples, {len(hists)} histogram sets, {len(final)} donor rows")
