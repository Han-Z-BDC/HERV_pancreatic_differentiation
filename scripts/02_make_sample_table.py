#!/usr/bin/env python
"""
Sample table and per-sample cell whitelists for the TE counting jobs.

Reads the diff-multiome sample manifest (Cell Ranger folder -> analysis name) and the shared
annotation table; keeps cells with include_downstream == True.
Writes:
  configs/samples.tsv             task, sample, cellranger_folder, gex_bam, atac_bam, n_cells
  $TE_OUT/whitelists/<sample>.txt barcodes (as in the BAM CB tag, e.g. AAAC...-1), one per line
Usage: python scripts/02_make_sample_table.py   (paths from config.sh via environment)
"""
import os
from pathlib import Path

import pandas as pd

repo = Path(__file__).resolve().parents[1]
cfg = {}
for line in (repo / "config.sh").read_text().splitlines():
    if "=" in line and not line.lstrip().startswith(("#", "mkdir")):
        k, v = line.split("=", 1)
        cfg[k.strip()] = os.path.expandvars(v.split("#")[0].strip().strip('"').replace("$", "$"))
for k in list(cfg):  # resolve ${VAR} references between config entries
    for kk, vv in cfg.items():
        cfg[k] = cfg[k].replace(f"${kk}", vv)

man = pd.read_csv(cfg["SAMPLE_MANIFEST"]).query("status == 'ok'").sort_values("new_name")
ann = pd.read_csv(cfg["ANNOTATION_CSV"], index_col=0, dtype=str)
ann = ann[ann.include_downstream.str.lower() == "true"]
wl_dir = Path(cfg["TE_OUT"]) / "whitelists"; wl_dir.mkdir(parents=True, exist_ok=True)

rows = []
for i, r in enumerate(man.itertuples(), start=1):
    outs = next(Path(d) / r.old_name / "outs" for d in cfg["CELLRANGER_DIRS"].split()
                if (Path(d) / r.old_name / "outs").is_dir())
    bc = [b.split("#", 1)[1] for b in ann.index if b.startswith(r.new_name + "#")]
    (wl_dir / f"{r.new_name}.txt").write_text("\n".join(bc) + "\n")
    rows.append({"task": i, "sample": r.new_name, "cellranger_folder": r.old_name,
                 "gex_bam": outs / "gex_possorted_bam.bam", "atac_bam": outs / "atac_possorted_bam.bam",
                 "n_cells": len(bc)})
tab = pd.DataFrame(rows)
(repo / "configs").mkdir(exist_ok=True)
tab.to_csv(repo / "configs" / "samples.tsv", sep="\t", index=False)
print(tab[["task", "sample", "cellranger_folder", "n_cells"]].to_string(index=False))
print("total cells:", tab.n_cells.sum())
