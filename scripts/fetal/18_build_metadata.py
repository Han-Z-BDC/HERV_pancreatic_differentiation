#!/usr/bin/env python
"""Per-cell metadata for the fetal atlas: study, stage, donor, material.

Written from the authoritative sources only -- the GEO demultiplexing labels for Olaniru's
pooled libraries and the library identity for Migliorini's one-donor-per-library design --
so no field here is inferred from a file name except the gestational age, which is all either
study published.
"""
import re
from pathlib import Path
import numpy as np, pandas as pd, h5py

B = Path("/scratch/alpine/hzhu10@xsede.org/fetal_pancreas")
STUDY = {"olaniru2023": dict(gse="GSE197064", material="cell", assay="scRNA-seq"),
         "migliorini2024": dict(gse="GSE230403", material="nucleus", assay="snRNA-seq")}
# pool composition, from the GEO supplementary file names (GSM5907965_121313_... etc.)
POOL_AGES = {"pool1": "12,13,13", "pool2": "12,14,19", "pool3": "13,14,19", "pool4": "15,18,20"}
BINS = [(12, 14, "12-14 PCW"), (15, 16, "15-16 PCW"), (18, 20, "18-20 PCW")]


def stage_bin(p):
    for lo, hi, lab in BINS:
        if lo <= p <= hi:
            return lab
    return None


rows = []
for ds, meta in STUDY.items():
    pc = pd.read_csv(B / ds / "qc_filtered" / "per_cell_flags.csv.gz")
    pc = pc[pc["keep"]].copy()
    h5d = B / ds / "qc_filtered" / "detect_doublets" / "h5ad"
    for p in sorted(h5d.glob("*.doublets_annotated.h5ad")):
        lib = p.name.split(".")[0]
        with h5py.File(p, "r") as f:
            bc = f["obs/_index"][:].astype(str)
            dbl = f["obs/doublet_filter_union"][:].astype(bool)
        sub = pc[pc["pool"] == lib].set_index("barcode")
        lab = sub["geo_label"].reindex(bc).to_numpy()
        keep = ~dbl
        df = pd.DataFrame({"barcode": bc, "library": lib, "donor_label": lab,
                           "is_doublet": dbl, "in_atlas": keep})
        df["study"] = ds
        df["study_gse"] = meta["gse"]
        df["material"] = meta["material"]
        df["assay"] = meta["assay"]
        # Two separate facts, conflated in an earlier version of this script:
        #   multiplexed  -- more than one donor in the library. TRUE for all four Olaniru
        #                   pools, pool4 included: it carries three donors (15/18/20 PCW).
        #   hto_available -- usable hashtag FASTQ we could demultiplex ourselves. FALSE for
        #                   pool4: its deposited hashtag run shares 0% of barcodes with its
        #                   GEX run after 10x 3'v3 barcode translation, and GSM5907968's own
        #                   published matrix has zero Antibody Capture rows. Its donor labels
        #                   come from the authors' GEO demultiplexing table instead, which is
        #                   complete -- so pool4 has full PCW information but no independent
        #                   check on the assignment.
        df["multiplexed"] = ds == "olaniru2023"
        df["hto_available"] = ds == "olaniru2023" and lib != "pool4"
        df["donor_label_source"] = ("GEO demultiplexing table" if ds == "olaniru2023"
                                    else "one donor per library")
        if ds == "olaniru2023":
            df["pcw"] = df["donor_label"].str.extract(r"^(\d+)PCW")[0].astype(float)
        else:
            df["pcw"] = float(re.search(r"(\d+)PCW", lib).group(1))
            df["donor_label"] = lib
        # donor_id is LIBRARY-SCOPED on purpose: an age label shared by two Olaniru pools
        # cannot be shown from the deposited data to be the same embryo (see the schema doc).
        df["donor_id"] = df["study"] + "|" + df["library"] + "|" + df["donor_label"].astype(str)
        df["stage"] = "PCW" + df["pcw"].astype("Int64").astype(str)
        df["stage_bin"] = df["pcw"].map(lambda x: stage_bin(x) if pd.notna(x) else None)
        rows.append(df)

m = pd.concat(rows, ignore_index=True)
cols = ["barcode", "library", "study", "study_gse", "assay", "material", "multiplexed",
        "hto_available", "donor_label", "donor_label_source", "donor_id", "pcw", "stage",
        "stage_bin", "is_doublet", "in_atlas"]
m = m[cols]
out = B / "qc_review" / "fetal_cell_metadata.csv.gz"
m.to_csv(out, index=False)

atlas = m[m["in_atlas"]]
print(f"wrote {out}  ({len(m):,} QC-passing cells, {len(atlas):,} in atlas after doublets)\n")
print("=== per library ===")
t = atlas.groupby(["study", "library", "material", "multiplexed", "hto_available"],
                  sort=False).agg(
    cells=("barcode", "size"), donors=("donor_id", "nunique"),
    ages=("pcw", lambda s: ",".join(str(int(x)) for x in sorted(s.unique()))))
print(t.to_string())
print("\n=== per donor ===")
d = atlas.groupby(["study", "library", "donor_label", "pcw"], sort=False).size()
print(d.to_string())
print("\n=== per stage (both studies pooled) ===")
s = atlas.groupby(["stage", "material"], sort=False).size().unstack(fill_value=0)
s["total"] = s.sum(axis=1)
print(s.sort_index(key=lambda i: [int(x[3:]) for x in i]).to_string())
print("\n=== per stage_bin ===")
print(atlas.groupby(["stage_bin", "material"]).size().unstack(fill_value=0).to_string())
print(f"\ndistinct donor_id: {atlas['donor_id'].nunique()}   "
      f"distinct donor_label: {atlas['donor_label'].nunique()}")
