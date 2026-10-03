#!/usr/bin/env python
"""
Tiered ERV candidate table: every RNA method on equal footing, plus read-through and ATAC evidence.

RNA calls come from 11_compare_scte_modes.py (IRescue = ir, scTE exclusive = ex, scTE nointron = no;
call = padj < 0.05 and log2FC > 1). For each subfamily x level (group or stage):
  candidate   up in IRescue OR scTE exclusive
  robust      the other method agrees in direction (log2FC > 0.5) and nointron is not intronic_only
  discordant  otherwise; `reasons` says why:
                read_through      nointron log2FC <= 0.5 (copies inside genes carry the signal)
                scte_only         IRescue flat; IRescue is sense-strand only, scTE unstranded ->
                                  candidate antisense transcription (or IRescue filtering)
                irescue_only      scTE exclusive flat; scTE drops reads on gene exons, IRescue keeps them
                absent_in_<m>     subfamily below the count filter of that method
ATAC: cut sites (05) and scTEATAC (05b) DESeq2 results, for the subfamily itself and for its LTRs.
  LTRs of an internal region are taken from the genome: RepeatMasker LTR-class elements within 50 bp
  of the -int fragments; an LTR subfamily is listed if it is >= 15% of those flanks (>= 20 flanks).
  atac_support = any of them up (padj < 0.05, log2FC > 1) in the same level by either ATAC method.
Small groups (< 200 cells in units) are flagged.

Outputs: $TE_OUT/qc/candidate_tiers_{group,stage}.csv, $TE_OUT/qc/erv_int_ltr_pairs.csv
"""
import os
import subprocess
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

TE_OUT = Path(os.environ["TE_OUT"]); D = TE_OUT / "deseq2"; qc = TE_OUT / "qc"
TE_ANNOT = Path(os.environ["TE_ANNOT_DIR"]); BEDTOOLS = Path(os.environ["TE_ENV"]) / "bin/bedtools"
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400); pd.set_option("display.max_colwidth", 60)


def int_ltr_pairs():
    """LTR subfamilies flanking each ERV internal subfamily in the genome."""
    rm = pd.read_csv(TE_ANNOT / "rmsk_hg38.bed.gz", sep="\t", header=None,
                     names=["chrom", "start", "end", "sub", "score", "strand", "family", "class"])
    rm = rm[rm["class"].str.rstrip("?") == "LTR"]
    is_int = rm["sub"].str.endswith("-int")
    tmp = qc / "tmp_pairs"; tmp.mkdir(exist_ok=True)
    rm[is_int].iloc[:, :6].to_csv(tmp / "int.bed", sep="\t", header=False, index=False)
    rm[~is_int].iloc[:, :6].to_csv(tmp / "ltr.bed", sep="\t", header=False, index=False)
    out = subprocess.run([str(BEDTOOLS), "closest", "-a", str(tmp / "int.bed"), "-b", str(tmp / "ltr.bed"),
                          "-k", "2", "-d", "-t", "all"], capture_output=True, text=True, check=True).stdout
    c = pd.read_csv(StringIO(out), sep="\t", header=None)
    c = c[(c[12] >= 0) & (c[12] <= 50)]
    n = c.groupby([3, 9]).size().rename("flanks").reset_index().rename(columns={3: "int", 9: "ltr"})
    n["share"] = n.flanks / n.groupby("int").flanks.transform("sum")
    pairs = n[(n.share >= 0.15) & (n.flanks >= 20)].sort_values(["int", "share"], ascending=[True, False])
    pairs.round(3).to_csv(qc / "erv_int_ltr_pairs.csv", index=False)
    for f in ["int.bed", "ltr.bed"]:
        (tmp / f).unlink()
    tmp.rmdir()
    return pairs.groupby("int").ltr.apply(list).to_dict()


pairs = int_ltr_pairs()
print("example pairings (genome-derived):", {k: pairs.get(k) for k in
      ["HERVH-int", "HERVS71-int", "HERV15-int", "HERVK-int", "HERVK11-int", "HERV9-int", "MER110-int",
       "MER92-int", "HERVIP10B3-int", "HERVIP10F-int", "PABL_A-int", "MLT1F2-int"]})


def atac_table(tag):
    a = {}
    for k, m in {"cut": "atac", "scte": "scte_atac"}.items():
        x = pd.read_csv(D / m / "results" / f"{tag}_vs_rest_all.csv").set_index(["subfamily", "level"])
        a[k] = x[["log2FoldChange", "padj"]]
    return a


def atac_evidence(sub, level, a):
    """Best ATAC evidence over the subfamily and its LTRs: (support, best_feature, lfc_cut, lfc_scte)."""
    feats = [sub] + pairs.get(sub, [])
    best = (False, "", np.nan, np.nan, -np.inf)
    for f in feats:
        lc = a["cut"].log2FoldChange.get((f, level), np.nan); pc = a["cut"].padj.get((f, level), np.nan)
        ls = a["scte"].log2FoldChange.get((f, level), np.nan); ps = a["scte"].padj.get((f, level), np.nan)
        up = (pc < 0.05 and lc > 1) or (ps < 0.05 and ls > 1)
        score = np.nanmax([lc, ls, -np.inf])
        if (up, score) > (best[0], best[4]):
            best = (up, f, lc, ls, score)
    return best[:4]


def tiers(tag):
    t = pd.read_csv(qc / f"scte_mode_comparison_{tag}.csv")
    t = t[(t["class"] == "LTR") & (t.up_ir | t.up_ex)].copy()
    agree = np.where(t.up_ir & t.up_ex, True,
                     np.where(t.up_ir, t.log2FoldChange_ex > 0.5, t.log2FoldChange_ir > 0.5))
    t["other_agrees"] = agree
    t["read_through"] = t.verdict.eq("intronic_only")
    t["tier"] = np.where(t.other_agrees & ~t.read_through, "robust", "discordant")

    def reasons(r):
        out = []
        if r.read_through:
            out.append("read_through")
        if not r.other_agrees:
            if r.up_ex and not r.up_ir:
                out.append("absent_in_irescue" if np.isnan(r.log2FoldChange_ir) else "scte_only")
            else:
                out.append("absent_in_scte" if np.isnan(r.log2FoldChange_ex) else "irescue_only")
        return ",".join(out)
    t["reasons"] = t.apply(reasons, axis=1)

    a = atac_table(tag)
    ev = [atac_evidence(s, l, a) for s, l in zip(t.subfamily, t.level)]
    t["atac_support"] = [e[0] for e in ev]
    t["atac_feature"] = [e[1] for e in ev]
    t["atac_lfc_cut"] = [e[2] for e in ev]
    t["atac_lfc_scte"] = [e[3] for e in ev]
    t["paired_ltrs"] = [";".join(pairs.get(s, [])) for s in t.subfamily]
    if tag == "group":
        cells = pd.read_csv(D / "rna/coldata.csv").groupby("group").n_cells.sum()
        t["group_cells"] = t.level.map(cells).values
        t["small_group"] = t.group_cells < 200
    cols = ["subfamily", "family", "level", "tier", "reasons", "baseMean_ir", "log2FoldChange_ir", "padj_ir",
            "log2FoldChange_ex", "padj_ex", "log2FoldChange_no", "padj_no", "verdict", "frac_ratio",
            "atac_support", "atac_feature", "atac_lfc_cut", "atac_lfc_scte", "paired_ltrs"] + \
           (["group_cells", "small_group"] if tag == "group" else [])
    t = t[cols].sort_values(["level", "tier", "log2FoldChange_ir"], ascending=[True, False, False])
    t.round(3).to_csv(qc / f"candidate_tiers_{tag}.csv", index=False)

    print(f"\n===== {tag}: {len(t)} ERV/LTR candidates; robust {(t.tier == 'robust').sum()} "
          f"({t[t.tier == 'robust'].atac_support.sum()} with ATAC support); discordant {(t.tier == 'discordant').sum()}")
    print("discordant reasons:", t[t.tier == "discordant"].reasons.value_counts().to_dict())
    show = ["subfamily", "level", "tier", "reasons", "log2FoldChange_ir", "log2FoldChange_ex", "log2FoldChange_no",
            "atac_support", "atac_feature", "atac_lfc_cut", "atac_lfc_scte"]
    big = t[(t.baseMean_ir.fillna(0) > 20) | (t.log2FoldChange_ex > 1.5)]
    print(big[show].round(2).to_string(index=False), flush=True)
    return t


tiers("group")
tiers("stage")
print("done", qc)
