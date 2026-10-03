#!/usr/bin/env python
"""
Name the copies behind each ERV candidate (14_candidate_tiers.py) from the IRescue locus counts (15).

1. Pseudobulk the per-cell locus counts into the 119 RNA units (sample x group, as 07):
   $TE_OUT/irescue_locus/locus_units.npz (+ _loci.txt, _units.csv). Expression = UMIs per million gene UMIs.
2. Host-gene context of every LTR-class copy, from the Cell Ranger ARC GTF (same genes as the RNA):
   intergenic / intronic / exonic, the gene name, and whether the copy is sense or antisense to it.
3. For each candidate subfamily x group: delta_i = CPM in the group - mean CPM over the other groups,
   per copy i. Copies are ranked by delta; reported are the copies that make up 80% of the summed
   positive delta, with their context and how many of the group's units show them (> 2x rest mean).
   n50 = number of copies giving 50% of the increase (1-2 = driven by single loci).

Outputs: $TE_OUT/qc/locus_candidates_{group}.csv (one row per candidate x copy),
         $TE_OUT/qc/locus_candidate_summary.csv (one row per candidate),
         $TE_OUT/qc/ltr_loci_context.csv.gz
"""
import os
import subprocess
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[1]
TE_OUT = Path(os.environ["TE_OUT"]); ANNOT = Path(os.environ["ANNOTATION_CSV"]); SCP_OUT = Path(os.environ["SCP_OUT"])
TE_ANNOT = Path(os.environ["TE_ANNOT_DIR"]); BEDTOOLS = str(Path(os.environ["TE_ENV"]) / "bin/bedtools")
qc = TE_OUT / "qc"; LOC = TE_OUT / "irescue_locus"
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)

# ---------- 1. pseudobulk locus counts into the RNA units
units_ref = pd.read_csv(TE_OUT / "deseq2/rna/coldata.csv", index_col=0)
units = units_ref.index
ann = pd.read_csv(ANNOT, index_col=0, dtype=str)
ann = ann[ann.include_downstream.str.lower() == "true"]
npz = LOC / "locus_units.npz"
if npz.exists():
    M = sp.load_npz(npz).tocsr(); loci = pd.Index(Path(LOC / "locus_units_loci.txt").read_text().split("\n"))
else:
    samples = pd.read_csv(REPO / "configs/samples.tsv", sep="\t")["sample"].tolist()
    loci = pd.Index(pd.read_csv(TE_ANNOT / "irescue/ltr_loci.bed", sep="\t", header=None, usecols=[3])[3])
    M = sp.csr_matrix((len(units), len(loci)), dtype=np.float64)
    for s in samples:
        d = LOC / s / "counts"
        m = scipy.io.mmread(d / "matrix.mtx.gz").tocsr().T          # cells x features
        f = pd.read_csv(d / "features.tsv.gz", sep="\t", header=None)[0].values
        b = s + "#" + pd.read_csv(d / "barcodes.tsv.gz", header=None)[0].astype(str)
        grp = ann["ATAC_WNN_annotation"].reindex(b).astype(str).values
        codes = units.get_indexer(s + "__" + pd.Series(grp))
        keep = codes >= 0
        ind = sp.csr_matrix((np.ones(keep.sum()), (codes[keep], np.where(keep)[0])), shape=(len(units), m.shape[0]))
        part = (ind @ m).tocoo()
        M = M + sp.csr_matrix((part.data, (part.row, loci.get_indexer(f)[part.col])), shape=M.shape)
        print(f"{s}: {m.shape[0]} cells, {m.shape[1]} expressed loci, {m.sum():,.0f} UMIs", flush=True)
    sp.save_npz(npz, M.tocsr()); Path(LOC / "locus_units_loci.txt").write_text("\n".join(loci))
    units_ref.to_csv(LOC / "locus_units_units.csv")
cpm = sp.diags(1e6 / units_ref.lib_size.values) @ M                 # units x loci, per million gene UMIs
groups = units_ref.group.values
G = pd.Index(sorted(set(groups)))
gi = sp.csr_matrix((np.ones(len(units)) / pd.Series(groups).map(pd.Series(groups).value_counts()).values,
                    (G.get_indexer(groups), np.arange(len(units)))), shape=(len(G), len(units)))
gcpm = (gi @ cpm).tocsr()                                            # group mean of unit CPMs, groups x loci

# ---------- 2. host-gene context of each copy
ctx_file = qc / "ltr_loci_context.csv.gz"
if ctx_file.exists():
    ctx = pd.read_csv(ctx_file, index_col=0)
else:
    gtf = TE_ANNOT / "scte/genes.gtf"
    rows = []
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if t[2] not in ("gene", "exon"):
                continue
            name = t[8].split('gene_name "')[1].split('"')[0]
            rows.append((t[0], int(t[3]) - 1, int(t[4]), name, t[2], t[6]))
    g = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene", "kind", "strand"])
    tmp = qc / "tmp_ctx"; tmp.mkdir(exist_ok=True)
    for k in ("gene", "exon"):
        g[g.kind == k].sort_values(["chrom", "start"])[["chrom", "start", "end", "gene", "kind", "strand"]] \
            .to_csv(tmp / f"{k}.bed", sep="\t", header=False, index=False)
    lb = TE_ANNOT / "irescue/ltr_loci.bed"

    def hits(k):
        out = subprocess.run([BEDTOOLS, "intersect", "-a", str(lb), "-b", str(tmp / f"{k}.bed"), "-wa", "-wb"],
                             capture_output=True, text=True, check=True).stdout
        h = pd.read_csv(StringIO(out), sep="\t", header=None, usecols=[3, 5, 9, 11], names=["locus", "te_strand", "gene", "g_strand"])
        h["orient"] = np.where(h.te_strand == h.g_strand, "sense", "antisense")
        return h.groupby("locus").agg(gene=("gene", lambda x: ",".join(sorted(set(x)))),
                                      orient=("orient", lambda x: ",".join(sorted(set(x)))))
    ex, ge = hits("exon"), hits("gene")
    for f in ("gene.bed", "exon.bed"):
        (tmp / f).unlink()
    tmp.rmdir()
    ctx = pd.DataFrame(index=loci)
    ctx["context"] = np.where(ctx.index.isin(ex.index), "exonic", np.where(ctx.index.isin(ge.index), "intronic", "intergenic"))
    ctx["host_gene"] = ge.gene.reindex(ctx.index)
    ctx["orientation"] = ge.orient.reindex(ctx.index)
    ctx.to_csv(ctx_file)
print("copy context:", ctx.context.value_counts().to_dict(), flush=True)

# ---------- 3. copies behind each candidate
cand = pd.read_csv(qc / "candidate_tiers_group.csv")
sub_of = pd.Series([l.split("|")[0] for l in loci], index=loci)
by_sub = pd.Series(np.arange(len(loci))).groupby(sub_of.values).apply(list)
rows, summ = [], []
for r in cand.itertuples():
    if r.subfamily not in by_sub.index:
        continue
    idx = np.array(by_sub[r.subfamily]); g = G.get_loc(r.level)
    gc = gcpm[:, idx].toarray()
    in_g = gc[g]; rest = np.delete(gc, g, axis=0).mean(0)
    delta = in_g - rest
    pos = np.where(delta > 0)[0]
    tot = delta[pos].sum()
    if tot <= 0:
        summ.append((r.subfamily, r.level, r.tier, r.reasons, 0, np.nan, np.nan, "", ""))
        continue
    order = pos[np.argsort(-delta[pos])]
    cum = np.cumsum(delta[order]) / tot
    n50 = int(np.searchsorted(cum, 0.5) + 1); n80 = int(np.searchsorted(cum, 0.8) + 1)
    u = np.where(groups == r.level)[0]
    ucpm = cpm[u][:, idx].toarray()
    for k, j in enumerate(order[:min(n80, 25)]):
        lid = loci[idx[j]]
        rows.append({"subfamily": r.subfamily, "level": r.level, "tier": r.tier, "reasons": r.reasons, "rank": k + 1,
                     "locus": lid.split("|")[1], "cpm_group": in_g[j], "cpm_rest": rest[j],
                     "log2fc": np.log2((in_g[j] + 0.5) / (rest[j] + 0.5)), "share_of_increase": delta[j] / tot,
                     "units_up": f"{int((ucpm[:, j] > 2 * rest[j]).sum())}/{len(u)}",
                     "context": ctx.context.get(lid), "host_gene": ctx.host_gene.get(lid), "orientation": ctx.orientation.get(lid)})
    top = [x for x in rows if x["subfamily"] == r.subfamily and x["level"] == r.level][:n50]
    ctx_share = pd.Series({c: sum(x["share_of_increase"] for x in rows if x["subfamily"] == r.subfamily
                                  and x["level"] == r.level and x["context"] == c) for c in ["intergenic", "intronic", "exonic"]})
    summ.append((r.subfamily, r.level, r.tier, r.reasons, len(pos), n50, n80,
                 "; ".join(f"{x['locus']} ({x['context']}{', ' + str(x['host_gene']) if x['host_gene'] == x['host_gene'] and x['host_gene'] else ''})" for x in top[:3]),
                 ", ".join(f"{c} {v:.0%}" for c, v in ctx_share.items() if v > 0)))
loc = pd.DataFrame(rows)
loc.round(4).to_csv(qc / "locus_candidates_group.csv", index=False)
s = pd.DataFrame(summ, columns=["subfamily", "level", "tier", "reasons", "copies_up", "n50", "n80", "top_copies", "context_of_top80"])
s.to_csv(qc / "locus_candidate_summary.csv", index=False)
print(s[s.tier == "robust"].to_string(index=False))
print("\ndiscordant:\n", s[s.tier != "robust"].to_string(index=False))
print("done", qc)
