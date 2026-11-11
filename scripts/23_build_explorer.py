#!/usr/bin/env python3
"""Build docs/HERV_explorer.html -- an interactive explorer over the FULL TE/ERV result set.

Unlike docs/HERV_progress_report.html (a narrative report on selected candidates), this emits
every subfamily x group/stage statistic from all four quantifications, plus the per-locus
breakdown for every called candidate, and lets the reader drill in themselves.

Inputs (all under results/):
  scte_modes/scte_mode_comparison_{group,stage}.csv   RNA: IRescue + scTE exclusive + scTE nointron
  atac_methods/atac_method_comparison_{group,stage}.csv  ATAC: Tn5 cut sites + scTEATAC
  atac_methods/atac_method_totals.csv                 per-subfamily ATAC library totals
  candidates/candidate_tiers_{group,stage}.csv        tier calls
  candidates/erv_int_ltr_pairs.csv                    genome-derived int<->LTR pairing
  qc/irescue_cpm_by_{group,stage}.csv                 absolute CPM (contrast-free trend)
  scte_modes/scte_intergenic_fraction_by_group.csv    intergenic read share
  locus/locus_candidates_group.csv, locus_followup.csv, locus_atac.csv,
  locus/locus_candidate_summary.csv                   per-copy results

Usage:  python3 scripts/23_build_explorer.py
"""
import csv, json, os, sys
from collections import defaultdict
from html.parser import HTMLParser

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R    = os.path.join(ROOT, "results")
TPL  = os.path.join(ROOT, "docs", "explorer_template.html")
OUT  = os.path.join(ROOT, "docs", "HERV_explorer.html")

# developmental order, taken from the CPM table header (biological order, not alphabetical)
GROUPS = ["hPSC","DE","GT","OTX2_endo","FG","liver","PP1","PP2","ductal","exocrine",
          "ENP1","ENP2","ENP3","early_SC_beta","early_EC_beta","early_SC_alpha","early_SC_delta",
          "early_SC_EC","late_SC_beta","late_EC_beta","late_SC_alpha","late_SC_delta","late_SC_EC","stroma"]
STAGES = ["D0","D4","D7","D9","D12","D15","D22","D45"]
LV = {"g": GROUPS, "s": STAGES}

def rd(p):
    f = os.path.join(R, p)
    if not os.path.exists(f):
        print("  ! missing", p, file=sys.stderr); return []
    with open(f, newline="") as fh:
        return list(csv.DictReader(fh))

def num(v, nd=None):
    """float or None; optionally rounded. Keeps 0.0 but drops '', NA, nan, inf."""
    if v is None: return None
    s = str(v).strip()
    if s in ("", "NA", "nan", "NaN", "None", "inf", "-inf", "Inf"): return None
    try: x = float(s)
    except ValueError: return None
    if x != x or x in (float("inf"), float("-inf")): return None
    if nd is not None:
        x = round(x, nd)
        if x == int(x) and abs(x) < 1e15: x = int(x) if nd == 0 else x
    return x

def sig2(v):
    """padj rounded to 2 significant digits (keeps JSON short, keeps readability)."""
    x = num(v)
    if x is None: return None
    if x == 0: return 0
    from math import floor, log10
    try: e = floor(log10(abs(x)))
    except ValueError: return None
    return round(x, max(0, -e + 1))

def boolish(v): return str(v).strip().lower() in ("true", "1", "yes")

print("reading result tables ...")

# ---------- 1. RNA: 3 quantifications ----------
rna = {"g": {}, "s": {}}
meta = {}
verdicts = set()
for lt, fn in (("g", "scte_modes/scte_mode_comparison_group.csv"),
               ("s", "scte_modes/scte_mode_comparison_stage.csv")):
    idx = {n: i for i, n in enumerate(LV[lt])}
    for row in rd(fn):
        sf, lvl = row["subfamily"], row["level"]
        if lvl not in idx: continue
        meta.setdefault(sf, {})
        if row.get("family"): meta[sf]["f"] = row["family"]
        if row.get("class"):  meta[sf]["c"] = row["class"]
        arr = rna[lt].setdefault(sf, [None] * len(idx))
        verdicts.add(row.get("verdict", ""))
        arr[idx[lvl]] = [
            num(row["log2FoldChange_ir"], 2), sig2(row["padj_ir"]),
            num(row["log2FoldChange_ex"], 2), sig2(row["padj_ex"]),
            num(row["log2FoldChange_no"], 2), sig2(row["padj_no"]),
            num(row.get("frac_ratio"), 2), (row.get("verdict") or None),
        ]
    # fill untested levels with an all-null record so the JS can index blindly
    for sf, arr in rna[lt].items():
        for i, v in enumerate(arr):
            if v is None: arr[i] = [None] * 8
print(f"  RNA  group={len(rna['g'])} stage={len(rna['s'])} subfamilies; verdicts={sorted(v for v in verdicts if v)}")

# ---------- 2. ATAC: 2 quantifications ----------
atac = {"g": {}, "s": {}}
agree = set()
for lt, fn in (("g", "atac_methods/atac_method_comparison_group.csv"),
               ("s", "atac_methods/atac_method_comparison_stage.csv")):
    idx = {n: i for i, n in enumerate(LV[lt])}
    for row in rd(fn):
        sf, lvl = row["subfamily"], row["level"]
        if lvl not in idx: continue
        meta.setdefault(sf, {})
        if row.get("family"): meta[sf].setdefault("f", row["family"])
        if row.get("class"):  meta[sf].setdefault("c", row["class"])
        arr = atac[lt].setdefault(sf, [None] * len(idx))
        agree.add(row.get("agreement", ""))
        arr[idx[lvl]] = [
            num(row["log2FoldChange_cut"], 2), sig2(row["padj_cut"]),
            num(row["log2FoldChange_scte"], 2), sig2(row["padj_scte"]),
            (row.get("agreement") or None),
        ]
    for sf, arr in atac[lt].items():
        for i, v in enumerate(arr):
            if v is None: arr[i] = [None] * 5
print(f"  ATAC group={len(atac['g'])} stage={len(atac['s'])} subfamilies; agreement={sorted(a for a in agree if a)}")

# ---------- 3. ATAC library totals ----------
for row in rd("atac_methods/atac_method_totals.csv"):
    sf = row.get("") or row.get("subfamily")
    if not sf: continue
    m = meta.setdefault(sf, {})
    m["cut"] = num(row.get("cut_sites"), 0)
    m["sct"] = num(row.get("scte_fragments"), 0)
    m["rt"]  = num(row.get("ratio"), 2)
    if row.get("family"): m.setdefault("f", row["family"])
    if row.get("class"):  m.setdefault("c", row["class"])

# ---------- 4. tier calls ----------
tiers = {"g": defaultdict(dict), "s": defaultdict(dict)}
ncand = 0
for lt, fn in (("g", "candidates/candidate_tiers_group.csv"),
               ("s", "candidates/candidate_tiers_stage.csv")):
    for row in rd(fn):
        sf = row["subfamily"]
        tiers[lt][sf][row["level"]] = {
            "tier": row["tier"], "rs": row.get("reasons") or None,
            "ir": num(row["log2FoldChange_ir"], 2), "pir": sig2(row["padj_ir"]),
            "ex": num(row["log2FoldChange_ex"], 2), "pex": sig2(row["padj_ex"]),
            "no": num(row["log2FoldChange_no"], 2), "pno": sig2(row["padj_no"]),
            "vd": row.get("verdict") or None, "fr": num(row.get("frac_ratio"), 2),
            "at": 1 if boolish(row.get("atac_support")) else 0,
            "af": row.get("atac_feature") or None,
            "ac": num(row.get("atac_lfc_cut"), 2), "as": num(row.get("atac_lfc_scte"), 2),
            "pl": row.get("paired_ltrs") or None,
            "gc": num(row.get("group_cells"), 0),
            "sg": 1 if boolish(row.get("small_group")) else 0,
        }
        ncand += 1
print(f"  tiers: {ncand} candidate calls over "
      f"{len(set(tiers['g']) | set(tiers['s']))} subfamilies")

# ---------- 5. CPM + intergenic fraction (contrast-free trends) ----------
def wide(path, order):
    out = {}
    for row in rd(path):
        sf = row.get("") or row.get("subfamily")
        if not sf: continue
        out[sf] = [num(row.get(c), 2) for c in order]
    return out
cpm = {"g": wide("qc/irescue_cpm_by_group.csv", GROUPS),
       "s": wide("qc/irescue_cpm_by_stage.csv", STAGES)}
igf = wide("scte_modes/scte_intergenic_fraction_by_group.csv", GROUPS)
print(f"  CPM group={len(cpm['g'])} stage={len(cpm['s'])}; intergenic fraction={len(igf)}")

# ---------- 6. per-locus ----------
# base table (all candidate loci) + followup (segdup / host expr) + ATAC per locus
loci = defaultdict(list)
extra, latac = {}, {}
for row in rd("locus/locus_followup.csv"):
    extra[(row["subfamily"], row["level"], row["locus"])] = row
for row in rd("locus/locus_atac.csv"):
    latac[(row["subfamily"], row["level"], row["locus"])] = row

LT_OF = {}
for lt in ("g", "s"):
    for n in LV[lt]: LT_OF[n] = lt

for row in rd("locus/locus_candidates_group.csv"):
    sf, lvl, lo = row["subfamily"], row["level"], row["locus"]
    lt = LT_OF.get(lvl)
    if lt is None: continue
    ex = extra.get((sf, lvl, lo), {})
    at = latac.get((sf, lvl, lo), {})
    loci[f"{sf}|{lvl}|{lt}"].append({
        "rk": num(row.get("rank"), 0), "lo": lo,
        "cg": num(row.get("cpm_group"), 2), "cr": num(row.get("cpm_rest"), 2),
        "fc": num(row.get("log2fc"), 2), "sh": num(row.get("share_of_increase"), 3),
        "uu": row.get("units_up") or None, "cx": row.get("context") or None,
        "hg": row.get("host_gene") or None, "or": row.get("orientation") or None,
        "sd": num(ex.get("segdup_identity"), 4),
        "he": ex.get("host_expr") or None, "np": ex.get("nearest_pc") or None,
        "ne": ex.get("nearest_expr") or None,
        "al": num(at.get("atac_log2fc"), 2), "au": at.get("atac_units_open") or None,
        "wl": at.get("window_ltrs") or None,
    })
for k in loci: loci[k].sort(key=lambda r: (r["rk"] if r["rk"] is not None else 1e9))

lsum = {}
for row in rd("locus/locus_candidate_summary.csv"):
    lt = LT_OF.get(row["level"])
    if lt is None: continue
    lsum[f'{row["subfamily"]}|{row["level"]}|{lt}'] = {
        "cu": num(row.get("copies_up"), 0), "n50": num(row.get("n50"), 0),
        "n80": num(row.get("n80"), 0), "ctx": row.get("context_of_top80") or None,
    }
print(f"  loci: {sum(len(v) for v in loci.values())} rows over {len(loci)} subfamily x level; "
      f"{len(lsum)} summaries")

# ---------- 7. int <-> LTR pairs ----------
pairs = defaultdict(list)
for row in rd("candidates/erv_int_ltr_pairs.csv"):
    pairs[row["int"]].append([row["ltr"], num(row.get("flanks"), 0), num(row.get("share"), 3)])
    pairs[row["ltr"]].append([row["int"], num(row.get("flanks"), 0), num(row.get("share"), 3)])
print(f"  pairs: {len(pairs)} subfamilies with a genome-derived partner")

# ---------- 8. family universe ----------
fam = sorted(set(rna["g"]) | set(rna["s"]) | set(atac["g"]) | set(atac["s"]))
meta = {k: v for k, v in meta.items() if k in fam}
print(f"  universe: {len(fam)} subfamilies")

HELP = """
<div class="prose">
<h2 style="font-size:19px;color:var(--fg);margin-top:0">How to read this explorer</h2>
<p>This page is the <b>complete</b> result set, not a shortlist. Every transposable-element
subfamily that passed filtering is here, with its statistics in all four independent
quantifications, across all 24 annotated cell groups and all 8 differentiation timepoints,
plus the per-genomic-copy breakdown wherever that was computed. Nothing is pre-selected
for you; the <span class="kbd">All results</span> view is the raw table.</p>

<h3>The four quantifications, and why there are four</h3>
<p>A TE "expression" number is only as good as the read-assignment rule behind it, and every
rule has a failure mode. So the same cells were counted four separate ways and the answers
are shown side by side. Agreement across methods is the actual evidence; a hit in one method
alone is a hypothesis about that method.</p>
<ul>
<li><b>IRescue (RNA)</b> — UMI-aware expectation-maximisation. Reads that map ambiguously to
many copies of the same subfamily get probabilistically assigned rather than discarded, so
this is the most sensitive RNA method. Run with <code>--strandedness forward</code>, so a read
only counts for a TE on the TE's own strand. Downside: EM can smear signal between related
subfamilies.</li>
<li><b>scTE exclusive</b> — counts a TE read only if it does not overlap a gene <i>exon</i>.
Conservative about mature mRNA contamination, but still counts intronic reads, so an ERV
sitting in an intron of a highly induced gene will light up from host pre-mRNA alone.</li>
<li><b>scTE nointron</b> — counts a TE read only if it falls outside gene bodies
<i>entirely</i>. This is the read-through control. If a subfamily is still up here, the
signal cannot be host pre-mRNA. If it collapses here but was up in the other two, the
"expression" was the host gene all along.</li>
<li><b>Tn5 cut sites (ATAC)</b> — counts only the two base pairs each fragment actually cut
(fragment start and end−1) falling inside a TE. This is the strict definition of accessibility
at the element: a fragment that merely spans a TE contributes nothing.</li>
<li><b>scTEATAC</b> — counts whole fragment spans overlapping a TE. More sensitive, lower
resolution, and prone to crediting a TE for a neighbouring open region. Shown alongside the
cut-site counts; the <code>agreement</code> column records whether they concur.</li>
</ul>
<p class="hint">The scTEATAC run needed a fix before it produced anything usable: the tool
reads the cell barcode from the read <i>name</i>, not the <code>CB</code> tag, and needs a
name-sorted BAM. Without both, all cells of a sample collapse into one or two pseudo-cells.
See <code>scripts/05b_scte_atac.sb</code>.</p>

<h3>Statistics</h3>
<p>All contrasts are pseudobulk DESeq2 on <code>sample &times; cell-group</code> units with
design <code>~ 0 + group</code>, testing each level against the mean of all others
(contrast +1 for the level, −1/(k−1) for each other). Size factors are <b>fixed to library
size</b> — total gene UMIs for RNA, total cut sites for ATAC — rather than estimated from
the TE matrix itself, because a global shift in TE load is a real biological signal here and
median-ratio normalisation would absorb it. There is deliberately no sample or batch term:
in this design stage, batch and cell type are confounded, so a batch term would remove the
effect of interest. <code>padj</code> is Benjamini–Hochberg across subfamilies within a
contrast. A <code>padj</code> of exactly <code>0</code> is underflow, not a perfect result.</p>

<h3>The read-through test</h3>
<p>Two columns encode it. <b>verdict</b> compares a subfamily's log2FC in scTE exclusive
against scTE nointron:</p>
<ul>
<li><span class="pill v-intergenic">intergenic</span> — still clearly up with gene bodies
excluded. The strongest class.</li>
<li><span class="pill v-weak">weak</span> — survives nointron but only just
(log2FC &gt; 0.5). <b>Check the nointron <code>padj</code> yourself</b>: this label does not
require significance, so some "weak" calls are not significant at all.</li>
<li><span class="pill v-intronic_only">intronic_only</span> — gone once gene bodies are
removed. Host-transcript read-through.</li>
<li><span class="pill v-no_intergenic">no_intergenic</span> — the subfamily has too little
intergenic signal to test.</li>
</ul>
<p><b>frac ratio</b> is the share of that subfamily's reads falling intergenically in the
called group, divided by the same share in all other groups. Above 1 means the element
becomes <i>more</i> independent of genes in this group, which is what an autonomously
transcribed ERV should do. Below 1 means the opposite.</p>
<p class="hint"><b>A known blind spot:</b> the nointron filter drops any read inside a gene
body regardless of strand. An ERV that sits <i>antisense</i> within a host intron therefore
fails the test even though host pre-mRNA could not produce its reads, because IRescue is
strand-aware. Check the <code>orientation</code> column in the per-locus tab before accepting a
read-through verdict.</p>

<h3>Tiers</h3>
<p>Both RNA methods are treated as equal footing — neither is the "main" one — because
picking either as primary silently discards real candidates (27 lost if IRescue leads,
43 if scTE leads).</p>
<ul>
<li><b>candidate</b>: significantly up in IRescue <i>or</i> in scTE exclusive.</li>
<li><b>robust</b>: a candidate where the other RNA method also has log2FC &gt; 0.5 and the
read-through verdict is not <code>intronic_only</code>.</li>
<li><b>discordant</b>: a candidate failing that, with the reason recorded —
<code>read_through</code>, <code>irescue_only</code>, <code>scte_only</code>,
<code>absent_in_irescue</code>.</li>
</ul>
<p>Discordant is not a synonym for false. <code>scte_only</code> with a strong nointron
log2FC, for instance, usually means IRescue's EM distributed the reads to a sibling subfamily.</p>

<h3>Per-locus tab</h3>
<p>Subfamily-level calls are averages over hundreds to thousands of genomic copies, and the
average can be produced by a single locus. Locus-mode IRescue (<code>-l</code>) was run over
725,763 LTR-class copies to resolve this.</p>
<ul>
<li><b>share of increase</b> — fraction of the subfamily's CPM increase in that group
contributed by this one copy. A copy above ~40% means the subfamily call is really a
single-locus event, and should be described that way.</li>
<li><b>n50 / n80</b> — how many copies account for 50% / 80% of the increase. <code>n50 = 1</code>
is one locus wearing a subfamily's name.</li>
<li><b>context / host gene / orientation</b> — whether the copy is intergenic, intronic or
exonic, in which gene, and sense or antisense to it.</li>
<li><b>segdup</b> — identity to a segmental duplication. Values above 0.98 are flagged: in
those regions multi-mapping is effectively unresolvable and the locus assignment should not
be trusted.</li>
<li><b>LTRs in window</b> — other LTR-class elements in the surrounding window, which is how
a provirus reveals itself (internal body flanked by its own LTRs).</li>
</ul>

<h3>What this page does not do</h3>
<p>It does not test causality. Accessibility plus transcription at an ERV is consistent with
an enhancer or an alternative promoter, but nothing here demonstrates that the element drives
a neighbouring gene. The nearest-gene columns in the per-locus tab are proximity only. CRISPRi
or deletion is the experiment that would settle it.</p>
<p>It also inherits the design's confounding: with batch tied to stage and cell type, a
subfamily specific to one group cannot be separated from a batch effect on that group by
statistics alone. Groups with few cells are flagged <code>!</code> in the Overview tab.</p>

<h3>Provenance</h3>
<p>Generated by <code>scripts/23_build_explorer.py</code> from the CSVs under
<code>results/</code>. Upstream scripts, in order: <code>01</code>–<code>09</code> build the
scTE indexes and run IRescue, scTE and the ATAC cut-site counting with pseudobulk DESeq2;
<code>10</code>–<code>13</code> merge and compare the modes and the two ATAC methods;
<code>14</code> assigns tiers and derives int&harr;LTR pairs from genome flanks;
<code>15</code>–<code>21</code> run and interpret locus-mode IRescue; <code>22</code> is the
pseudotime/trajectory analysis. The narrative report with figures is
<code>docs/HERV_progress_report.html</code>.</p>
</div>
"""

DATA = {
    "groups": GROUPS, "stages": STAGES, "fam": fam, "meta": meta,
    "rna": rna, "atac": atac, "cpm": cpm, "igf": igf,
    "tiers": {lt: dict(tiers[lt]) for lt in ("g", "s")},
    "loci": dict(loci), "lsum": lsum, "pairs": dict(pairs),
    "ncand": ncand, "help": HELP,
}

print("rendering ...")
with open(TPL) as fh: tpl = fh.read()
if "__DATA__" not in tpl: sys.exit("template has no __DATA__ placeholder")
blob = json.dumps(DATA, separators=(",", ":"), allow_nan=False)
html = tpl.replace("__DATA__", blob)

# sanity: the JSON must survive a round trip, and the HTML must parse
json.loads(blob)
class V(HTMLParser):
    def __init__(s):
        super().__init__(convert_charrefs=True); s.stack = []; s.bad = []
        s.void = {"br","img","meta","link","input","hr","col","source","area","base","wbr"}
    def handle_starttag(s, t, a):
        if t not in s.void: s.stack.append(t)
    def handle_endtag(s, t):
        if t in s.void: return
        if not s.stack: s.bad.append(f"stray </{t}>"); return
        if s.stack[-1] == t: s.stack.pop()
        elif t in s.stack:
            while s.stack and s.stack[-1] != t: s.bad.append(f"unclosed <{s.stack.pop()}>")
            if s.stack: s.stack.pop()
        else: s.bad.append(f"stray </{t}>")
v = V(); v.feed(html)
if v.bad or v.stack:
    print("  ! HTML problems:", (v.bad + [f"unclosed <{t}>" for t in v.stack])[:12], file=sys.stderr)
else:
    print("  HTML balanced")

with open(OUT, "w") as fh: fh.write(html)
mb = os.path.getsize(OUT) / 1e6
print(f"wrote {OUT}  ({mb:.2f} MB, data {len(blob)/1e6:.2f} MB)")
