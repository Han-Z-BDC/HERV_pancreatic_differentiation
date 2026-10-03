#!/usr/bin/env python3
"""
Build docs/HERV_progress_report.html from docs/report_template.html (standard library only; runs on the Mac).
Placeholders: <!--SCRIPT:path--> (file inlined, HTML-escaped), <!--IMG:path--> (base64 data URI),
<!--VERDICT_ROWS-->, <!--TIER_SUMMARY_ROWS-->, <!--TIER_FULL_ROWS--> (tables from results/ CSVs).
Usage: python3 docs/build_report.py
"""
import base64
import csv
import html
import re
from collections import defaultdict
from pathlib import Path

repo = Path(__file__).resolve().parents[1]
t = (repo / "docs/report_template.html").read_text()
t = re.sub(r"<!--SCRIPT:(.+?)-->", lambda m: html.escape((repo / m.group(1)).read_text().rstrip()), t)


def img(m):
    p = repo / m.group(1)
    mime = "image/png" if p.suffix == ".png" else "image/jpeg"
    return f"data:{mime};base64," + base64.b64encode(p.read_bytes()).decode()


t = re.sub(r"<!--IMG:(.+?)-->", img, t)


def num(x, nd=2):
    try:
        v = float(x)
        return "" if v != v else f"{v:.{nd}f}"
    except ValueError:
        return ""


# read-through verdict table (scTE modes)
cls = {"weak": "warn", "intronic_only": "bad", "no_intergenic": "info"}
rows = []
with open(repo / "results/scte_modes/scte_mode_comparison_group.csv") as f:
    for r in csv.DictReader(f):
        if r["class"] != "LTR" or r["verdict"] in ("", "intergenic"):
            continue
        if not r["baseMean_ex"] or float(r["baseMean_ex"]) <= 20:
            continue
        rows.append((r["verdict"], r["level"], f"<tr><td>{r['subfamily']}</td><td>{r['level']}</td>"
                     + "".join(f'<td class="num">{num(r[k])}</td>' for k in ["log2FoldChange_ir", "log2FoldChange_ex",
                               "log2FoldChange_no", "frac_overall", "frac_level", "frac_ratio"])
                     + f'<td><span class="pill {cls[r["verdict"]]}">{r["verdict"]}</span></td></tr>'))
order = {"intronic_only": 0, "no_intergenic": 1, "weak": 2}
rows.sort(key=lambda x: (order[x[0]], x[1].lower()))
t = t.replace("<!--VERDICT_ROWS-->", "\n".join(r[2] for r in rows))

# candidate tiers
GROUPS = ["hPSC", "DE", "GT", "OTX2_endo", "FG", "liver", "PP1", "PP2", "ENP1", "ENP2", "ENP3", "early_SC_beta",
          "early_EC_beta", "early_SC_alpha", "early_SC_delta", "early_SC_EC", "late_SC_beta", "late_EC_beta",
          "late_SC_alpha", "late_SC_delta", "late_SC_EC", "ductal", "exocrine", "stroma"]
cand = list(csv.DictReader(open(repo / "results/candidates/candidate_tiers_group.csv")))
by = defaultdict(lambda: {"atac": [], "rna": [], "disc": []})
full = []
for r in cand:
    bm = float(r["baseMean_ir"]) if r["baseMean_ir"] else 0
    ex = float(r["log2FoldChange_ex"]) if r["log2FoldChange_ex"] else 0
    lfc = num(r["log2FoldChange_ir"], 1) or num(r["log2FoldChange_ex"], 1)
    atac = r["atac_support"] == "True"
    feat = r["atac_feature"] if r["atac_feature"] != r["subfamily"] else ""
    if bm > 20 or ex > 1.5:
        if r["tier"] == "robust" and atac:
            by[r["level"]]["atac"].append(f'<span class="pill good">{r["subfamily"]} {lfc}{" → " + feat if feat else ""}</span>')
        elif r["tier"] == "robust":
            by[r["level"]]["rna"].append(f'<span class="pill info">{r["subfamily"]} {lfc}</span>')
        else:
            by[r["level"]]["disc"].append(f'<span class="pill bad">{r["subfamily"]} ({r["reasons"]})</span>')
    tier_pill = "good" if r["tier"] == "robust" else "bad"
    full.append((GROUPS.index(r["level"]) if r["level"] in GROUPS else 99, r["tier"] != "robust",
                 f"<tr><td>{r['subfamily']}</td><td>{r['level']}</td><td><span class=\"pill {tier_pill}\">{r['tier']}</span></td>"
                 f"<td>{r['reasons']}</td>"
                 + "".join(f'<td class="num">{num(r[k])}</td>' for k in ["log2FoldChange_ir", "log2FoldChange_ex", "log2FoldChange_no"])
                 + f"<td>{r['atac_feature']}{' ✓' if atac else ''}</td>"
                 + "".join(f'<td class="num">{num(r[k])}</td>' for k in ["atac_lfc_cut", "atac_lfc_scte"]) + "</tr>"))
summ = []
for g in GROUPS:
    if g in by:
        d = by[g]
        summ.append(f"<tr><td>{g}</td><td>{' '.join(d['atac'])}</td><td>{' '.join(d['rna'])}</td><td>{' '.join(d['disc'])}</td></tr>")
t = t.replace("<!--TIER_SUMMARY_ROWS-->", "\n".join(summ))
t = t.replace("<!--TIER_FULL_ROWS-->", "\n".join(x[2] for x in sorted(full, key=lambda x: (x[0], x[1]))))

assert "<!--" not in t, re.findall(r"<!--.{0,40}", t)
out = repo / "docs/HERV_progress_report.html"
out.write_text(t)
print(out, len(t) // 1024, "KB;", len(rows), "verdict rows;", len(full), "candidate rows")
