#!/usr/bin/env python
"""Build an interactive HTML QC review page from 07_qc_metrics_only.py's per-cell table.

Nothing is filtered anywhere in this pipeline step either -- the page exists so the cutoffs can be
chosen by moving them and watching the cost, then applied once by setting configs/qc_RNA.yaml.
The thresholds it computes are the same rules as the lab pipeline's scripts/1.1_rna_qc.py:
  counts/genes : one-sided MAD, median(log1p(x)) - n_mads * MAD(log1p(x))
  %MT          : median + n_mads * MAD, clamped to [floor, cap]   (or a fixed value, or off)
  %ribo        : plain fixed cutoff, off by default

Usage: python 08_qc_report.py --per-cell qc_per_cell.csv.gz --summary qc_metrics_summary.csv --out report.html
"""
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = {"min_counts": 500, "min_genes": 250, "max_mt_pct": 25.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-cell", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    df = pd.read_csv(a.per_cell)
    summ = pd.read_csv(a.summary)
    kind_of = dict(zip(summ["sample"], summ["kind"]))

    samples = []
    for sid, g in df.groupby("sample", sort=False):
        samples.append({
            "id": str(sid),
            "kind": str(kind_of.get(sid, "")),
            "n": int(len(g)),
            # log1p space is what the MAD rule operates in; keep 3 dp, it is plenty
            "lc": [round(float(x), 3) for x in np.log1p(g["total_counts"].to_numpy())],
            "lg": [round(float(x), 3) for x in np.log1p(g["n_genes_by_counts"].to_numpy())],
            "mt": [round(float(x), 2) for x in g["pct_counts_mt"].to_numpy()],
            "rb": [round(float(x), 2) for x in g["pct_counts_ribo"].to_numpy()],
        })
    print(f"  {len(samples)} samples, {sum(s['n'] for s in samples):,} cells")

    data = {"samples": samples, "paper": PAPER,
            "generated": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")}
    html = TEMPLATE.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    Path(a.out).write_text(html)
    print(f"  wrote {a.out} ({Path(a.out).stat().st_size/1e6:.2f} MB)")


TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Fetal QC Review</title>
<style>
:root{--bg:#fff;--bg2:#f6f7f9;--bg3:#eceef2;--fg:#16181d;--fg2:#5b6472;--fg3:#8b94a3;
 --bd:#dcdfe6;--bd2:#c3c8d2;--accent:#2f5fd0;--accent2:#e8eefc;--keep:#2f6fb5;--drop:#d1495b;
 --warn:#b26a00;--ok:#1d7a43;--mono:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
 --sans:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
@media(prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#13151a;--bg2:#1a1d24;--bg3:#232733;
 --fg:#e6e9ef;--fg2:#a2abba;--fg3:#737d8c;--bd:#2d323d;--bd2:#3d4452;--accent:#6f9bf0;--accent2:#1e2940;
 --keep:#5aa8e0;--drop:#e8705f;--warn:#e0a44a;--ok:#4fc07a}}
:root[data-theme=dark]{--bg:#13151a;--bg2:#1a1d24;--bg3:#232733;--fg:#e6e9ef;--fg2:#a2abba;--fg3:#737d8c;
 --bd:#2d323d;--bd2:#3d4452;--accent:#6f9bf0;--accent2:#1e2940;--keep:#5aa8e0;--drop:#e8705f;--warn:#e0a44a;--ok:#4fc07a}
*{box-sizing:border-box}html,body{margin:0;padding:0}
body{background:var(--bg);color:var(--fg);font:14px/1.55 var(--sans);-webkit-font-smoothing:antialiased}
header{position:sticky;top:0;z-index:50;background:var(--bg);border-bottom:1px solid var(--bd);padding:11px 18px;
 display:flex;align-items:center;gap:14px;flex-wrap:wrap}
header h1{margin:0;font-size:15px;font-weight:650}
.sub{color:var(--fg3);font-size:12px}
.spacer{flex:1}
button,select,input[type=number]{font:inherit}
#themeBtn{border:1px solid var(--bd);background:var(--bg2);color:var(--fg2);border-radius:6px;padding:5px 9px;cursor:pointer;font-size:12px}
main{display:flex;align-items:flex-start}
#ctl{width:310px;flex:0 0 310px;border-right:1px solid var(--bd);background:var(--bg2);padding:14px;
 position:sticky;top:50px;max-height:calc(100vh - 50px);overflow-y:auto}
#body{flex:1;min-width:0;padding:16px 20px 60px}
h2{font-size:14px;font-weight:650;margin:0 0 10px}
h3{font-size:12.5px;font-weight:650;margin:20px 0 7px;color:var(--fg2);text-transform:uppercase;letter-spacing:.05em}
h3:first-child{margin-top:0}
.grp{border:1px solid var(--bd);border-radius:8px;padding:11px;background:var(--bg);margin-bottom:11px}
.grp>.t{font-size:12px;font-weight:650;margin-bottom:8px}
.row{display:flex;align-items:center;gap:8px;margin:7px 0;font-size:12.5px}
.row label{flex:1;color:var(--fg2)}
.row .v{font-family:var(--mono);font-size:12px;font-weight:600;min-width:46px;text-align:right}
input[type=range]{width:100%;accent-color:var(--accent)}
input[type=number]{width:62px;padding:3px 5px;border:1px solid var(--bd2);border-radius:5px;background:var(--bg);color:var(--fg);font-size:12px}
.seg{display:flex;gap:2px;background:var(--bg3);border-radius:6px;padding:2px;margin:6px 0}
.seg button{flex:1;border:0;background:transparent;color:var(--fg2);font-size:11.5px;padding:4px 6px;border-radius:4px;cursor:pointer}
.seg button.on{background:var(--bg);color:var(--fg);font-weight:650}
.hint{font-size:11px;color:var(--fg3);line-height:1.5;margin-top:6px}
.tabs{display:flex;gap:1px;border-bottom:1px solid var(--bd);margin-bottom:14px;flex-wrap:wrap}
.tabs button{border:0;border-bottom:2px solid transparent;background:transparent;color:var(--fg2);font-size:12.5px;
 padding:7px 11px;cursor:pointer;margin-bottom:-1px;font-family:var(--mono)}
.tabs button.on{color:var(--accent);border-bottom-color:var(--accent);font-weight:650}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{padding:5px 9px;text-align:right;border-bottom:1px solid var(--bd);white-space:nowrap}
th{background:var(--bg3);font-size:10.5px;text-transform:uppercase;letter-spacing:.04em;color:var(--fg2);font-weight:650}
td:first-child,th:first-child{text-align:left}
tbody tr:hover{background:var(--bg3)}
.tw{border:1px solid var(--bd);border-radius:8px;overflow:hidden;background:var(--bg2)}
.mono{font-family:var(--mono)}
.drop{color:var(--drop);font-weight:650}
.ok{color:var(--ok)}
.charts{display:grid;grid-template-columns:repeat(auto-fit,minmax(290px,1fr));gap:14px;margin-bottom:16px}
.card{border:1px solid var(--bd);border-radius:8px;padding:10px;background:var(--bg2)}
.card .ct{font-size:11.5px;font-weight:650;color:var(--fg2);margin-bottom:6px;display:flex;justify-content:space-between}
.card .ct em{font-style:normal;color:var(--drop);font-family:var(--mono);font-weight:600}
canvas{display:block;width:100%;height:auto}
pre{background:var(--bg3);border:1px solid var(--bd);border-radius:7px;padding:11px;font-family:var(--mono);
 font-size:11.5px;overflow-x:auto;line-height:1.5;color:var(--fg)}
.legend{display:flex;gap:14px;font-size:11.5px;color:var(--fg2);margin-top:8px;flex-wrap:wrap}
.legend i{display:inline-block;width:11px;height:11px;border-radius:2px;vertical-align:-1px;margin-right:4px}
.note{font-size:12px;color:var(--fg3);max-width:76ch;margin:0 0 12px}
@media(max-width:860px){main{flex-direction:column}#ctl{width:100%;flex:none;position:static;max-height:none;border-right:0;border-bottom:1px solid var(--bd)}}
</style></head><body>
<header><h1>Fetal pancreas — QC cutoff review</h1>
<span class="sub" id="hsub"></span><div class="spacer"></div>
<button id="themeBtn">Theme</button></header>
<main>
<aside id="ctl">
  <h2>Thresholds</h2>
  <div class="grp"><div class="t">Counts &amp; genes</div>
    <div class="row"><label>n_mads</label><span class="v" id="vNm">3.0</span></div>
    <input type="range" id="nm" min="1" max="6" step="0.5" value="3">
    <div class="hint">One-sided MAD on log1p, as in <code>is_mad_outlier()</code>. Flags cells
      <i>below</i> median − n·MAD only; high counts are a doublet question, not this step's.</div>
  </div>
  <div class="grp"><div class="t">% MT</div>
    <div class="seg" id="mtMode">
      <button data-m="mad" class="on">MAD</button><button data-m="fixed">Fixed</button><button data-m="off">Off</button>
    </div>
    <div id="mtMad">
      <div class="row"><label>n_mads</label><span class="v" id="vMtNm">5.0</span></div>
      <input type="range" id="mtNm" min="1" max="8" step="0.5" value="5">
      <div class="row"><label>floor %</label><input type="number" id="mtFloor" value="3" step="0.5" min="0"></div>
      <div class="row"><label>cap %</label><input type="number" id="mtCap" value="15" step="0.5" min="0"></div>
    </div>
    <div id="mtFixedBox" style="display:none">
      <div class="row"><label>max %MT</label><input type="number" id="mtFixed" value="10" step="0.5" min="0"></div>
    </div>
    <div class="hint">median + n·MAD, clamped to [floor, cap]. Upper bound only.</div>
  </div>
  <div class="grp"><div class="t">% ribo</div>
    <div class="seg" id="rbMode"><button data-m="off" class="on">Off</button><button data-m="fixed">Fixed</button></div>
    <div id="rbBox" style="display:none">
      <div class="row"><label>max %ribo</label><input type="number" id="rbMax" value="50" step="1" min="0"></div>
    </div>
    <div class="hint">Off by default — high ribosomal content is often real biology in progenitors.</div>
  </div>
  <div class="grp"><div class="t">Reference</div>
    <div class="row"><label><input type="checkbox" id="showPaper" checked> show paper cutoff</label></div>
    <div class="hint">Olaniru 2023: &gt;500 UMI, &gt;250 genes, &lt;25% MT.</div>
  </div>
</aside>
<section id="body">
  <div class="tabs" id="tabs"></div>
  <div id="view"></div>
</section>
</main>
<script>
const D = __DATA__;
const $ = s => document.querySelector(s);
const S = {nm:3, mtMode:'mad', mtNm:5, mtFloor:3, mtCap:15, mtFixed:10, rbMode:'off', rbMax:50, paper:true, tab:'ALL'};

/* ---------- stats ---------- */
function median(a){const b=Float64Array.from(a);b.sort();const n=b.length;return n%2?b[(n-1)/2]:(b[n/2-1]+b[n/2])/2}
function madOf(a,m){const d=new Float64Array(a.length);for(let i=0;i<a.length;i++)d[i]=Math.abs(a[i]-m);return median(d)}
const cache = new Map();
function stats(s){
  if(cache.has(s.id)) return cache.get(s.id);
  const o={};
  for(const k of ['lc','lg','mt','rb']){const m=median(s[k]);o[k]={med:m,mad:madOf(s[k],m)}}
  cache.set(s.id,o); return o;
}
function cuts(s){
  const st=stats(s);
  const lcCut = st.lc.mad===0?null:st.lc.med-S.nm*st.lc.mad;
  const lgCut = st.lg.mad===0?null:st.lg.med-S.nm*st.lg.mad;
  let mtCut=null;
  if(S.mtMode==='mad'){ mtCut=st.mt.med+S.mtNm*st.mt.mad;
    if(S.mtFloor!==''&&!isNaN(S.mtFloor)) mtCut=Math.max(mtCut,+S.mtFloor);
    if(S.mtCap!==''&&!isNaN(S.mtCap)) mtCut=Math.min(mtCut,+S.mtCap); }
  else if(S.mtMode==='fixed') mtCut=+S.mtFixed;
  const rbCut = S.rbMode==='fixed'? +S.rbMax : null;
  return {lcCut,lgCut,mtCut,rbCut};
}
function evaluate(s){
  const c=cuts(s), n=s.n;
  const fail=new Uint8Array(n); let fLc=0,fLg=0,fMt=0,fRb=0;
  for(let i=0;i<n;i++){
    let f=0;
    if(c.lcCut!==null&&s.lc[i]<c.lcCut){f=1;fLc++}
    if(c.lgCut!==null&&s.lg[i]<c.lgCut){f=1;fLg++}
    if(c.mtCut!==null&&s.mt[i]>c.mtCut){f=1;fMt++}
    if(c.rbCut!==null&&s.rb[i]>c.rbCut){f=1;fRb++}
    fail[i]=f;
  }
  let lost=0; for(let i=0;i<n;i++) lost+=fail[i];
  // paper rule, for comparison
  let pLost=0;
  const pc=Math.log1p(D.paper.min_counts), pg=Math.log1p(D.paper.min_genes);
  for(let i=0;i<n;i++) if(s.lc[i]<=pc||s.lg[i]<=pg||s.mt[i]>=D.paper.max_mt_pct) pLost++;
  return {c,fail,lost,fLc,fLg,fMt,fRb,pLost};
}

/* ---------- drawing ---------- */
const CSS = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
function hist(canvas, vals, fail, cut, opts={}){
  const dpr=Math.min(2,devicePixelRatio||1), W=canvas.clientWidth, H=opts.h||150;
  canvas.width=W*dpr; canvas.height=H*dpr;
  const x=canvas.getContext('2d'); x.scale(dpr,dpr); x.clearRect(0,0,W,H);
  const ml=34,mr=8,mt=8,mb=22, iw=W-ml-mr, ih=H-mt-mb;
  let lo=Infinity,hi=-Infinity; for(const v of vals){if(v<lo)lo=v;if(v>hi)hi=v}
  if(!(hi>lo)) hi=lo+1;
  const NB=opts.bins||60, bw=(hi-lo)/NB;
  const keep=new Float64Array(NB), drop=new Float64Array(NB);
  for(let i=0;i<vals.length;i++){
    let b=Math.floor((vals[i]-lo)/bw); if(b<0)b=0; if(b>=NB)b=NB-1;
    (fail[i]?drop:keep)[b]++;
  }
  let mx=0; for(let b=0;b<NB;b++) mx=Math.max(mx,keep[b]+drop[b]);
  if(mx===0)mx=1;
  const px=v=>ml+iw*(v-lo)/(hi-lo), py=v=>mt+ih*(1-v/mx);
  const cK=CSS('--keep')||'#2f6fb5', cD=CSS('--drop')||'#d1495b';
  for(let b=0;b<NB;b++){
    const x0=ml+iw*b/NB, w=Math.max(1,iw/NB-0.6);
    if(keep[b]){x.fillStyle=cK; x.fillRect(x0,py(keep[b]),w,ih-(py(keep[b])-mt))}
    if(drop[b]){x.fillStyle=cD; const h=ih*drop[b]/mx; x.fillRect(x0,py(keep[b]+drop[b]),w,h)}
  }
  x.strokeStyle=CSS('--bd2')||'#c3c8d2'; x.lineWidth=1;
  x.beginPath(); x.moveTo(ml,mt+ih); x.lineTo(W-mr,mt+ih); x.stroke();
  if(cut!==null&&cut>=lo&&cut<=hi){
    x.strokeStyle=cD; x.lineWidth=1.8; x.setLineDash([4,3]);
    x.beginPath(); x.moveTo(px(cut),mt); x.lineTo(px(cut),mt+ih); x.stroke(); x.setLineDash([]);
  }
  if(opts.paper!==undefined&&S.paper&&opts.paper>=lo&&opts.paper<=hi){
    x.strokeStyle=CSS('--fg3')||'#888'; x.lineWidth=1.3; x.setLineDash([2,3]);
    x.beginPath(); x.moveTo(px(opts.paper),mt); x.lineTo(px(opts.paper),mt+ih); x.stroke(); x.setLineDash([]);
  }
  x.fillStyle=CSS('--fg3')||'#888'; x.font='9px '+CSS('--sans'); x.textAlign='center';
  for(let i=0;i<=4;i++){const v=lo+(hi-lo)*i/4; x.fillText(opts.fmt?opts.fmt(v):v.toFixed(1), px(v), H-7)}
  x.textAlign='right'; x.fillText(String(Math.round(mx)), ml-4, mt+8);
}
function scatter(canvas, s, fail){
  const dpr=Math.min(2,devicePixelRatio||1), W=canvas.clientWidth, H=300;
  canvas.width=W*dpr; canvas.height=H*dpr;
  const x=canvas.getContext('2d'); x.scale(dpr,dpr); x.clearRect(0,0,W,H);
  const ml=38,mr=10,mt=10,mb=28, iw=W-ml-mr, ih=H-mt-mb;
  let xlo=Infinity,xhi=-Infinity,ylo=Infinity,yhi=-Infinity;
  for(let i=0;i<s.n;i++){if(s.lc[i]<xlo)xlo=s.lc[i];if(s.lc[i]>xhi)xhi=s.lc[i];
    if(s.lg[i]<ylo)ylo=s.lg[i];if(s.lg[i]>yhi)yhi=s.lg[i]}
  const px=v=>ml+iw*(v-xlo)/(xhi-xlo||1), py=v=>mt+ih*(1-(v-ylo)/(yhi-ylo||1));
  const cK=CSS('--keep')||'#2f6fb5', cD=CSS('--drop')||'#d1495b';
  // draw kept first so dropped stay visible on top
  for(const pass of [0,1]){
    x.fillStyle = pass?cD:cK; x.globalAlpha = pass?0.75:0.33;
    for(let i=0;i<s.n;i++){ if(fail[i]!==pass) continue;
      x.fillRect(px(s.lc[i])-0.9, py(s.lg[i])-0.9, 1.8, 1.8); }
  }
  x.globalAlpha=1;
  const c=cuts(s);
  x.strokeStyle=cD; x.lineWidth=1.4; x.setLineDash([4,3]);
  if(c.lcCut!==null){x.beginPath();x.moveTo(px(c.lcCut),mt);x.lineTo(px(c.lcCut),mt+ih);x.stroke()}
  if(c.lgCut!==null){x.beginPath();x.moveTo(ml,py(c.lgCut));x.lineTo(W-mr,py(c.lgCut));x.stroke()}
  x.setLineDash([]);
  x.strokeStyle=CSS('--bd2'); x.lineWidth=1;
  x.beginPath();x.moveTo(ml,mt+ih);x.lineTo(W-mr,mt+ih);x.moveTo(ml,mt);x.lineTo(ml,mt+ih);x.stroke();
  x.fillStyle=CSS('--fg3'); x.font='9.5px '+CSS('--sans'); x.textAlign='center';
  x.fillText('log1p(total counts)', ml+iw/2, H-6);
  x.save(); x.translate(11, mt+ih/2); x.rotate(-Math.PI/2); x.fillText('log1p(n genes)',0,0); x.restore();
}

/* ---------- render ---------- */
function summaryTable(){
  let tot=0, lost=0, plost=0;
  const rows = D.samples.map(s=>{
    const e=evaluate(s); tot+=s.n; lost+=e.lost; plost+=e.pLost;
    const cc=e.c.lcCut===null?'–':Math.round(Math.expm1(e.c.lcCut)), gc=e.c.lgCut===null?'–':Math.round(Math.expm1(e.c.lgCut));
    return `<tr><td class="mono">${s.id}</td><td style="text-align:left;font-size:11px;color:var(--fg3)">${s.kind}</td>
      <td class="mono">${s.n.toLocaleString()}</td>
      <td class="mono">${cc}</td><td class="mono">${gc}</td>
      <td class="mono">${e.c.mtCut===null?'off':e.c.mtCut.toFixed(1)+'%'}</td>
      <td class="mono">${e.fLc.toLocaleString()}</td><td class="mono">${e.fLg.toLocaleString()}</td>
      <td class="mono">${e.fMt.toLocaleString()}</td><td class="mono">${e.c.rbCut===null?'–':e.fRb.toLocaleString()}</td>
      <td class="mono drop">${e.lost.toLocaleString()}</td>
      <td class="mono drop">${(100*e.lost/s.n).toFixed(1)}%</td>
      <td class="mono" style="color:var(--fg3)">${e.pLost.toLocaleString()} (${(100*e.pLost/s.n).toFixed(1)}%)</td>
      <td class="mono ok">${(s.n-e.lost).toLocaleString()}</td></tr>`;
  }).join('');
  return `<div class="tw"><table><thead><tr>
    <th>sample</th><th style="text-align:left">kind</th><th>cells</th>
    <th>UMI cut</th><th>gene cut</th><th>%MT cut</th>
    <th>fail UMI</th><th>fail gene</th><th>fail MT</th><th>fail ribo</th>
    <th>lost</th><th>lost %</th><th>paper rule</th><th>kept</th></tr></thead>
    <tbody>${rows}</tbody>
    <tfoot><tr style="font-weight:650;background:var(--bg3)">
      <td>TOTAL</td><td></td><td class="mono">${tot.toLocaleString()}</td>
      <td colspan="6"></td>
      <td class="mono drop">${lost.toLocaleString()}</td><td class="mono drop">${(100*lost/tot).toFixed(1)}%</td>
      <td class="mono" style="color:var(--fg3)">${plost.toLocaleString()} (${(100*plost/tot).toFixed(1)}%)</td>
      <td class="mono ok">${(tot-lost).toLocaleString()}</td></tr></tfoot></table></div>`;
}
function yamlBlock(){
  const mt = S.mtMode==='off'
    ? `mt_pct_strategy: fixed\nmt_pct_max: null            # MT filtering off`
    : S.mtMode==='fixed'
      ? `mt_pct_strategy: fixed\nmt_pct_max: ${(+S.mtFixed).toFixed(1)}`
      : `mt_pct_strategy: mad\nmt_pct_n_mads: ${(+S.mtNm)}\nmt_pct_floor: ${(+S.mtFloor).toFixed(1)}\nmt_pct_cap: ${(+S.mtCap).toFixed(1)}`;
  const rb = S.rbMode==='off' ? 'ribo_pct_max: null' : `ribo_pct_max: ${(+S.rbMax).toFixed(1)}`;
  return `# configs/qc_RNA.yaml -- settings currently selected on this page\nfiltering_strategy: mad\nn_mads: ${S.nm}\n\n${mt}\n\n${rb}\n`;
}
function render(){
  const v=$('#view');
  if(S.tab==='ALL'){
    v.innerHTML = `<p class="note">Move the controls on the left; every number and chart below updates.
      Nothing is filtered — this page only shows what each choice would cost. The
      <b>paper rule</b> column is Olaniru's fixed cutoffs, for reference.</p>
      <h3>All samples at the current settings</h3>${summaryTable()}
      <h3>Apply this choice</h3>
      <p class="note">Paste into <code>configs/qc_RNA.yaml</code> in the single-cell-pipeline repo.</p>
      <pre id="yml">${yamlBlock()}</pre>`;
    return;
  }
  const s=D.samples.find(x=>x.id===S.tab); if(!s) return;
  const e=evaluate(s), c=e.c;
  v.innerHTML = `<p class="note"><b class="mono">${s.id}</b> — ${s.kind}, ${s.n.toLocaleString()} cells called by Cell Ranger.
      <span class="drop">${e.lost.toLocaleString()} (${(100*e.lost/s.n).toFixed(1)}%)</span> would be removed at the current settings.</p>
    <div class="charts">
      <div class="card"><div class="ct"><span>log1p(total counts)</span><em>${e.fLc.toLocaleString()} fail</em></div><canvas id="c1"></canvas></div>
      <div class="card"><div class="ct"><span>log1p(n genes)</span><em>${e.fLg.toLocaleString()} fail</em></div><canvas id="c2"></canvas></div>
      <div class="card"><div class="ct"><span>% MT counts</span><em>${c.mtCut===null?'off':e.fMt.toLocaleString()+' fail'}</em></div><canvas id="c3"></canvas></div>
      <div class="card"><div class="ct"><span>% ribo counts</span><em>${c.rbCut===null?'not filtered':e.fRb.toLocaleString()+' fail'}</em></div><canvas id="c4"></canvas></div>
    </div>
    <div class="card"><div class="ct"><span>counts vs genes</span><span style="color:var(--fg3);font-weight:400">dashed = current cutoffs</span></div>
      <canvas id="sc"></canvas></div>
    <div class="legend"><span><i style="background:var(--keep)"></i>kept</span>
      <span><i style="background:var(--drop)"></i>removed</span>
      ${S.paper?'<span style="color:var(--fg3)">dotted grey = paper cutoff</span>':''}</div>
    <h3>Thresholds in linear units</h3>
    <div class="tw"><table><tbody>
      <tr><td>min total counts</td><td class="mono">${c.lcCut===null?'–':Math.round(Math.expm1(c.lcCut)).toLocaleString()}</td></tr>
      <tr><td>min genes</td><td class="mono">${c.lgCut===null?'–':Math.round(Math.expm1(c.lgCut)).toLocaleString()}</td></tr>
      <tr><td>max %MT</td><td class="mono">${c.mtCut===null?'off':c.mtCut.toFixed(2)+'%'}</td></tr>
      <tr><td>max %ribo</td><td class="mono">${c.rbCut===null?'off (diagnostic only)':c.rbCut.toFixed(1)+'%'}</td></tr>
    </tbody></table></div>`;
  hist($('#c1'), s.lc, e.fail, c.lcCut, {paper:Math.log1p(D.paper.min_counts)});
  hist($('#c2'), s.lg, e.fail, c.lgCut, {paper:Math.log1p(D.paper.min_genes)});
  hist($('#c3'), s.mt, e.fail, c.mtCut, {paper:D.paper.max_mt_pct, fmt:v=>v.toFixed(0)+'%'});
  hist($('#c4'), s.rb, e.fail, c.rbCut, {fmt:v=>v.toFixed(0)+'%'});
  scatter($('#sc'), s, e.fail);
}

/* ---------- wiring ---------- */
function tabs(){
  const t=$('#tabs'); t.innerHTML='';
  const mk=(id,lab)=>{const b=document.createElement('button');b.textContent=lab;b.className=id===S.tab?'on':'';
    b.onclick=()=>{S.tab=id;tabs();render()};t.append(b)};
  mk('ALL','All samples');
  D.samples.forEach(s=>mk(s.id,s.id));
}
function bind(){
  $('#nm').oninput=e=>{S.nm=+e.target.value;$('#vNm').textContent=S.nm.toFixed(1);render()};
  $('#mtNm').oninput=e=>{S.mtNm=+e.target.value;$('#vMtNm').textContent=S.mtNm.toFixed(1);render()};
  for(const [id,k] of [['mtFloor','mtFloor'],['mtCap','mtCap'],['mtFixed','mtFixed'],['rbMax','rbMax']])
    $('#'+id).oninput=e=>{S[k]=e.target.value;render()};
  $('#mtMode').onclick=e=>{const b=e.target.closest('button');if(!b)return;
    [...$('#mtMode').children].forEach(c=>c.classList.toggle('on',c===b));
    S.mtMode=b.dataset.m;
    $('#mtMad').style.display=S.mtMode==='mad'?'':'none';
    $('#mtFixedBox').style.display=S.mtMode==='fixed'?'':'none';
    render()};
  $('#rbMode').onclick=e=>{const b=e.target.closest('button');if(!b)return;
    [...$('#rbMode').children].forEach(c=>c.classList.toggle('on',c===b));
    S.rbMode=b.dataset.m; $('#rbBox').style.display=S.rbMode==='fixed'?'':'none'; render()};
  $('#showPaper').onchange=e=>{S.paper=e.target.checked;render()};
  $('#themeBtn').onclick=()=>{const cur=document.documentElement.getAttribute('data-theme');
    const dark=cur?cur==='dark':matchMedia('(prefers-color-scheme:dark)').matches;
    document.documentElement.setAttribute('data-theme',dark?'light':'dark'); render()};
  addEventListener('resize',()=>{clearTimeout(window._rt);window._rt=setTimeout(render,120)});
}
$('#hsub').textContent = `${D.samples.length} samples · ${D.samples.reduce((a,s)=>a+s.n,0).toLocaleString()} cells · nothing filtered · ${D.generated}`;
tabs(); bind(); render();
</script></body></html>
"""

if __name__ == "__main__":
    main()
