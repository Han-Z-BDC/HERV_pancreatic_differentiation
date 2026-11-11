#!/bin/bash
# Download raw reads of Olaniru et al. 2023 Cell Metab (human fetal pancreas, 10x 3' v3 scRNA + HTO,
# 8-20 PCW; GEO GSE197064, BioProject PRJNA808751) from ENA, for TE / ERV replication.
#
# Default = DRY RUN: writes the manifest (runs, sample titles, files, sizes, md5) and checks free space.
# Nothing is downloaded until you add --go.
#   bash scripts/fetal/01_download_olaniru2023.sh          # manifest + space check only
#   sbatch scripts/fetal/01_download_olaniru2023.sb --go   # download (acompile has internet)
#
# Prefers ENA "submitted" files when they are 10x BAMs (keep CB/UB tags -> straight into IRescue);
# otherwise the FASTQs. Resumable (wget -c) and md5-checked.
set -euo pipefail
GO=${1:-}
ACC=PRJNA808751
OUT=${FETAL_OUT:-/scratch/alpine/hzhu10@xsede.org/fetal_pancreas}/olaniru2023
mkdir -p "$OUT/raw"
MAN="$OUT/manifest.tsv"

curl -sS --retry 5 --retry-delay 30 \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&format=tsv&fields=run_accession,sample_accession,sample_title,experiment_title,library_strategy,read_count,fastq_ftp,fastq_bytes,fastq_md5,submitted_ftp,submitted_bytes,submitted_md5,submitted_format" \
  > "$MAN"
echo "manifest: $MAN ($(($(wc -l < "$MAN") - 1)) runs)"
column -t -s $'\t' <(cut -f1,3,5,6,13 "$MAN") | head -40

# total sizes (GB) for each option
python3 - "$MAN" <<'EOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
def tot(col):
    return sum(sum(int(x) for x in r[col].split(";") if x) for r in rows) / 1e9
print(f"FASTQ total:     {tot('fastq_bytes'):8.1f} GB")
print(f"submitted total: {tot('submitted_bytes'):8.1f} GB  (formats: {sorted({r['submitted_format'] for r in rows})})")
EOF
echo "free space at $OUT:"; df -h "$OUT" | tail -1
command -v curc-quota >/dev/null && curc-quota || true

if [[ "$GO" != "--go" ]]; then
  echo "dry run only; rerun with --go to download"; exit 0
fi

# choose per run: submitted BAM if available, else FASTQ
python3 - "$MAN" > "$OUT/download_list.tsv" <<'EOF'
import csv, sys
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    sub = [x for x in r["submitted_ftp"].split(";") if x]
    if sub and all(x.endswith(".bam") for x in sub):
        urls, md5s = sub, r["submitted_md5"].split(";")
    else:
        urls, md5s = [x for x in r["fastq_ftp"].split(";") if x], r["fastq_md5"].split(";")
    for u, m in zip(urls, md5s):
        print(r["run_accession"], "https://" + u, m, sep="\t")
EOF
echo "$(wc -l < "$OUT/download_list.tsv") files to fetch"

# ENA throttles each connection (~1 MB/s), not the total: fetch up to $NPAR files at once
NPAR=${NPAR:-8}
fetch() {
  run=$1; url=$2; md5=$3
  d="$OUT/raw/$run"; mkdir -p "$d"; f="$d/$(basename "$url")"
  if [[ -s "$f" ]] && [[ "$(md5sum "$f" | cut -d' ' -f1)" == "$md5" ]]; then echo "ok   $f"; return; fi
  wget -c -q --tries=10 --waitretry=60 -O "$f" "$url"
  [[ "$(md5sum "$f" | cut -d' ' -f1)" == "$md5" ]] && echo "done $f" || echo "MD5 MISMATCH $f"
}
export -f fetch; export OUT
xargs -P "$NPAR" -L 1 bash -c 'fetch "$@"' _ < "$OUT/download_list.tsv"
du -sh "$OUT/raw"
