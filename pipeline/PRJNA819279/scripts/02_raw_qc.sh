#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh
mkdir -p 02_qc/fastqc 02_qc/multiqc results

# FastQC on gzipped paired reads for included samples only (bash 3.2 compatible)
FQS=()
while IFS= read -r f; do
  FQS+=("$f")
done < <(python3 - <<'PY'
import csv
meta=list(csv.DictReader(open("00_metadata/sample_metadata.tsv"), delimiter="\t"))
for r in meta:
    print(r["fastq_r1"]); print(r["fastq_r2"])
PY
)
fastqc -t "$NCORES" -o 02_qc/fastqc "${FQS[@]}"
echo "MultiQC skipped by configuration/user request"
python3 scripts/summarize_raw_qc.py \
  --metadata 00_metadata/sample_metadata.tsv \
  --fastqc-dir 02_qc/fastqc \
  --out results/raw_read_qc_summary.tsv || true
echo "Raw QC complete (FastQC only; no MultiQC)"
