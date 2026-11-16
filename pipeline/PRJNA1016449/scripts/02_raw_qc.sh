#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh
mkdir -p 02_qc/fastqc 02_qc/multiqc results

# FastQC on all paired reads
fastqc -t "$NCORES" -o 02_qc/fastqc 01_raw_data/*.fastq

# MultiQC skipped (user request: FastQC is sufficient)
echo "MultiQC skipped by configuration/user request"

# Summarize read counts / lengths from available FastQC zips (no sample removal)
python3 scripts/summarize_raw_qc.py \
  --metadata 00_metadata/sample_metadata.tsv \
  --fastqc-dir 02_qc/fastqc \
  --out results/raw_read_qc_summary.tsv || true

echo "Raw QC complete (FastQC only; no MultiQC)"
