#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh
mkdir -p 02_qc/fastqc results

# Absolute FASTQ directory (sibling folder; do not rely on symlink resolution)
FQDIR="$(cd "$ROOT/.." && pwd)/PRJNA838778_fastq"
echo "FQDIR=$FQDIR"
ls "$FQDIR"/SRR*_1.fastq.gz >/dev/null
N=$(ls "$FQDIR"/SRR*_1.fastq.gz | wc -l | tr -d ' ')
echo "FastQC inputs: $N"
# Use array to handle spaces in path
FILES=()
for f in "$FQDIR"/SRR*_1.fastq.gz; do FILES+=("$f"); done
fastqc -t "$NCORES" -o 02_qc/fastqc "${FILES[@]}"

echo "MultiQC skipped by configuration/user request"
python3 scripts/summarize_raw_qc.py \
  --metadata 00_metadata/sample_metadata.tsv \
  --fastqc-dir 02_qc/fastqc \
  --out results/raw_read_qc_summary.tsv || true
echo "Raw QC complete (FastQC only; no MultiQC)"
