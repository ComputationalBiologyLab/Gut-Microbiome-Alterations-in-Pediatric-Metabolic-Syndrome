#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh
mkdir -p 04_dada2 results

# Variable-length SE reads (~212–430 bp). trunc_len=0 disables length truncation
# (trunc_len=400 previously dropped all reads from multiple samples).
TRUNC_LEN="${TRUNC_LEN:-0}"
TRIM_LEFT="${TRIM_LEFT:-0}"
MAXEE="${MAXEE:-2}"
TRUNCQ="${TRUNCQ:-2}"

qiime_run qiime dada2 denoise-single \
  --i-demultiplexed-seqs 03_preprocessing/trimmed.qza \
  --p-trim-left "$TRIM_LEFT" \
  --p-trunc-len "$TRUNC_LEN" \
  --p-max-ee "$MAXEE" \
  --p-trunc-q "$TRUNCQ" \
  --p-pooling-method independent \
  --p-chimera-method consensus \
  --p-n-threads "$NCORES" \
  --o-table 04_dada2/table.qza \
  --o-representative-sequences 04_dada2/rep_seqs.qza \
  --o-denoising-stats 04_dada2/denoising_stats.qza \
  --verbose

qiime_run qiime metadata tabulate \
  --m-input-file 04_dada2/denoising_stats.qza \
  --o-visualization 04_dada2/denoising_stats.qzv

qiime_run qiime feature-table summarize \
  --i-table 04_dada2/table.qza \
  --m-sample-metadata-file 00_metadata/qiime_metadata.tsv \
  --o-visualization 04_dada2/table_summary.qzv

qiime_run qiime feature-table tabulate-seqs \
  --i-data 04_dada2/rep_seqs.qza \
  --o-visualization 04_dada2/rep_seqs.qzv

qiime_run qiime tools export \
  --input-path 04_dada2/denoising_stats.qza \
  --output-path 04_dada2/denoising_stats_export

cp 04_dada2/denoising_stats_export/stats.tsv results/dada2_read_tracking.tsv

cat > 04_dada2/dada2_parameters.txt <<EOF
method=qiime2 dada2 denoise-single
trim_left=${TRIM_LEFT}
trunc_len=${TRUNC_LEN}
max_ee=${MAXEE}
trunc_q=${TRUNCQ}
pooling_method=independent
chimera_method=consensus
n_threads=${NCORES}
feature_type=ASV
layout=single-end
note=trunc_len=0 required because deposited SE reads are variable length (often <400bp)
EOF

echo "DADA2 complete (single-end)"
