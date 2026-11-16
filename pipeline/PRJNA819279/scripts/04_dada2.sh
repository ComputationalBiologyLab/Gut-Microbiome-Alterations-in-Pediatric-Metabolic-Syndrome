#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh
mkdir -p 04_dada2

# Parameters from config (MiSeq 2x300 V3V4 after primer removal)
TRUNC_F="${TRUNC_F:-260}"
TRUNC_R="${TRUNC_R:-200}"
MAXEE_F="${MAXEE_F:-2}"
MAXEE_R="${MAXEE_R:-2}"
TRUNCQ="${TRUNCQ:-2}"

qiime_run qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 03_preprocessing/trimmed.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f "$TRUNC_F" \
  --p-trunc-len-r "$TRUNC_R" \
  --p-max-ee-f "$MAXEE_F" \
  --p-max-ee-r "$MAXEE_R" \
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

# Export read-tracking table
qiime_run qiime tools export \
  --input-path 04_dada2/denoising_stats.qza \
  --output-path 04_dada2/denoising_stats_export

cp 04_dada2/denoising_stats_export/stats.tsv results/dada2_read_tracking.tsv

cat > 04_dada2/dada2_parameters.txt <<EOF
method=qiime2 dada2 denoise-paired
trim_left_f=0
trim_left_r=0
trunc_len_f=${TRUNC_F}
trunc_len_r=${TRUNC_R}
max_ee_f=${MAXEE_F}
max_ee_r=${MAXEE_R}
trunc_q=${TRUNCQ}
pooling_method=independent
chimera_method=consensus
n_threads=${NCORES}
feature_type=ASV
EOF

echo "DADA2 complete"
