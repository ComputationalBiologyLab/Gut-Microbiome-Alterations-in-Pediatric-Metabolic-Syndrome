#!/usr/bin/env bash
# Single-end import + cutadapt (pass-through if primers already removed)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh
mkdir -p 03_preprocessing

MANIFEST="00_metadata/fastq_manifest.tsv"
# Paper primers: 515F / 806R (V4). Deposited reads typically already trimmed.
FRONT="GTGCCAGCMGCCGCGGTAA"
# Reverse complement of GGACTACVSGGGTATCTAAT
ADAPTER="ATTAGATACCCSBGTAGTCC"

if [[ ! -f 03_preprocessing/demux.qza ]]; then
  qiime_run qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "$MANIFEST" \
    --input-format SingleEndFastqManifestPhred33V2 \
    --output-path 03_preprocessing/demux.qza
fi

qiime_run qiime demux summarize \
  --i-data 03_preprocessing/demux.qza \
  --o-visualization 03_preprocessing/demux_summary.qzv

qiime_run qiime cutadapt trim-single \
  --i-demultiplexed-sequences 03_preprocessing/demux.qza \
  --p-front "$FRONT" \
  --p-adapter "$ADAPTER" \
  --p-error-rate 0.15 \
  --p-indels \
  --p-match-adapter-wildcards \
  --p-match-read-wildcards \
  --p-minimum-length 50 \
  --p-cores "$NCORES" \
  --o-trimmed-sequences 03_preprocessing/trimmed.qza \
  --verbose

qiime_run qiime demux summarize \
  --i-data 03_preprocessing/trimmed.qza \
  --o-visualization 03_preprocessing/trimmed_summary.qzv

cat > 03_preprocessing/trim_parameters.txt <<EOF
method=qiime2 cutadapt trim-single
forward_primer_515F=${FRONT}
reverse_primer_806R_rc_as_adapter=${ADAPTER}
error_rate=0.15
indels=true
match_adapter_wildcards=true
match_read_wildcards=true
minimum_length=50
discard_untrimmed=false
expected_region=V4
note=Deposited SRA FASTQs already lack front primers; cutadapt is pass-through for most reads.
layout=single-end
EOF

echo "Primer trimming complete (single-end)"
