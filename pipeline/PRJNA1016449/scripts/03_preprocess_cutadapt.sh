#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh
mkdir -p 03_preprocessing

MANIFEST="00_metadata/fastq_manifest.tsv"
FWD="ACTCCTACGGGAGGCAGCAG"
REV="GGACTACHVGGGTWTCTAAT"

if [[ ! -f 03_preprocessing/demux.qza ]]; then
  qiime_run qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST" \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path 03_preprocessing/demux.qza
fi

qiime_run qiime demux summarize \
  --i-data 03_preprocessing/demux.qza \
  --o-visualization 03_preprocessing/demux_summary.qzv

# Remove primers (+ any 5' bases upstream of primer match = barcodes/linkers)
qiime_run qiime cutadapt trim-paired \
  --i-demultiplexed-sequences 03_preprocessing/demux.qza \
  --p-front-f "$FWD" \
  --p-front-r "$REV" \
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

# Persist exact trim parameters
cat > 03_preprocessing/trim_parameters.txt <<EOF
method=qiime2 cutadapt trim-paired
forward_primer_338F=${FWD}
reverse_primer_806R=${REV}
error_rate=0.15
indels=true
match_adapter_wildcards=true
match_read_wildcards=true
minimum_length=50
discard_untrimmed=false
expected_region=V3V4
barcode_note=Demultiplexed SRA reads contain ~8-15nt 5' barcode/linker before primer; cutadapt removes upstream+primer when front adapter matches.
EOF

echo "Primer trimming complete"
