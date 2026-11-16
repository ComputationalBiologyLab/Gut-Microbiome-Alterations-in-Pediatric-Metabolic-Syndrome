#!/usr/bin/env bash
# Master pipeline for PRJNA1016449 ONLY (no cross-study meta-analysis).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source "$ROOT/config/env.sh"

mkdir -p logs results 02_qc/fastqc 02_qc/multiqc 03_preprocessing 04_dada2 \
  05_taxonomy 06_phylogeny 07_diversity 08_differential_abundance 09_effect_sizes \
  10_functional_prediction 11_figures 12_tables

log() { echo "[$(date '+%F %T')] $*"; }

log "=== 00 Metadata already prepared ==="
[[ -f 00_metadata/sample_metadata.tsv ]] || { echo "Missing metadata"; exit 1; }

log "=== 01 Verify raw data ==="
bash scripts/01_verify_raw_data.sh | tee logs/01_verify_raw_data.log

log "=== 02 FastQC + MultiQC ==="
bash scripts/02_raw_qc.sh | tee logs/02_raw_qc.log

log "=== 03 Import + cutadapt primer removal ==="
bash scripts/03_preprocess_cutadapt.sh | tee logs/03_preprocess_cutadapt.log

log "=== 04 DADA2 ==="
bash scripts/04_dada2.sh | tee logs/04_dada2.log

log "=== 05 Taxonomy (SILVA 138) ==="
bash scripts/05_taxonomy.sh | tee logs/05_taxonomy.log

log "=== 06 Filter contaminants + phylogeny ==="
bash scripts/06_filter_and_phylogeny.sh | tee logs/06_filter_and_phylogeny.log

log "=== 07 Diversity ==="
bash scripts/07_diversity.sh | tee logs/07_diversity.log

log "=== 08–16 Differential abundance, effect sizes, tables, figures, report ==="
bash scripts/08_downstream_r.sh | tee logs/08_downstream_r.log

log "=== DONE: PRJNA1016449 single-study analysis ==="
