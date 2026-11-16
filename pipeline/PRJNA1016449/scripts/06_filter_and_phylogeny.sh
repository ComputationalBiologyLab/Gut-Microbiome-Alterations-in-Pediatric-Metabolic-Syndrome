#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh
mkdir -p 06_phylogeny results

# Keep raw ASV table export
qiime_run qiime tools export \
  --input-path 04_dada2/table.qza \
  --output-path 04_dada2/table_export
qiime_run biom convert \
  -i 04_dada2/table_export/feature-table.biom \
  -o results/ASV_table_raw.tsv \
  --to-tsv

# Filter mitochondria / chloroplast
qiime_run qiime taxa filter-table \
  --i-table 04_dada2/table.qza \
  --i-taxonomy 05_taxonomy/taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table 06_phylogeny/table_no_mito_chloro.qza

qiime_run qiime taxa filter-seqs \
  --i-sequences 04_dada2/rep_seqs.qza \
  --i-taxonomy 05_taxonomy/taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences 06_phylogeny/rep_seqs_no_mito_chloro.qza

# Prevalence / abundance filter for analysis table
# ASV must appear in >=5% of samples (ceil) and have >=10 total reads
# Computed dynamically in Python helper for exact sample count
N_SAMPLES=$(tail -n +2 00_metadata/qiime_metadata.tsv | wc -l | tr -d ' ')
MIN_SAMPLES=$(python3 - <<PY
import math
n=int("$N_SAMPLES")
print(max(2, math.ceil(0.05*n)))
PY
)

qiime_run qiime feature-table filter-features \
  --i-table 06_phylogeny/table_no_mito_chloro.qza \
  --p-min-samples "$MIN_SAMPLES" \
  --p-min-frequency 10 \
  --o-filtered-table 06_phylogeny/table_filtered.qza

qiime_run qiime feature-table filter-seqs \
  --i-data 06_phylogeny/rep_seqs_no_mito_chloro.qza \
  --i-table 06_phylogeny/table_filtered.qza \
  --o-filtered-data 06_phylogeny/rep_seqs_filtered.qza

qiime_run qiime tools export \
  --input-path 06_phylogeny/table_filtered.qza \
  --output-path 06_phylogeny/table_filtered_export
qiime_run biom convert \
  -i 06_phylogeny/table_filtered_export/feature-table.biom \
  -o results/ASV_table_filtered.tsv \
  --to-tsv

# Phylogenetic tree (MAFFT + FastTree) for UniFrac
qiime_run qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 06_phylogeny/rep_seqs_filtered.qza \
  --p-n-threads "$NCORES" \
  --o-alignment 06_phylogeny/aligned_rep_seqs.qza \
  --o-masked-alignment 06_phylogeny/masked_aligned_rep_seqs.qza \
  --o-tree 06_phylogeny/unrooted_tree.qza \
  --o-rooted-tree 06_phylogeny/rooted_tree.qza

qiime_run qiime tools export \
  --input-path 06_phylogeny/rooted_tree.qza \
  --output-path 06_phylogeny/tree_export
cp 06_phylogeny/tree_export/tree.nwk results/rooted_tree.nwk

cat > 06_phylogeny/filtering_parameters.txt <<EOF
remove_mitochondria=true
remove_chloroplast=true
min_samples_asv=${MIN_SAMPLES}  # ceil(0.05 * n_samples)
min_frequency_asv=10
raw_table=results/ASV_table_raw.tsv
filtered_table=results/ASV_table_filtered.tsv
phylogeny=MAFFT alignment + FastTree via qiime phylogeny align-to-tree-mafft-fasttree
EOF

echo "Filtering + phylogeny complete"
