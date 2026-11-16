#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh
mkdir -p 07_diversity

# Choose rarefaction depth: 20th percentile of sample depths (or min if higher),
# but never below 5000 unless data force it. Do not drop samples silently —
# samples below depth are excluded from rarefied diversity and listed.
DEPTH=$(qiime_run python - <<'PY'
import biom, numpy as np
t = biom.load_table('06_phylogeny/table_filtered_export/feature-table.biom')
depths = t.sum(axis='sample')
p20 = int(np.percentile(depths, 20))
d = max(5000, p20) if p20 >= 5000 else max(1000, int(np.percentile(depths, 10)))
# If too many samples would be lost, relax to 10th percentile
n_keep = sum(1 for x in depths if x >= d)
if n_keep < 0.8 * len(depths):
    d = max(1000, int(np.percentile(depths, 10)))
print(d)
open('07_diversity/rarefaction_depth.txt','w').write(str(d)+'\n')
print('depths min/median/max', int(depths.min()), int(np.median(depths)), int(depths.max()), file=__import__('sys').stderr)
print('chosen depth', d, 'n_keep', sum(1 for x in depths if x >= d), file=__import__('sys').stderr)
PY
)

echo "Rarefaction depth: $DEPTH"

qiime_run qiime diversity core-metrics-phylogenetic \
  --i-table 06_phylogeny/table_filtered.qza \
  --i-phylogeny 06_phylogeny/rooted_tree.qza \
  --p-sampling-depth "$DEPTH" \
  --m-metadata-file 00_metadata/qiime_metadata.tsv \
  --p-n-jobs-or-threads "$NCORES" \
  --output-dir 07_diversity/core_metrics

# Alpha group significance for mets_status
for MET in observed_features shannon; do
  qiime_run qiime diversity alpha-group-significance \
    --i-alpha-diversity 07_diversity/core_metrics/${MET}_vector.qza \
    --m-metadata-file 00_metadata/qiime_metadata.tsv \
    --o-visualization 07_diversity/${MET}_group_significance.qzv || true
done

# Beta: PERMANOVA for mets_status
for DIST in bray_curtis jaccard unweighted_unifrac weighted_unifrac; do
  qiime_run qiime diversity beta-group-significance \
    --i-distance-matrix 07_diversity/core_metrics/${DIST}_distance_matrix.qza \
    --m-metadata-file 00_metadata/qiime_metadata.tsv \
    --m-metadata-column mets_status \
    --p-method permanova \
    --p-permutations 999 \
    --o-visualization 07_diversity/${DIST}_permanova_mets_status.qzv || true
done

# Export distance matrices / ordinations for R
for DIST in bray_curtis jaccard unweighted_unifrac weighted_unifrac; do
  qiime_run qiime tools export \
    --input-path 07_diversity/core_metrics/${DIST}_distance_matrix.qza \
    --output-path 07_diversity/${DIST}_export
done

# Genus collapse for DA / abundance tables
qiime_run qiime taxa collapse \
  --i-table 06_phylogeny/table_filtered.qza \
  --i-taxonomy 05_taxonomy/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table 07_diversity/genus_table.qza

qiime_run qiime tools export \
  --input-path 07_diversity/genus_table.qza \
  --output-path 07_diversity/genus_table_export
qiime_run biom convert \
  -i 07_diversity/genus_table_export/feature-table.biom \
  -o results/genus_abundance.tsv \
  --to-tsv

# Family collapse for sensitivity
qiime_run qiime taxa collapse \
  --i-table 06_phylogeny/table_filtered.qza \
  --i-taxonomy 05_taxonomy/taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table 07_diversity/family_table.qza
qiime_run qiime tools export \
  --input-path 07_diversity/family_table.qza \
  --output-path 07_diversity/family_table_export
qiime_run biom convert \
  -i 07_diversity/family_table_export/feature-table.biom \
  --to-tsv \
  -o results/family_abundance.tsv

echo "Diversity complete"
