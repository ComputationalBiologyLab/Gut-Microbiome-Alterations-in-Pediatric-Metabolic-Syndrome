#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh
mkdir -p 05_taxonomy ref

CLASSIFIER="ref/silva-138-99-nb-classifier.qza"
if [[ ! -s "$CLASSIFIER" ]]; then
  echo "Downloading SILVA 138 NB classifier (sklearn 1.4.2)..."
  curl -L --fail --retry 5 -o "$CLASSIFIER" \
    "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
fi

qiime_run qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads 04_dada2/rep_seqs.qza \
  --p-confidence 0.7 \
  --p-n-jobs "$NCORES" \
  --o-classification 05_taxonomy/taxonomy.qza

qiime_run qiime metadata tabulate \
  --m-input-file 05_taxonomy/taxonomy.qza \
  --o-visualization 05_taxonomy/taxonomy.qzv

qiime_run qiime taxa barplot \
  --i-table 04_dada2/table.qza \
  --i-taxonomy 05_taxonomy/taxonomy.qza \
  --m-metadata-file 00_metadata/qiime_metadata.tsv \
  --o-visualization 05_taxonomy/taxa_barplot.qzv

# Export taxonomy TSV
qiime_run qiime tools export \
  --input-path 05_taxonomy/taxonomy.qza \
  --output-path 05_taxonomy/taxonomy_export
cp 05_taxonomy/taxonomy_export/taxonomy.tsv results/ASV_taxonomy.tsv

cat > 05_taxonomy/taxonomy_parameters.txt <<EOF
database=SILVA
version=138
identity=99%
classifier=Naive Bayes q2-feature-classifier classify-sklearn
classifier_artifact=${CLASSIFIER}
confidence=0.7
EOF

echo "Taxonomy complete"
