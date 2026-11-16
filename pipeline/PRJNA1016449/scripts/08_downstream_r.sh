#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh

# Use QIIME2 env R (has ANCOMBC + phyloseq). MaAsLin2 skipped — one DA method is enough.
export PATH="$QIIME_ENV_PREFIX/bin:$PATH"
export R_LIBS_USER=
export R_LIBS_SITE=
unset R_HOME || true

echo "Rscript=$(which Rscript)"
Rscript -e 'stopifnot(requireNamespace("ANCOMBC", quietly=TRUE)); cat("ANCOMBC", as.character(packageVersion("ANCOMBC")), "\n")'

Rscript scripts/downstream_analysis.R 2>&1 | tee logs/downstream_analysis.log
echo "Downstream R analysis complete (ANCOM-BC2 only)"
