#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/env.sh

export PATH="$QIIME_ENV_PREFIX/bin:$PATH"
export R_LIBS_USER=
export R_LIBS_SITE=
unset R_HOME || true

echo "Rscript=$(which Rscript)"
Rscript -e 'stopifnot(requireNamespace("ANCOMBC", quietly=TRUE)); cat("ANCOMBC", as.character(packageVersion("ANCOMBC")), "\n")'

Rscript scripts/run_ancombc2_only.R 2>&1 | tee logs/ancombc2_only.log
Rscript scripts/finish_pipeline.R 2>&1 | tee logs/finish_pipeline.log
echo "Downstream R analysis complete (ANCOM-BC2 only)"
