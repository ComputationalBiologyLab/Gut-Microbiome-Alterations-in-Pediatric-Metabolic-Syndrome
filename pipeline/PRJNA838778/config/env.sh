#!/usr/bin/env bash
# Shared environment activation helpers
# Locate the conda installation (override with MINICONDA_ROOT).
export MINICONDA_ROOT="${MINICONDA_ROOT:-$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")}"
export CONDA_EXE="${CONDA_EXE:-$MINICONDA_ROOT/bin/conda}"
export QIIME_ENV="${QIIME_ENV:-qiime2-amplicon-2024.10}"
export R_ENV="${R_ENV:-bio}"
export NCORES="${NCORES:-4}"

export QIIME_ENV_PREFIX="${QIIME_ENV_PREFIX:-$MINICONDA_ROOT/envs/$QIIME_ENV}"
export R_ENV_PREFIX="${R_ENV_PREFIX:-$MINICONDA_ROOT/envs/$R_ENV}"

# run_dada.R uses `#!/usr/bin/env Rscript`: make sure the conda environment's R is used
# and not a system-wide R.
export PATH="$QIIME_ENV_PREFIX/bin:$PATH"
export R_LIBS_USER=""
export R_LIBS_SITE=""
unset R_HOME || true

# Call tools from the environment prefix directly (avoids the PATH reset of `conda run`)
qiime_run() {
  env \
    PATH="$QIIME_ENV_PREFIX/bin:$PATH" \
    R_LIBS_USER= \
    R_LIBS_SITE= \
    R_HOME="$QIIME_ENV_PREFIX/lib/R" \
    "$@"
}

r_run() {
  env PATH="$R_ENV_PREFIX/bin:$PATH" "$@"
}
