#!/usr/bin/env bash
set -euo pipefail

ENV_TARBALL="${ENV_TARBALL:-pymc_portable.tar.gz}"

# Unpack the portable env into the job sandbox
mkdir -p env
tar -xzf "$ENV_TARBALL" -C env

# Use the env's interpreter directly (no 'activate' needed)
if [[ -x env/bin/python ]]; then
  PY="env/bin/python"
  [[ -x env/bin/conda-unpack ]] && env/bin/conda-unpack || true
elif [[ -x env/pymc-portable/bin/python ]]; then
  PY="env/pymc-portable/bin/python"
  [[ -x env/pymc-portable/bin/conda-unpack ]] && env/pymc-portable/bin/conda-unpack || true
else
  echo "ERROR: Could not find Python inside extracted env." >&2
  ls -R env >&2 || true
  exit 2
fi
# 👉 Disable C compilation in PyTensor (portable; no ld needed)
export PYTENSOR_FLAGS='cxx=,optimizer=fast_compile,on_opt_error=warn'
# ── IMPORTANT: give PyTensor a home/cache inside the sandbox ──
export HOME="${HOME:-$PWD}"                    # ensure a usable HOME
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$PWD/.cache}"
#export PYTENSOR_FLAGS="base_compiledir=$PWD/.pytensor_compiledir"  # where to JIT-compile
mkdir -p "$PWD/.pytensor_compiledir" "$XDG_CACHE_HOME"

# Headless matplotlib + sane BLAS threading
export MPLBACKEND=Agg
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

# Outputs
mkdir -p out

# (Optional) sanity check: will now succeed without HOME errors
$PY - <<'PY'
import pymc as pm, numpy as np, pynmrstar
print("PyMC:", pm.__version__)
PY

# Run your job
$PY  lacs.py $1 --data-id $2 --out /home/nmrbox/kbaskaran/lacs/bayes/$2 --method bayes
$PY  lacs.py $1 --data-id $2 --out /home/nmrbox/kbaskaran/lacs/bayes/$2 --method theilsen
$PY  lacs.py $1 --data-id $2 --out /home/nmrbox/kbaskaran/lacs/bayes/$2 --method ransac
$PY  lacs.py $1 --data-id $2 --out /home/nmrbox/kbaskaran/lacs/bayes/$2 --method quantile
$PY  lacs.py $1 --data-id $2 --out /home/nmrbox/kbaskaran/lacs/bayes/$2 --method tukey
