#!/bin/zsh
# DiffractoKit References: the orix oracle for orientation symmetry (IPF sectors/colour key, disorientation).
# orix 0.15.0 from PyPI into its own venv (Python 3.12, the same interpreter as the py4DSTEM venv). MTEX is MATLAB: not used.
# Run:  ../../tools/memcap.sh -- ./provision-orix.sh   (from this folder). Writes orix-venv/, orix-venv-freeze.txt.
set -euo pipefail
REF="$(cd "$(dirname "$0")" && pwd)"
VENV="$REF/orix-venv"
PY="$HOME/miniconda3/envs/py4dstem/bin/python3.12"
[[ -x "$VENV/bin/python" ]] || "$PY" -m venv "$VENV"
"$VENV/bin/python" -m pip install -q --upgrade pip
"$VENV/bin/python" -m pip install -q "orix==0.15.0"
"$VENV/bin/python" -m pip freeze > "$REF/orix-venv-freeze.txt"
"$VENV/bin/python" -I -B -c 'import orix, numpy, scipy; print("orix", orix.__version__, "numpy", numpy.__version__, "scipy", scipy.__version__)'
echo "provision-orix: ready"
