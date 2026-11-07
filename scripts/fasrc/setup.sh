#!/usr/bin/env bash
# One-time setup on a FASRC login node (needs internet access; do NOT run heavy work on login nodes).
#   bash scripts/fasrc/setup.sh            # core toolkit (PARI via cypari2, python-flint)
#   WITH_SAGE=1 bash scripts/fasrc/setup.sh  # also Sage (passagemath), needed only for prove_array.sbatch
set -euo pipefail
cd "$(dirname "$0")/../.."
# Load a Python >= 3.10.  Module names change over time: `module avail python` lists them.
module load python 2>/dev/null || true
python3 -c 'import sys; assert sys.version_info >= (3, 10), sys.version' || { echo "need Python >= 3.10 (try: module avail python)"; exit 1; }
VENV=${BSD_VENV:-$HOME/bsd-venv}
[ -d "$VENV" ] || python3 -m venv "$VENV"
source "$VENV/bin/activate"
pip install --upgrade pip
pip install -r requirements.txt
if [ "${WITH_SAGE:-0}" = "1" ]; then pip install passagemath-standard; fi
# Download all of Cremona's database blocks (conductors 1..499999, ~650 MB) so compute nodes need no internet.
python -m bsd fetch --min 0 --max 499999
python -m pytest -q tests/test_known_values.py
echo "setup done: activate with  source $VENV/bin/activate"
