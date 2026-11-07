#!/usr/bin/env bash
# OPTIONAL: prove the full BSD formula (Sage prove_BSD: Kolyvagin/Kato bounds, 2-descent) for curves of
# analytic rank <= 1, with provable saturation of Mordell-Weil bases.  Needs Sage: WITH_SAGE=1 bash scripts/fasrc/setup.sh
# One task per 1000-conductor slice; default slices cover conductors 1..19999 (array 0-19).
#   sbatch scripts/fasrc/prove_array.sbatch
#SBATCH -J bsd-prove
#SBATCH -p shared,serial_requeue
#SBATCH --requeue
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 0-12:00
#SBATCH --array=0-19
#SBATCH -o results/fasrc/logs/prove_%a_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
source "${BSD_VENV:-$HOME/bsd-venv}/bin/activate"
mkdir -p results/fasrc/logs
ID=${SLURM_ARRAY_TASK_ID:-0}
LO=$(( ID * 1000 )); HI=$(( LO + 999 )); [ "$LO" -eq 0 ] && LO=1
OUT=$(printf "results/fasrc/prove_%06d_%06d" "$LO" "$HI")
python -m bsd sweep --min "$LO" --max "$HI" --workers "${SLURM_CPUS_PER_TASK:-1}" \
    --sage --sage-timeout 600 --compact --stack-gb 2 --resume --out "$OUT"
