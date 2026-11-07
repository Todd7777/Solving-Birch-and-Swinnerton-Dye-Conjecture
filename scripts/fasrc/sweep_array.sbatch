#!/usr/bin/env bash
# Numerical verification of the BSD formula for every curve in Cremona's database (conductors < 500000),
# one 10000-conductor block per array task.  Resumable: a requeued/preempted task continues where it stopped.
#   sbatch scripts/fasrc/sweep_array.sbatch                 # all 50 blocks
#   sbatch --array=0 --export=ALL,LIMIT=2000 scripts/fasrc/sweep_array.sbatch   # quick test
#SBATCH -J bsd-sweep
#SBATCH -p shared,serial_requeue
#SBATCH --requeue
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 0-06:00
#SBATCH --array=0-49
#SBATCH -o results/fasrc/logs/sweep_%a_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
source "${BSD_VENV:-$HOME/bsd-venv}/bin/activate"
mkdir -p results/fasrc/logs
ID=${SLURM_ARRAY_TASK_ID:-0}
LO=$(( ID * 10000 )); HI=$(( LO + 9999 )); [ "$LO" -eq 0 ] && LO=1
EXTRA=""; [ -n "${LIMIT:-}" ] && EXTRA="--limit $LIMIT"
OUT=$(printf "results/fasrc/sweep_%06d_%06d" "$LO" "$HI")
echo "task $ID: conductors $LO..$HI on $(hostname) with ${SLURM_CPUS_PER_TASK:-1} workers"
python -m bsd sweep --min "$LO" --max "$HI" --workers "${SLURM_CPUS_PER_TASK:-1}" \
    --compact --no-modsym --stack-gb 2 --resume --out "$OUT" $EXTRA
echo "task $ID finished: $(wc -l < "$OUT.jsonl") rows"
