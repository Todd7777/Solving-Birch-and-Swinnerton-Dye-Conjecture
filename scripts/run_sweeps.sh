#!/usr/bin/env bash
# Run (or resume) the numerical BSD sweeps over Cremona's database for all conductors < 20000.
# Idempotent: every stage resumes from its JSONL file, so it is safe to re-run after an interruption.
#   WORKERS=2 scripts/run_sweeps.sh
set -u
cd "$(dirname "$0")/.."
W=${WORKERS:-2}
run_range () {  # $1 = min, $2 = max, $3 = output stem
  # wait for any other sweep on the same output to finish first
  while pgrep -f "bsd sweep .*--out $3( |$)" > /dev/null; do sleep 30; done
  python3 -m bsd sweep --min "$1" --max "$2" --workers "$W" --out "$3" --resume >> "$3.log" 2>&1
  python3 scripts/repair_records.py "$3.jsonl" --workers "$W" >> "$3.log" 2>&1
  python3 scripts/compact_results.py "$3.jsonl" >> "$3.log" 2>&1
}
run_range 1 1000 results/sweep_0001_1000
run_range 1001 10000 results/sweep_1001_10000
run_range 10001 19999 results/sweep_10001_19999
echo "all sweeps done"
