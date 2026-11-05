#!/usr/bin/env bash
# Chain the systematic sweeps over Cremona's database (conductor ranges),
# each waiting for the previous one to finish so the machine is not oversubscribed.
set -u
cd "$(dirname "$0")/.."
W=${WORKERS:-3}
while pgrep -f "bsd sweep --min 1 --max 1000" > /dev/null; do sleep 20; done
python3 -m bsd sweep --min 1001 --max 10000 --workers $W --out results/sweep_1001_10000 --resume > results/sweep_1001_10000.log 2>&1
python3 -m bsd sweep --min 10001 --max 19999 --workers $W --out results/sweep_10001_19999 --resume > results/sweep_10001_19999.log 2>&1
