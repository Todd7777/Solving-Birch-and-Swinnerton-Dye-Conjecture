# Running the heavy computations on FASRC (Harvard's Cannon cluster)

The toolkit's sweeps are embarrassingly parallel, so a cluster turns days into minutes. A single
16-core task verifies one 10,000-conductor block of Cremona's database, about 60,000 curves, in roughly 10–30
minutes. The 50-task job array below covers **all ~3.06 million curves of conductor < 500,000**.

These steps run under *your* FASRC account. The session that produced this repository cannot SSH
anywhere: it has HTTPS-only network access and must never hold your password or two-factor codes.

## 1. One-time setup (login node; uses the network, light CPU)

```bash
ssh <your-username>@login.rc.fas.harvard.edu
git clone -b claude/magical-cannon-cfb1sm https://github.com/Todd7777/Solving-Birch-and-Swinnerton-Dye-Conjecture.git
cd Solving-Birch-and-Swinnerton-Dye-Conjecture
bash scripts/fasrc/setup.sh            # venv in ~/bsd-venv, pip install, download Cremona's data, quick tests
```

If `setup.sh` cannot find Python ≥ 3.10, run `module avail python`, `module load` a suitable one, and rerun.

## 2. Smoke test (one block, 2,000 curves, ~5 minutes)

```bash
sbatch --array=0 --export=ALL,LIMIT=2000 scripts/fasrc/sweep_array.sbatch
squeue -u $USER                                     # watch it
tail -f results/fasrc/logs/sweep_0_*.out
```

## 3. The full database (50 tasks × 16 cores)

```bash
sbatch scripts/fasrc/sweep_array.sbatch
```

Tasks are resumable and `--requeue`able, so the preemptible `serial_requeue` partition is fine. Resubmitting
continues where a task stopped. Edit the `#SBATCH -p/-c/--mem/-t` lines to match your allocation.

## 4. Summarise and send the results back

```bash
python scripts/fasrc/aggregate.py sweep   # writes results/fasrc/sweep_summary.{json,md}, sweep_anomalies.jsonl, sweep_rank_ge2.csv.gz
git add results/fasrc/sweep_summary.* results/fasrc/sweep_anomalies.jsonl results/fasrc/sweep_rank_ge2.csv.gz
git commit -m "FASRC sweep of Cremona's database (conductor < 500000)" && git push
```

The per-curve files (`results/fasrc/sweep_*.jsonl`, ~2 GB in total) are ignored by git. Keep them on the cluster.

## Optional: proofs of the full BSD formula with Sage

`prove_array.sbatch` runs Sage's `prove_BSD` (Kolyvagin and Kato bounds, 2-descent) plus provable saturation on every
curve of analytic rank ≤ 1, in 1,000-conductor slices (default: conductors < 20,000). It needs
`WITH_SAGE=1 bash scripts/fasrc/setup.sh` (a large install). Summarise with `python scripts/fasrc/aggregate.py prove`.

## What the sweep certifies, per curve

* Numerical verification of the full BSD formula: $\#Ш_{\rm an}$ is an integer and a perfect square. The output also records
  agreement with Cremona's independently computed values (eclib/mwrank vs PARI here).
* The rank part of BSD:
  * as a **theorem** for analytic rank ≤ 1 (Gross–Zagier–Kolyvagin; numerical non-vanishing of the
    first non-zero derivative);
  * **verified curve by curve** for analytic rank 2 and 3: the 2-descent upper bound is met by independent points, and
    Gross–Zagier–Kolyvagin plus parity give the matching lower bound for $\operatorname{ord}_{s=1}L(E,s)$;
  * **open** for analytic rank ≥ 4, where no method proves the vanishing of $L''(E,1)$.
* This mode skips modular symbols (`--no-modsym`). For root number $+1$ and rank 0 the non-vanishing of $L(E,1)$ is
  therefore numerical, and the status line says so.
