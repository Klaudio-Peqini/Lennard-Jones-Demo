#!/usr/bin/env bash
set -euo pipefail

# Small parameter sweep that mimics a computational chemistry campaign.
# Adjust lists and resource sizes for the real cluster.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
make -C cpp

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
N="${N:-512}"
STEPS="${STEPS:-800}"
TEMPS="${TEMPS:-0.6 0.8 1.0}"
RHOS="${RHOS:-0.70 0.85}"
BASE="results/sweep_N${N}_T${OMP_NUM_THREADS}"
mkdir -p "$BASE"

for RHO in $RHOS; do
  for TEMP in $TEMPS; do
    OUT="$BASE/rho_${RHO}_temp_${TEMP}"
    echo "Running rho=${RHO}, T=${TEMP}, N=${N}, threads=${OMP_NUM_THREADS}"
    ./cpp/lj_md --n "$N" --rho "$RHO" --temperature "$TEMP" --steps "$STEPS" \
      --sample-every 20 --traj-every "$STEPS" --outdir "$OUT" >/dev/null
  done
done

python3 scripts/summarize_sweep.py --base "$BASE"
echo "Sweep summary: $BASE/sweep_summary.csv"
