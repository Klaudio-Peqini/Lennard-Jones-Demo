#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
make -C cpp
mkdir -p results/benchmark

N="${N:-1024}"
STEPS="${STEPS:-500}"
THREAD_LIST="${THREAD_LIST:-1 2 4 8 16}"
CSV="results/benchmark/cpp_openmp_scaling_N${N}_steps${STEPS}.csv"
echo "threads,n,steps,wall_seconds" > "$CSV"

for T in $THREAD_LIST; do
  export OMP_NUM_THREADS="$T"
  OUT="results/benchmark/cpp_N${N}_T${T}"
  ./cpp/lj_md --n "$N" --steps "$STEPS" --sample-every "$STEPS" --traj-every "$STEPS" --outdir "$OUT" >/dev/null
  WALL=$(awk -F= '/wall_seconds/ {print $2}' "$OUT/metadata.txt")
  echo "$T,$N,$STEPS,$WALL" | tee -a "$CSV"
done

python3 scripts/plot_scaling.py --csv "$CSV"
echo "Benchmark table: $CSV"
