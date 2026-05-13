#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

make -C cpp
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
./cpp/lj_md \
  --n 512 \
  --rho 0.8442 \
  --temperature 0.728 \
  --dt 0.005 \
  --steps 1000 \
  --sample-every 10 \
  --traj-every 100 \
  --outdir "results/cpp_N512_T${OMP_NUM_THREADS}"

python3 python/plot_md.py --input "results/cpp_N512_T${OMP_NUM_THREADS}/thermo.csv" --label "C++/OpenMP, N=512, threads=${OMP_NUM_THREADS}"
python3 python/rdf_from_xyz.py --xyz "results/cpp_N512_T${OMP_NUM_THREADS}/trajectory.xyz"
