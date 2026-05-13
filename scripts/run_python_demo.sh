#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

python3 python/lj_md.py \
  --n 256 \
  --rho 0.8442 \
  --temperature 0.728 \
  --dt 0.005 \
  --steps 1000 \
  --sample-every 10 \
  --traj-every 100 \
  --outdir results/python_N256

python3 python/plot_md.py --input results/python_N256/thermo.csv --label "Python/NumPy, N=256"
python3 python/rdf_from_xyz.py --xyz results/python_N256/trajectory.xyz
