#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

p = argparse.ArgumentParser(description="Plot OpenMP scaling from benchmark.csv")
p.add_argument("--csv", required=True)
args = p.parse_args()
path = Path(args.csv)
df = pd.read_csv(path)
df["speedup"] = df["wall_seconds"].iloc[0] / df["wall_seconds"]
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(df["threads"], df["speedup"], marker="o", label="Measured speedup")
ax.plot(df["threads"], df["threads"], linestyle="--", label="Ideal speedup")
ax.set_xlabel("OpenMP threads")
ax.set_ylabel("Speedup relative to 1 thread")
ax.set_title("C++/OpenMP Lennard-Jones MD scaling")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()
out = path.with_suffix(".png")
fig.savefig(out, dpi=180)
print(f"Scaling plot written to: {out}")
