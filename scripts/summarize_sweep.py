#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

p = argparse.ArgumentParser(description="Summarize LJ-MD parameter sweep outputs")
p.add_argument("--base", required=True, help="base directory containing rho_*_temp_* folders")
args = p.parse_args()
base = Path(args.base)
rows = []
pattern = re.compile(r"rho_([0-9.]+)_temp_([0-9.]+)")
for case_dir in sorted(base.glob("rho_*_temp_*")):
    thermo = case_dir / "thermo.csv"
    meta = case_dir / "metadata.txt"
    if not thermo.exists():
        continue
    m = pattern.match(case_dir.name)
    rho = float(m.group(1)) if m else float("nan")
    temp0 = float(m.group(2)) if m else float("nan")
    df = pd.read_csv(thermo)
    tail = df.iloc[len(df)//2:]
    wall = None
    if meta.exists():
        for line in meta.read_text().splitlines():
            if line.startswith("wall_seconds="):
                wall = float(line.split("=", 1)[1])
    rows.append({
        "case": case_dir.name,
        "rho_input": rho,
        "temperature_input": temp0,
        "temperature_mean_second_half": tail["temperature"].mean(),
        "pressure_mean_second_half": tail["pressure"].mean(),
        "total_energy_mean_second_half": tail["etot"].mean(),
        "total_energy_std_second_half": tail["etot"].std(),
        "wall_seconds": wall,
    })
summary = pd.DataFrame(rows)
out = base / "sweep_summary.csv"
summary.to_csv(out, index=False)
print(summary.to_string(index=False))
print(f"Written: {out}")
