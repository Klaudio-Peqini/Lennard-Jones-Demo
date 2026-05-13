#!/usr/bin/env python3
"""Plot thermodynamic diagnostics from the Lennard-Jones MD demo."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_column(df: pd.DataFrame, x: str, y: str, ylabel: str, outpath: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(df[x], df[y], linewidth=1.8)
    ax.set_xlabel("Reduced time")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description="Plot LJ-MD thermodynamics")
    p.add_argument("--input", type=str, required=True, help="path to thermo.csv")
    p.add_argument("--outdir", type=str, default=None, help="directory for plots; default: same as input")
    p.add_argument("--label", type=str, default="Lennard-Jones MD", help="plot title prefix")
    args = p.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir) if args.outdir else input_path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path)
    plot_column(df, "time", "etot", "Total energy", outdir / "total_energy.png", f"{args.label}: energy conservation")
    plot_column(df, "time", "temperature", "Temperature", outdir / "temperature.png", f"{args.label}: instantaneous temperature")
    plot_column(df, "time", "pressure", "Pressure", outdir / "pressure.png", f"{args.label}: pressure")

    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(df["time"], df["ekin"], label="Kinetic")
    ax.plot(df["time"], df["epot"], label="Potential")
    ax.plot(df["time"], df["etot"], label="Total")
    ax.set_xlabel("Reduced time")
    ax.set_ylabel("Energy")
    ax.set_title(f"{args.label}: energy components")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "energy_components.png", dpi=180)
    plt.close(fig)

    print(f"Plots written to: {outdir}")


if __name__ == "__main__":
    main()
