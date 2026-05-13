#!/usr/bin/env python3
"""Compute and plot the radial distribution function g(r) from the last frame of trajectory.xyz."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_last_xyz_frame(path: Path):
    lines = path.read_text().splitlines()
    i = 0
    last_symbols = []
    last_pos = None
    last_box = None
    while i < len(lines):
        n = int(lines[i].strip())
        comment = lines[i + 1]
        box = None
        for item in comment.split():
            if item.startswith("box="):
                box = float(item.split("=", 1)[1])
        data = lines[i + 2 : i + 2 + n]
        pos = []
        symbols = []
        for row in data:
            s, x, y, z = row.split()[:4]
            symbols.append(s)
            pos.append([float(x), float(y), float(z)])
        last_symbols = symbols
        last_pos = np.asarray(pos)
        last_box = box
        i += n + 2
    if last_pos is None or last_box is None:
        raise ValueError("No complete XYZ frame with box metadata was found.")
    return last_symbols, last_pos, last_box


def compute_rdf(pos: np.ndarray, box: float, bins: int, rmax: float | None = None):
    n = len(pos)
    if rmax is None:
        rmax = box / 2.0
    distances = []
    for i in range(n - 1):
        rij = pos[i + 1 :] - pos[i]
        rij -= box * np.rint(rij / box)
        distances.extend(np.linalg.norm(rij, axis=1))
    distances = np.asarray(distances)
    hist, edges = np.histogram(distances, bins=bins, range=(0.0, rmax))
    r = 0.5 * (edges[:-1] + edges[1:])
    dr = edges[1] - edges[0]
    rho = n / box**3
    shell_volume = 4.0 * np.pi * r**2 * dr
    ideal_pairs = 0.5 * n * rho * shell_volume
    g = hist / ideal_pairs
    return r, g


def main() -> None:
    p = argparse.ArgumentParser(description="Compute RDF g(r) from an LJ trajectory")
    p.add_argument("--xyz", required=True, help="trajectory.xyz")
    p.add_argument("--bins", type=int, default=100)
    p.add_argument("--outdir", default=None)
    args = p.parse_args()

    xyz = Path(args.xyz)
    outdir = Path(args.outdir) if args.outdir else xyz.parent
    outdir.mkdir(parents=True, exist_ok=True)

    _, pos, box = read_last_xyz_frame(xyz)
    r, g = compute_rdf(pos, box, args.bins)
    df = pd.DataFrame({"r": r, "g_r": g})
    df.to_csv(outdir / "rdf.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(r, g, linewidth=1.8)
    ax.set_xlabel("r / σ")
    ax.set_ylabel("g(r)")
    ax.set_title("Radial distribution function from final MD frame")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "rdf.png", dpi=180)
    plt.close(fig)

    print(f"RDF written to: {outdir / 'rdf.csv'}")
    print(f"Plot written to: {outdir / 'rdf.png'}")


if __name__ == "__main__":
    main()
