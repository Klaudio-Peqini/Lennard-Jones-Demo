#!/usr/bin/env python3
"""
Lennard-Jones molecular dynamics demo for computational chemistry / HPC teaching.

Model: monoatomic argon-like fluid in reduced Lennard-Jones units.
Integrator: velocity-Verlet with periodic boundary conditions and minimum-image convention.
Outputs: thermo.csv and trajectory.xyz.

Reduced units:
    sigma = epsilon = mass = k_B = 1
Typical liquid argon-like state point:
    rho* = 0.8442, T* = 0.728, dt* = 0.005
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class MDState:
    pos: np.ndarray
    vel: np.ndarray
    force: np.ndarray
    box: float
    epot: float
    virial: float


def initialize_positions(n: int, rho: float) -> tuple[np.ndarray, float]:
    """Place particles on a simple cubic lattice inside a cubic periodic box."""
    box = (n / rho) ** (1.0 / 3.0)
    cells = math.ceil(n ** (1.0 / 3.0))
    spacing = box / cells
    coords = []
    for ix in range(cells):
        for iy in range(cells):
            for iz in range(cells):
                if len(coords) < n:
                    coords.append([(ix + 0.5) * spacing, (iy + 0.5) * spacing, (iz + 0.5) * spacing])
    return np.asarray(coords, dtype=np.float64), box


def initialize_velocities(n: int, temperature: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    vel = rng.normal(0.0, math.sqrt(temperature), size=(n, 3))
    vel -= vel.mean(axis=0)  # remove center-of-mass drift
    kin = 0.5 * np.sum(vel * vel)
    current_temp = 2.0 * kin / (3.0 * (n - 1))
    vel *= math.sqrt(temperature / current_temp)
    return vel


def compute_forces(pos: np.ndarray, box: float, cutoff: float) -> tuple[np.ndarray, float, float]:
    """Vectorized all-pairs LJ force with periodic boundaries.

    Returns forces, potential energy, and virial contribution. This is clear and compact,
    suitable for teaching. For very large N, domain decomposition / neighbor lists would be used.
    """
    n = pos.shape[0]
    rij = pos[:, None, :] - pos[None, :, :]
    rij -= box * np.rint(rij / box)
    r2 = np.sum(rij * rij, axis=2)

    mask = np.triu((r2 < cutoff * cutoff) & (r2 > 0.0), k=1)
    pair_rij = rij[mask]
    pair_r2 = r2[mask]

    inv_r2 = 1.0 / pair_r2
    inv_r6 = inv_r2 ** 3
    inv_r12 = inv_r6 * inv_r6

    # LJ potential, shifted so U(rc)=0.
    inv_rc2 = 1.0 / (cutoff * cutoff)
    inv_rc6 = inv_rc2 ** 3
    inv_rc12 = inv_rc6 * inv_rc6
    u_shift = 4.0 * (inv_rc12 - inv_rc6)
    epot = np.sum(4.0 * (inv_r12 - inv_r6) - u_shift)

    # Force on i from j: F_ij = 24(2/r^14 - 1/r^8) r_ij
    fij = (24.0 * inv_r2 * (2.0 * inv_r12 - inv_r6))[:, None] * pair_rij
    forces = np.zeros_like(pos)
    ii, jj = np.where(mask)
    np.add.at(forces, ii, fij)
    np.add.at(forces, jj, -fij)

    # Virial: sum over pairs r_ij · F_ij
    virial = np.sum(np.einsum("ij,ij->i", pair_rij, fij))
    return forces, float(epot), float(virial)


def kinetic_energy(vel: np.ndarray) -> float:
    return float(0.5 * np.sum(vel * vel))


def instantaneous_temperature(ekin: float, n: int) -> float:
    return float(2.0 * ekin / (3.0 * (n - 1)))


def write_xyz(handle, pos: np.ndarray, step: int, box: float, symbol: str = "Ar") -> None:
    handle.write(f"{pos.shape[0]}\n")
    handle.write(f"step={step} box={box:.8f}\n")
    for x, y, z in pos:
        handle.write(f"{symbol} {x:.8f} {y:.8f} {z:.8f}\n")


def run_md(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pos, box = initialize_positions(args.n, args.rho)
    vel = initialize_velocities(args.n, args.temperature, args.seed)
    force, epot, virial = compute_forces(pos, box, args.cutoff)

    thermo_path = outdir / "thermo.csv"
    traj_path = outdir / "trajectory.xyz"
    meta_path = outdir / "metadata.txt"

    start = time.perf_counter()
    with thermo_path.open("w", newline="") as thermofile, traj_path.open("w") as xyzfile:
        writer = csv.writer(thermofile)
        writer.writerow(["step", "time", "ekin", "epot", "etot", "temperature", "pressure", "wall_seconds"])

        for step in range(args.steps + 1):
            ekin = kinetic_energy(vel)
            temp = instantaneous_temperature(ekin, args.n)
            volume = box ** 3
            pressure = (args.n * temp + virial / 3.0) / volume
            elapsed = time.perf_counter() - start

            if step % args.sample_every == 0:
                writer.writerow([step, step * args.dt, ekin, epot, ekin + epot, temp, pressure, elapsed])

            if step % args.traj_every == 0:
                write_xyz(xyzfile, pos, step, box)

            if step == args.steps:
                break

            # Velocity-Verlet integration.
            vel += 0.5 * args.dt * force
            pos += args.dt * vel
            pos %= box
            force, epot, virial = compute_forces(pos, box, args.cutoff)
            vel += 0.5 * args.dt * force

    total = time.perf_counter() - start
    with meta_path.open("w") as f:
        f.write("Lennard-Jones molecular dynamics demo\n")
        f.write("Implementation: Python/NumPy vectorized all-pairs\n")
        f.write(f"N={args.n}\n")
        f.write(f"rho={args.rho}\n")
        f.write(f"temperature={args.temperature}\n")
        f.write(f"dt={args.dt}\n")
        f.write(f"steps={args.steps}\n")
        f.write(f"cutoff={args.cutoff}\n")
        f.write(f"box={box}\n")
        f.write(f"seed={args.seed}\n")
        f.write(f"wall_seconds={total}\n")

    print(f"Done. Results written to: {outdir}")
    print(f"Thermodynamics: {thermo_path}")
    print(f"Trajectory:      {traj_path}")
    print(f"Wall time:       {total:.3f} s")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Lennard-Jones MD demo in reduced units")
    p.add_argument("--n", type=int, default=256, help="number of particles")
    p.add_argument("--rho", type=float, default=0.8442, help="reduced density")
    p.add_argument("--temperature", type=float, default=0.728, help="initial reduced temperature")
    p.add_argument("--dt", type=float, default=0.005, help="reduced time step")
    p.add_argument("--steps", type=int, default=1000, help="number of MD steps")
    p.add_argument("--cutoff", type=float, default=2.5, help="LJ cutoff radius")
    p.add_argument("--sample-every", type=int, default=10, help="write thermodynamics every this many steps")
    p.add_argument("--traj-every", type=int, default=100, help="write XYZ frame every this many steps")
    p.add_argument("--seed", type=int, default=12345, help="random seed for velocities")
    p.add_argument("--outdir", type=str, default="results/python_lj", help="output directory")
    return p.parse_args()


if __name__ == "__main__":
    run_md(parse_args())
