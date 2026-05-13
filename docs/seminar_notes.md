# Seminar notes: using the HPC cluster for computational chemistry

## Demonstration title

**From an interatomic potential to HPC molecular dynamics: a Python and C++/OpenMP demonstration**

## Core message for the chemistry audience

The cluster is not only a place where code runs faster. It is a place where one can systematically explore molecular models: different temperatures, densities, compositions, force-field parameters, and trajectory lengths. A single desktop run gives one trajectory; an HPC cluster turns the problem into a reproducible computational experiment.

## What is being simulated?

A simple argon-like fluid. Each atom interacts with all other atoms inside a cutoff radius through the Lennard-Jones potential. The repulsive part prevents atoms from collapsing into each other, while the attractive part models dispersion-like interactions.

## Equations shown during the demo

Potential:

```text
U(r) = 4 epsilon [ (sigma/r)^12 - (sigma/r)^6 ]
```

Force:

```text
F_ij = 24 epsilon / r_ij^2 [ 2 (sigma/r_ij)^12 - (sigma/r_ij)^6 ] r_ij
```

Velocity-Verlet algorithm:

```text
v(t + dt/2) = v(t) + (dt/2) F(t)/m
r(t + dt)   = r(t) + dt v(t + dt/2)
F(t + dt)   = force from new positions
v(t + dt)   = v(t + dt/2) + (dt/2) F(t + dt)/m
```

Temperature from kinetic energy:

```text
T = 2 K / [3 (N - 1) k_B]
```

Pressure from the virial expression:

```text
P = [N k_B T + W/3] / V
```

## What students should observe

- The total energy is approximately conserved in an NVE trajectory.
- Kinetic and potential energies fluctuate, but their sum is much more stable.
- The radial distribution function has peaks, showing short-range structure.
- Increasing particle number makes the computational cost grow rapidly.
- The OpenMP implementation can use several CPU cores, but the speedup is not perfectly linear.

## Suggested live commands

```bash
python3 -m venv chem-md-env
source chem-md-env/bin/activate
pip install -r requirements.txt
bash scripts/run_python_demo.sh
make -C cpp
OMP_NUM_THREADS=8 bash scripts/run_cpp_demo.sh
THREAD_LIST="1 2 4 8" N=512 STEPS=300 bash scripts/benchmark.sh
```

## Discussion prompts

- Why is the force calculation the expensive part?
- Why does the Python version help us understand the algorithm but the C++ version helps us use the cluster?
- What would change if the atoms were charged?
- Why are neighbor lists important in production molecular dynamics?
- How would we design a parameter scan over temperature and density?
