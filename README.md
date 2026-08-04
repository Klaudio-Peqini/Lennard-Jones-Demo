# HPC Computational Chemistry Demo: Lennard-Jones Molecular Dynamics          

This package is a compact seminar demonstration for a chemistry department on how an HPC cluster can be used for molecular simulation. It contains the same physical model implemented in **Python/NumPy** and **C++/OpenMP**, plus plotting, benchmarking, and HTCondor submission scripts.

## 1. Scientific case

The model is a monoatomic argon-like fluid described by the Lennard-Jones potential

$$
U(r)=4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^6\right].
$$

The simulation uses reduced Lennard-Jones units:

$$
\sigma=\varepsilon=m=k_B=1.
$$

The default state point is liquid-like:

- reduced density: `rho = 0.8442`
- reduced temperature: `T = 0.728`
- time step: `dt = 0.005`
- cutoff: `rc = 2.5 sigma`

This is a useful teaching case because it connects directly to computational chemistry ideas:

- interatomic potential energy surface;
- forces from the gradient of the potential;
- molecular dynamics integration;
- energy conservation;
- temperature and pressure from microscopic motion;
- radial distribution function `g(r)` as a structural observable;
- HPC scaling by increasing particle number and thread count.

## 2. Package structure

```text
hpc_chem_md_demo/
├── python/
│   ├── lj_md.py             # Python/NumPy MD implementation
│   ├── plot_md.py           # Thermodynamic plots
│   └── rdf_from_xyz.py      # Radial distribution function from final frame
├── cpp/
│   ├── lj_md.cpp            # C++17/OpenMP MD implementation
│   └── Makefile             # Build file
├── scripts/
│   ├── run_python_demo.sh   # One-command Python run + plots
│   ├── run_cpp_demo.sh      # One-command C++ run + plots
│   ├── benchmark.sh         # OpenMP scaling benchmark
│   ├── plot_scaling.py      # Scaling plot
│   ├── htcondor_python.sub  # HTCondor submit file for Python
│   └── htcondor_cpp.sub     # HTCondor submit file for C++/OpenMP
├── configs/
│   └── argon_256.json       # Human-readable reference parameters
├── results/                 # Outputs are written here
└── requirements.txt
```

## 3. Prepare the Python environment

On the HPC login node:

```bash
python3 -m venv chem-md-env
source chem-md-env/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

## 4. Run the Python demonstration

```bash
bash scripts/run_python_demo.sh
```

This produces:

```text
results/python_N256/thermo.csv
results/python_N256/trajectory.xyz
results/python_N256/total_energy.png
results/python_N256/temperature.png
results/python_N256/pressure.png
results/python_N256/energy_components.png
results/python_N256/rdf.csv
results/python_N256/rdf.png
```

For a faster live run:

```bash
python3 python/lj_md.py --n 108 --steps 300 --outdir results/python_quick
python3 python/plot_md.py --input results/python_quick/thermo.csv
python3 python/rdf_from_xyz.py --xyz results/python_quick/trajectory.xyz
```

## 5. Run the C++/OpenMP demonstration

Build:

```bash
make -C cpp
```

Run with 8 OpenMP threads:

```bash
export OMP_NUM_THREADS=8
./cpp/lj_md --n 512 --steps 1000 --outdir results/cpp_N512_T8
python3 python/plot_md.py --input results/cpp_N512_T8/thermo.csv --label "C++/OpenMP, N=512, threads=8"
python3 python/rdf_from_xyz.py --xyz results/cpp_N512_T8/trajectory.xyz
```

Or run the complete wrapper:

```bash
OMP_NUM_THREADS=8 bash scripts/run_cpp_demo.sh
```

## 6. Run a small scaling benchmark

```bash
THREAD_LIST="1 2 4 8 16" N=1024 STEPS=500 bash scripts/benchmark.sh
```

The benchmark produces:

```text
results/benchmark/cpp_openmp_scaling_N1024_steps500.csv
results/benchmark/cpp_openmp_scaling_N1024_steps500.png
```

This is the most useful part for the HPC seminar: keep `N` and `steps` fixed, change the number of threads, and discuss speedup, saturation, memory traffic, and parallel overhead.


## 7. Run a small temperature-density sweep

This demonstrates the real HPC use case: not one simulation, but a small computational campaign over several thermodynamic state points.

```bash
OMP_NUM_THREADS=8 N=512 STEPS=800 TEMPS="0.6 0.8 1.0" RHOS="0.70 0.85" bash scripts/sweep_temperature_density.sh
```

The sweep writes one folder per state point and a compact summary table:

```text
results/sweep_N512_T8/sweep_summary.csv
```

Use this in the seminar to explain why job arrays and parameter scans are natural in computational chemistry.

## 8. Submit to HTCondor

From the package root:

```bash
condor_submit scripts/htcondor_python.sub
condor_submit scripts/htcondor_cpp.sub
```

Check jobs:

```bash
condor_q
```

After completion:

```bash
ls -lh results/
```

The submit files are intentionally minimal. Adjust `request_cpus`, `request_memory`, and the wrapper scripts according to your local cluster policy.

## 9. Suggested live seminar flow

1. **Physical problem**: atoms interact through an empirical potential.
2. **Numerical method**: velocity-Verlet integration of Newton's equations.
3. **Python run**: easy to read, good for pedagogy, produces scientific observables.
4. **C++/OpenMP run**: same physics, much better for HPC scaling.
5. **Output interpretation**: total energy should remain nearly conserved, temperature fluctuates, and `g(r)` reveals liquid-like short-range order.
6. **HPC message**: larger `N`, longer trajectories, multiple state points, and parameter scans are where the cluster becomes necessary.

## 10. Limitations and natural extensions

This package deliberately uses an all-pairs cutoff calculation for clarity. Production molecular dynamics codes normally add:

- neighbor lists / cell lists;
- long-range electrostatics for charged or polar systems;
- thermostats and barostats;
- molecular constraints;
- multi-component force fields;
- MPI domain decomposition;
- GPU acceleration.

For chemistry students, good next extensions are:

- replacing the monoatomic fluid with a binary mixture;
- adding a thermostat to sample the canonical ensemble;
- computing diffusion from the mean-square displacement;
- computing a simple vibrational spectrum from velocity autocorrelation;
- comparing Python, C++, OpenMP, MPI, and GPU versions.
