/*
 Lennard-Jones molecular dynamics demo for computational chemistry / HPC teaching.

 Model: monoatomic argon-like fluid in reduced Lennard-Jones units.
 Integrator: velocity-Verlet with periodic boundary conditions and minimum-image convention.
 Parallelism: OpenMP over pair interactions, using per-thread force buffers.
 Outputs: thermo.csv and trajectory.xyz.

 Build:
     make
 Run:
     OMP_NUM_THREADS=8 ./lj_md --n 512 --steps 2000 --outdir ../results/cpp_N512_T8
*/

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
};

struct Args {
    int n = 256;
    double rho = 0.8442;
    double temperature = 0.728;
    double dt = 0.005;
    int steps = 1000;
    double cutoff = 2.5;
    int sample_every = 10;
    int traj_every = 100;
    unsigned seed = 12345;
    std::string outdir = "results/cpp_lj";
};

static double sqr(double x) { return x * x; }

static double nearest_integer(double x) {
    return std::nearbyint(x);
}

void print_help() {
    std::cout << "Lennard-Jones MD demo in reduced units\n"
              << "Options:\n"
              << "  --n N                 number of particles [256]\n"
              << "  --rho R               reduced density [0.8442]\n"
              << "  --temperature T       initial reduced temperature [0.728]\n"
              << "  --dt DT               reduced time step [0.005]\n"
              << "  --steps S             number of MD steps [1000]\n"
              << "  --cutoff RC           cutoff radius [2.5]\n"
              << "  --sample-every K      write thermodynamics every K steps [10]\n"
              << "  --traj-every K        write XYZ every K steps [100]\n"
              << "  --seed SEED           random seed [12345]\n"
              << "  --outdir DIR          output directory [results/cpp_lj]\n";
}

Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value for " + name);
            }
            return argv[++i];
        };
        if (a == "--help" || a == "-h") {
            print_help();
            std::exit(0);
        } else if (a == "--n") args.n = std::stoi(need(a));
        else if (a == "--rho") args.rho = std::stod(need(a));
        else if (a == "--temperature") args.temperature = std::stod(need(a));
        else if (a == "--dt") args.dt = std::stod(need(a));
        else if (a == "--steps") args.steps = std::stoi(need(a));
        else if (a == "--cutoff") args.cutoff = std::stod(need(a));
        else if (a == "--sample-every") args.sample_every = std::stoi(need(a));
        else if (a == "--traj-every") args.traj_every = std::stoi(need(a));
        else if (a == "--seed") args.seed = static_cast<unsigned>(std::stoul(need(a)));
        else if (a == "--outdir") args.outdir = need(a);
        else throw std::runtime_error("Unknown argument: " + a);
    }
    if (args.n < 2) throw std::runtime_error("N must be at least 2.");
    if (args.rho <= 0.0 || args.temperature <= 0.0 || args.dt <= 0.0 || args.cutoff <= 0.0) {
        throw std::runtime_error("rho, temperature, dt and cutoff must be positive.");
    }
    return args;
}

std::vector<Vec3> initialize_positions(int n, double rho, double& box) {
    box = std::pow(static_cast<double>(n) / rho, 1.0 / 3.0);
    int cells = static_cast<int>(std::ceil(std::pow(static_cast<double>(n), 1.0 / 3.0)));
    double spacing = box / cells;
    std::vector<Vec3> pos;
    pos.reserve(n);
    for (int ix = 0; ix < cells; ++ix) {
        for (int iy = 0; iy < cells; ++iy) {
            for (int iz = 0; iz < cells; ++iz) {
                if (static_cast<int>(pos.size()) < n) {
                    pos.push_back({(ix + 0.5) * spacing, (iy + 0.5) * spacing, (iz + 0.5) * spacing});
                }
            }
        }
    }
    return pos;
}

std::vector<Vec3> initialize_velocities(int n, double temperature, unsigned seed) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> normal(0.0, std::sqrt(temperature));
    std::vector<Vec3> vel(n);
    Vec3 mean{0.0, 0.0, 0.0};
    for (int i = 0; i < n; ++i) {
        vel[i] = {normal(rng), normal(rng), normal(rng)};
        mean.x += vel[i].x; mean.y += vel[i].y; mean.z += vel[i].z;
    }
    mean.x /= n; mean.y /= n; mean.z /= n;
    double ekin = 0.0;
    for (int i = 0; i < n; ++i) {
        vel[i].x -= mean.x; vel[i].y -= mean.y; vel[i].z -= mean.z;
        ekin += 0.5 * (sqr(vel[i].x) + sqr(vel[i].y) + sqr(vel[i].z));
    }
    double current_temp = 2.0 * ekin / (3.0 * (n - 1));
    double scale = std::sqrt(temperature / current_temp);
    for (auto& v : vel) { v.x *= scale; v.y *= scale; v.z *= scale; }
    return vel;
}

void zero_forces(std::vector<Vec3>& f) {
    for (auto& v : f) { v.x = v.y = v.z = 0.0; }
}

void compute_forces(const std::vector<Vec3>& pos, std::vector<Vec3>& force, double box, double cutoff,
                    double& epot, double& virial) {
    const int n = static_cast<int>(pos.size());
    const double rc2 = cutoff * cutoff;
    const double inv_rc2 = 1.0 / rc2;
    const double inv_rc6 = inv_rc2 * inv_rc2 * inv_rc2;
    const double inv_rc12 = inv_rc6 * inv_rc6;
    const double u_shift = 4.0 * (inv_rc12 - inv_rc6);

    zero_forces(force);
    epot = 0.0;
    virial = 0.0;

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif
    std::vector<std::vector<Vec3>> local_forces(nthreads, std::vector<Vec3>(n));
    std::vector<double> local_epot(nthreads, 0.0);
    std::vector<double> local_virial(nthreads, 0.0);

#pragma omp parallel
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        auto& lf = local_forces[tid];
        double e = 0.0;
        double w = 0.0;

#pragma omp for schedule(static)
        for (int i = 0; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double dx = pos[i].x - pos[j].x;
                double dy = pos[i].y - pos[j].y;
                double dz = pos[i].z - pos[j].z;
                dx -= box * nearest_integer(dx / box);
                dy -= box * nearest_integer(dy / box);
                dz -= box * nearest_integer(dz / box);
                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 < rc2 && r2 > 0.0) {
                    const double inv_r2 = 1.0 / r2;
                    const double inv_r6 = inv_r2 * inv_r2 * inv_r2;
                    const double inv_r12 = inv_r6 * inv_r6;
                    e += 4.0 * (inv_r12 - inv_r6) - u_shift;
                    const double pref = 24.0 * inv_r2 * (2.0 * inv_r12 - inv_r6);
                    const double fx = pref * dx;
                    const double fy = pref * dy;
                    const double fz = pref * dz;
                    lf[i].x += fx; lf[i].y += fy; lf[i].z += fz;
                    lf[j].x -= fx; lf[j].y -= fy; lf[j].z -= fz;
                    w += dx * fx + dy * fy + dz * fz;
                }
            }
        }
        local_epot[tid] = e;
        local_virial[tid] = w;
    }

    for (int t = 0; t < nthreads; ++t) {
        epot += local_epot[t];
        virial += local_virial[t];
        for (int i = 0; i < n; ++i) {
            force[i].x += local_forces[t][i].x;
            force[i].y += local_forces[t][i].y;
            force[i].z += local_forces[t][i].z;
        }
    }
}

double kinetic_energy(const std::vector<Vec3>& vel) {
    double ekin = 0.0;
#pragma omp parallel for reduction(+ : ekin) schedule(static)
    for (int i = 0; i < static_cast<int>(vel.size()); ++i) {
        ekin += 0.5 * (sqr(vel[i].x) + sqr(vel[i].y) + sqr(vel[i].z));
    }
    return ekin;
}

void write_xyz(std::ofstream& out, const std::vector<Vec3>& pos, int step, double box) {
    out << pos.size() << "\n";
    out << "step=" << step << " box=" << std::fixed << std::setprecision(8) << box << "\n";
    for (const auto& p : pos) {
        out << "Ar " << std::fixed << std::setprecision(8)
            << p.x << " " << p.y << " " << p.z << "\n";
    }
}

int main(int argc, char** argv) {
    try {
        Args args = parse_args(argc, argv);
        std::filesystem::create_directories(args.outdir);

        double box = 0.0;
        auto pos = initialize_positions(args.n, args.rho, box);
        auto vel = initialize_velocities(args.n, args.temperature, args.seed);
        std::vector<Vec3> force(args.n);
        double epot = 0.0;
        double virial = 0.0;
        compute_forces(pos, force, box, args.cutoff, epot, virial);

        std::ofstream thermo(args.outdir + "/thermo.csv");
        std::ofstream xyz(args.outdir + "/trajectory.xyz");
        std::ofstream meta(args.outdir + "/metadata.txt");
        thermo << "step,time,ekin,epot,etot,temperature,pressure,wall_seconds\n";

        auto t0 = std::chrono::steady_clock::now();
        for (int step = 0; step <= args.steps; ++step) {
            double ekin = kinetic_energy(vel);
            double temp = 2.0 * ekin / (3.0 * (args.n - 1));
            double volume = box * box * box;
            double pressure = (args.n * temp + virial / 3.0) / volume;
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - t0).count();

            if (step % args.sample_every == 0) {
                thermo << step << "," << step * args.dt << ","
                       << std::setprecision(12) << ekin << "," << epot << "," << (ekin + epot) << ","
                       << temp << "," << pressure << "," << elapsed << "\n";
            }
            if (step % args.traj_every == 0) {
                write_xyz(xyz, pos, step, box);
            }
            if (step == args.steps) break;

#pragma omp parallel for schedule(static)
            for (int i = 0; i < args.n; ++i) {
                vel[i].x += 0.5 * args.dt * force[i].x;
                vel[i].y += 0.5 * args.dt * force[i].y;
                vel[i].z += 0.5 * args.dt * force[i].z;
                pos[i].x += args.dt * vel[i].x;
                pos[i].y += args.dt * vel[i].y;
                pos[i].z += args.dt * vel[i].z;
                pos[i].x -= box * std::floor(pos[i].x / box);
                pos[i].y -= box * std::floor(pos[i].y / box);
                pos[i].z -= box * std::floor(pos[i].z / box);
            }
            compute_forces(pos, force, box, args.cutoff, epot, virial);
#pragma omp parallel for schedule(static)
            for (int i = 0; i < args.n; ++i) {
                vel[i].x += 0.5 * args.dt * force[i].x;
                vel[i].y += 0.5 * args.dt * force[i].y;
                vel[i].z += 0.5 * args.dt * force[i].z;
            }
        }
        auto t1 = std::chrono::steady_clock::now();
        double total = std::chrono::duration<double>(t1 - t0).count();

        meta << "Lennard-Jones molecular dynamics demo\n";
        meta << "Implementation: C++17/OpenMP all-pairs\n";
        meta << "N=" << args.n << "\n";
        meta << "rho=" << args.rho << "\n";
        meta << "temperature=" << args.temperature << "\n";
        meta << "dt=" << args.dt << "\n";
        meta << "steps=" << args.steps << "\n";
        meta << "cutoff=" << args.cutoff << "\n";
        meta << "box=" << box << "\n";
        meta << "seed=" << args.seed << "\n";
#ifdef _OPENMP
        meta << "openmp_threads=" << omp_get_max_threads() << "\n";
#else
        meta << "openmp_threads=1\n";
#endif
        meta << "wall_seconds=" << total << "\n";

        std::cout << "Done. Results written to: " << args.outdir << "\n";
        std::cout << "Thermodynamics: " << args.outdir << "/thermo.csv\n";
        std::cout << "Trajectory:      " << args.outdir << "/trajectory.xyz\n";
        std::cout << "Wall time:       " << std::fixed << std::setprecision(3) << total << " s\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
