#include "bdmc/hydrodynamics.hpp"
#include "bdmc/reference.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double phi = .10;
constexpr double gap_floor = 1.05;
constexpr double epsilon = 24., cutoff = 1.65;
constexpr double move_width = .28;
constexpr int discard_sweeps = 2000, block_sweeps = 160, pressure_stride = 4, blocks = 6;
constexpr std::array<double,2> volume_steps{.002, .004};

double distance(const std::vector<double>& q, std::size_t i, std::size_t j) {
    const double x = q[3*i]-q[3*j];
    const double y = q[3*i+1]-q[3*j+1];
    const double z = q[3*i+2]-q[3*j+2];
    return std::sqrt(x*x+y*y+z*z);
}

class EquilibriumCube {
public:
    EquilibriumCube(std::size_t count, double length, std::uint64_t seed)
        : coordinates_(3*count), count_(count), length_(length), rng_(seed) {
        int side = std::ceil(std::cbrt(double(count)));
        while (side*side*side < static_cast<int>(count)) ++side;
        std::vector<int> sites(side*side*side);
        for (int i = 0; i < static_cast<int>(sites.size()); ++i) sites[i] = i;
        std::shuffle(sites.begin(), sites.end(), rng_);
        for (std::size_t i = 0; i < count_; ++i) {
            const int k = sites[i];
            const std::array<int,3> grid{k%side, k/side%side, k/(side*side)};
            for (std::size_t d = 0; d < grid.size(); ++d)
                coordinates_[3*i+d] = length_*(.5+grid[d])/side;
        }
    }

    void sweep() {
        std::uniform_int_distribution<std::size_t> particle(0, count_-1);
        for (std::size_t trial = 0; trial < count_; ++trial) {
            const auto i = particle(rng_);
            const double before = local_energy(i);
            std::array<double,3> old{};
            bool inside = true;
            for (std::size_t d = 0; d < old.size(); ++d) {
                old[d] = coordinates_[3*i+d];
                coordinates_[3*i+d] += move_width*(2*uniform_(rng_)-1);
                inside = inside && coordinates_[3*i+d] > 0. && coordinates_[3*i+d] < length_;
            }
            ++attempts_;
            // Reflecting geometric boundaries enter the equilibrium move as rejections.
            if (inside && uniform_(rng_) < std::exp(before-local_energy(i))) {
                ++accepted_;
            } else {
                std::copy(old.begin(), old.end(), coordinates_.begin()+3*i);
            }
        }
    }

    const std::vector<double>& coordinates() const { return coordinates_; }

    double acceptance() const { return static_cast<double>(accepted_)/attempts_; }

    std::array<double,2> compression_weights() const {
        std::array<double,2> energy{};
        static const std::array<double,2> scale{
            std::exp(-volume_steps[0]/3), std::exp(-volume_steps[1]/3)};
        for (std::size_t i = 0; i < count_; ++i)
            for (std::size_t j = i+1; j < count_; ++j) {
                const double r = distance(coordinates_, i, j);
                if (r*scale.back() >= cutoff) continue;
                const double original = bdmc::repulsive_pair_energy(r, epsilon, cutoff);
                for (std::size_t e = 0; e < energy.size(); ++e)
                    energy[e] += r*scale[e] < gap_floor ? std::numeric_limits<double>::infinity()
                        : bdmc::repulsive_pair_energy(r*scale[e], epsilon, cutoff)-original;
            }
        return {std::exp(-energy[0]), std::exp(-energy[1])};
    }

private:
    double local_energy(std::size_t bead) const {
        double energy = 0.;
        for (std::size_t j = 0; j < count_; ++j) {
            if (j == bead) continue;
            const double r = distance(coordinates_, bead, j);
            if (r < gap_floor) return std::numeric_limits<double>::infinity();
            if (r < cutoff) energy += bdmc::repulsive_pair_energy(r, epsilon, cutoff);
        }
        return energy;
    }

    std::vector<double> coordinates_;
    const std::size_t count_;
    const double length_;
    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> uniform_{0.,1.};
    std::size_t accepted_ = 0, attempts_ = 0;
};

std::array<double,3> responses(const std::vector<double>& q, double volume, double pi) {
    const std::size_t n = q.size()/3, dim = q.size();
    bdmc::HydrodynamicSpec rpy;
    rpy.n_particles = n;
    rpy.radius = .5;
    rpy.viscosity = 1/(3*pi);
    const auto mobility = bdmc::hydrodynamic_mobility(q, rpy);
    double trace = 0.;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t d = 0; d < 3; ++d)
                trace += mobility[(3*i+d)*dim+3*j+d];
    bdmc::StokesianSpec sd;
    sd.radius = .5;
    sd.viscosity = rpy.viscosity;
    const std::vector<double> zero(dim,0.), no_strain(9,0.);
    double common = 0.;
    for (std::size_t d = 0; d < 3; ++d) {
        auto force = zero;
        for (std::size_t i = d; i < dim; i += 3) force[i] = 1.;
        const auto fluid = bdmc::solve_stokesian(q, force, zero, zero, zero, no_strain, sd);
        for (std::size_t i = d; i < dim; i += 3) common += fluid.translation[i];
    }
    // An affine xy strain and its imposed ambient velocities, with the
    // stresslet contraction reported as a resolved high-frequency term only.
    auto strain = no_strain, ambient = zero;
    strain[1] = strain[3] = 1/std::sqrt(2.);
    for (std::size_t i = 0; i < dim; i += 3) {
        ambient[i] = strain[1]*q[i+1];
        ambient[i+1] = strain[3]*q[i];
    }
    const auto shear = bdmc::solve_stokesian(q, zero, zero, ambient, zero, strain, sd);
    double stresslet = 0.;
    for (std::size_t i = 0; i < shear.stresslet.size(); ++i)
        stresslet += shear.stresslet[i]*strain[i%9];
    return {trace/(3*n), common/(3*n), stresslet/(2*sd.viscosity*volume)};
}
} // namespace

int main(int argc, char* argv[]) try {
    if (argc != 3) throw std::invalid_argument("usage: bdmc_hi_size PARTICLES SEED");
    const std::string text = argv[1], seed_text = argv[2];
    std::size_t used = 0;
    const unsigned long particles = std::stoul(text, &used);
    if (used != text.size() || !(particles == 8 || particles == 16 || particles == 32 || particles == 64))
        throw std::invalid_argument("supported sizes are 8, 16, 32 and 64");
    const auto seed = std::stoull(seed_text, &used);
    if (used != seed_text.size() || seed_text.empty() || seed_text.front() == '-')
        throw std::invalid_argument("seed must be a nonnegative integer");
    const double pi = std::acos(-1.);
    const double volume = particles*pi/(6*phi), length = std::cbrt(volume);
    EquilibriumCube system(particles, length, seed);
    for (int i = 0; i < discard_sweeps; ++i) system.sweep();
    std::cout << "N,seed,block,L_over_sigma,phi,pure_K1,rpy_open_H0,sd_open_H0,sd_stresslet_eta_increment,"
                 "mean_weight_eps_002,mean_weight_eps_004,acceptance\n"
              << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10-1);
    for (int block = 0; block < blocks; ++block) {
        std::array<double,2> total{};
        for (int i = 0; i < block_sweeps; ++i) {
            system.sweep();
            if (i%pressure_stride == 0) {
                const auto weight = system.compression_weights();
                for (std::size_t j = 0; j < total.size(); ++j) total[j] += weight[j];
            }
        }
        const auto hydro = responses(system.coordinates(), volume, pi);
        for (const double value : hydro)
            if (!std::isfinite(value) || value <= 0.) throw std::runtime_error("invalid hydrodynamic response");
        std::cout << particles << ',' << seed << ',' << block << ',' << length << ',' << phi << ',' << 1.
                  << ',' << hydro[0] << ',' << hydro[1] << ',' << hydro[2]
                  << ',' << total[0]/(block_sweeps/pressure_stride)
                  << ',' << total[1]/(block_sweeps/pressure_stride) << ',' << system.acceptance() << std::endl;
    }
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
}
