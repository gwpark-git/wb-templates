#include "bdmc/hydrodynamics.hpp"
#include "bdmc/reference.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double diameter = 1.;
constexpr double radius = diameter / 2.;
constexpr double r_min = 1.05, r_max = 2.50;
constexpr double epsilon = 24., force_cutoff = 1.65;

struct PairMobility {
    double relative_parallel, relative_transverse;
    double center_parallel, center_transverse;
};

double scaled_close(double value, double expected) {
    return std::abs(value - expected) / std::max(1., std::abs(expected));
}

// In sigma=kBT=D0=1 units, diffusivity equals the mobility response.
PairMobility rpy(double r, const bdmc::HydrodynamicSpec& spec) {
    const auto matrix = bdmc::hydrodynamic_mobility({0.,0.,0.,r,0.,0.}, spec);
    const double a = matrix[3], b = matrix[10];
    const double x = radius/r;
    if (scaled_close(a, 1.5*x-x*x*x) > 2e-14 ||
        scaled_close(b, .75*x+.5*x*x*x) > 2e-14)
        throw std::runtime_error("RPY pair mobility does not match its analytic free-space limit");
    return {2*(matrix[0]-a), 2*(matrix[7]-b),
            (matrix[0]+a)/2, (matrix[7]+b)/2};
}

PairMobility stokesian(const std::vector<double>& positions, const bdmc::StokesianSpec& spec) {
    const std::vector<double> zero(6,0.), strain(9,0.);
    std::array<double,2> relative{}, center{};
    for (int axis = 0; axis < 2; ++axis) {
        auto force = zero;
        force[axis] = -1.;
        force[axis+3] = 1.;
        const auto anti = bdmc::solve_stokesian(positions, force, zero, zero, zero, strain, spec);
        if (scaled_close(anti.translation[axis], -anti.translation[axis+3]) > 1e-12)
            throw std::runtime_error("antisymmetric force produced a COM drift");
        relative[axis] = anti.translation[axis+3] - anti.translation[axis];
        force[axis] = 1.;
        const auto common = bdmc::solve_stokesian(positions, force, zero, zero, zero, strain, spec);
        if (scaled_close(common.translation[axis], common.translation[axis+3]) > 1e-12)
            throw std::runtime_error("common force produced a relative drift");
        // The center responds to total force = 2.
        center[axis] = (common.translation[axis]+common.translation[axis+3])/4;
    }
    return {relative[0], relative[1], center[0], center[1]};
}

// The trace in a five-dimensional strain basis is the isotropic average
// of the stresslet operator. Each basis tensor has E:E=1.
double stresslet_average(const std::vector<double>& positions, const bdmc::StokesianSpec& spec) {
    const std::vector<double> zero(positions.size(), 0.);
    std::array<std::array<double,9>,5> strain{};
    strain[0][0] = 1/std::sqrt(2.);
    strain[0][4] = -strain[0][0];
    strain[1][0] = strain[1][4] = 1/std::sqrt(6.);
    strain[1][8] = -2*strain[1][0];
    strain[2][1] = strain[2][3] = strain[0][0];
    strain[3][2] = strain[3][6] = strain[0][0];
    strain[4][5] = strain[4][7] = strain[0][0];
    double projection = 0.;
    for (const auto& basis : strain) {
        const std::vector<double> deformation(basis.begin(), basis.end());
        std::vector<double> ambient(positions.size(), 0.);
        for (std::size_t i = 0; i < positions.size(); i += 3)
            for (std::size_t a = 0; a < 3; ++a)
                for (std::size_t b = 0; b < 3; ++b)
                    ambient[i+a] += basis[3*a+b]*positions[i+b];
        const auto response = bdmc::solve_stokesian(positions, zero, zero, ambient, zero, deformation, spec);
        for (std::size_t i = 0; i < response.stresslet.size(); ++i)
            projection += response.stresslet[i]*basis[i % basis.size()];
    }
    return projection/strain.size();
}

void write(PairMobility mobility) {
    for (double value : {mobility.relative_parallel, mobility.relative_transverse,
                         mobility.center_parallel, mobility.center_transverse}) {
        if (!std::isfinite(value) || value <= 0.)
            throw std::runtime_error("invalid pair mobility in the admitted shell");
        std::cout << ',' << value;
    }
}
} // namespace

int main(int argc, char* argv[]) try {
    if (argc != 2) throw std::invalid_argument("usage: bdmc_hi_equilibrium_grid CELLS");
    std::size_t used = 0;
    const std::string text = argv[1];
    const int cells = std::stoi(text, &used);
    if (used != text.size() || cells < 16 || cells > 2048)
        throw std::invalid_argument("CELLS must be an integer from 16 to 2048");
    const double pi = std::acos(-1.);
    bdmc::HydrodynamicSpec bd;
    bd.n_particles = 2;
    bd.radius = radius;
    bd.viscosity = 1/(3*pi);
    bdmc::StokesianSpec sd;
    sd.radius = radius;
    sd.viscosity = bd.viscosity;
    const double isolated = stresslet_average({0.,0.,0.}, sd);
    const auto far = stokesian({0.,0.,0.,1000.,0.,0.}, sd);
    if (scaled_close(isolated, 20*pi*sd.viscosity*std::pow(radius,3)/3) > 1e-12 ||
        scaled_close(far.relative_parallel, 2.) > 1e-3 ||
        scaled_close(far.center_parallel, .5) > 1e-3)
        throw std::runtime_error("Stokes drag/stresslet dilute reference failed");
    std::cout << "r_over_sigma,u_over_kBT,pair_virial_over_kBT,"
                 "rpy_relative_parallel,rpy_relative_transverse,rpy_center_parallel,rpy_center_transverse,"
                 "sd_relative_parallel,sd_relative_transverse,sd_center_parallel,sd_center_transverse,"
                 "sd_isotropic_stresslet_response\n"
              << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10-1);
    // Even indices are faces; odd indices are cell centers of the radial finite-volume solver.
    for (int i = 0; i <= 2*cells; ++i) {
        const double r = r_min+(r_max-r_min)*i/(2*cells);
        const std::vector<double> positions{0.,0.,0.,r,0.,0.};
        std::cout << r << ',' << bdmc::repulsive_pair_energy(r, epsilon, force_cutoff)
                  << ',' << r*bdmc::repulsive_radial_force(r, epsilon, force_cutoff);
        write(rpy(r, bd));
        write(stokesian(positions, sd));
        std::cout << ',' << stresslet_average(positions, sd) << '\n';
    }
    if (!std::cout) throw std::runtime_error("could not write mobility grid");
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
}
