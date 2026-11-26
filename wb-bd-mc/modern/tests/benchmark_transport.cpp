#include "bdmc/reference.hpp"
#include "bdmc/structure.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t frames = 2000;
constexpr std::size_t max_lag = 50;  // 0.4 diameter Brownian times.
constexpr double cadence = 0.008;

struct State {
    double phi_z;
    double eta;
};

State sample(bdmc::RepulsiveSpec spec, double phi, bool viscosity) {
    const double pi = std::acos(-1.);
    const double volume = spec.n_particles * pi / (6 * phi);
    spec.box_length = std::cbrt(volume);
    bdmc::RepulsiveSimulation simulation(spec);
    simulation.advance(2000);  // t=4 at dt=0.002.
    std::vector<double> modes;
    if (viscosity) modes.reserve(6 * frames);
    double trace = 0.;
    for (std::size_t i = 0; i < frames; ++i) {
        simulation.advance(4);
        const auto w = simulation.observe().pair_virial;
        trace += w[0] + w[4] + w[8];
        if (viscosity) {
            const double xy = 1. / (volume * std::sqrt(5.));
            const double normal = 1. / (volume * std::sqrt(30.));
            // Symmetric deviatoric stress autocorrelation divided by ten.
            modes.insert(modes.end(), {w[1]*xy, w[2]*xy, w[5]*xy, (w[0]-w[4])*normal,
                                       (w[0]-w[8])*normal, (w[4]-w[8])*normal});
        }
    }
    const double z = 1. + trace / (3. * spec.n_particles * frames);
    if (!viscosity) return {phi*z, 1.};
    const auto corr = bdmc::autocorrelation(modes, frames, 6, max_lag, true);
    double integral = 0.;
    for (std::size_t lag = 0; lag <= max_lag; ++lag) {
        const auto begin = corr.values.begin() + 6*lag;
        integral += (lag == 0 || lag == max_lag ? 0.5 : 1.) * std::accumulate(begin, begin + 6, 0.);
    }
    // Stokes–Einstein conversion from reduced configurational integral to eta/eta_s.
    return {phi*z, 1. + 3.*pi*volume*cadence*integral};
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--capabilities") {
            for (const auto& item : bdmc::execution_backends()) {
                std::cout << item.name << '\t' << item.precision << '\t' << item.compiled << '\t'
                          << item.available << '\t' << item.reason << '\n';
            }
            return 0;
        }
        if (argc != 5) throw std::invalid_argument("usage: bdmc_transport_benchmark serial_brute|serial_cell|metal N phi seed");
        const std::string label(argv[1]);
        if (label != "serial_brute" && label != "serial_cell" && label != "metal")
            throw std::invalid_argument("unknown requested backend");
        bdmc::RepulsiveSpec spec;
        spec.n_particles = std::stoull(argv[2]);
        const double phi = std::stod(argv[3]);
        if ((spec.n_particles != 27 && spec.n_particles != 64 && spec.n_particles != 216) ||
            !(phi > 0.006 && phi <= 0.494)) throw std::invalid_argument("unsupported benchmark state");
        spec.seed = std::stoull(argv[4]);
        spec.epsilon = 12.;
        spec.dt = 0.002;
        spec.execution_backend = label == "metal" ? "metal" : "serial";
        spec.force_search = label == "serial_cell" ? bdmc::SpatialBackend::CellList : bdmc::SpatialBackend::BruteForce;
        const auto begin = std::chrono::steady_clock::now();
        const double lo = phi - 0.006, hi = std::min(phi + 0.006, 0.494);
        const State left = sample(spec, lo, false);
        const State center = sample(spec, phi, true);
        const State right = phi == hi ? center : sample(spec, hi, false);
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        const double d = (right.phi_z - left.phi_z) / (hi - lo);
        if (!(std::isfinite(d) && d > 0. && std::isfinite(center.eta) && center.eta > 0.))
            throw std::runtime_error("transport estimator produced inadmissible values");
        std::cout << std::setprecision(17) << elapsed << ',' << center.phi_z/phi << ','
                  << d << ',' << center.eta << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
