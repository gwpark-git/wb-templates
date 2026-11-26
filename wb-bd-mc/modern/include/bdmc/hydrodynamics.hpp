#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace bdmc {

enum class HydrodynamicMobility { Oseen, RPY };
std::string hydrodynamic_mobility_name(HydrodynamicMobility);
HydrodynamicMobility parse_hydrodynamic_mobility(const std::string&);

// Dense, equal-sphere, free-space reference; all distances remain unwrapped.
// epsilon is a direct energy and sigma a force cutoff, independent of radius.
struct HydrodynamicSpec {
    std::size_t n_particles = 1;
    double radius = 1.0, viscosity = 1.0, dt = 0.001, thermal_energy = 1.0;
    HydrodynamicMobility mobility = HydrodynamicMobility::RPY;
    double epsilon = 0.0, sigma = 1.0, shear_rate = 0.0;
    std::uint64_t seed = 0;
    std::vector<double> external_forces, background_velocity;
    void validate() const;
};

struct HydrodynamicFrame {
    std::vector<double> positions, forces;
    double potential_energy = 0.0;
    std::array<double,9> pair_virial{};
    std::uint64_t step = 0;
    double time = 0.0;
};

// Assembly alone does not assert positive definiteness. The stepper checks it.
std::vector<double> hydrodynamic_mobility(const std::vector<double>& positions,
                                         const HydrodynamicSpec&);

// This wrapper owns a separate GSL stream and delegates dense mobility/noise to
// the retained native HI reference. It does not inject HI into Simulation's
// scalar-mobility integrator, and matching seeds do not align those RNG streams.
class HydrodynamicSimulation {
public:
    HydrodynamicSimulation(HydrodynamicSpec, std::vector<double> initial_positions);
    ~HydrodynamicSimulation();
    HydrodynamicSpec spec() const { return spec_; }
    void advance(std::uint64_t steps);
    HydrodynamicFrame observe() const;
    std::vector<HydrodynamicFrame> sample(std::uint64_t steps, std::uint64_t stride=1,
                                          bool include_initial=true);
private:
    struct RngState;
    HydrodynamicFrame evaluate(std::vector<double>, std::uint64_t) const;
    void advance_unlocked(std::uint64_t);
    const HydrodynamicSpec spec_;
    std::unique_ptr<RngState> rng_;
    HydrodynamicFrame state_;
    mutable std::mutex mutex_;
};

struct StokesianSpec {
    double radius = 1.0, viscosity = 1.0;
    bool lubrication = true;
    double lubrication_cutoff = 3.99;
    void validate() const;
};

struct StokesianResponse {
    std::vector<double> translation, rotation, stresslet;
    double relative_residual = 0.0, minimum_gap_over_radius = 0.0;
};

// Athermal instantaneous response, N<=64. Ambient strain is 9 or 9N entries,
// symmetric and traceless. Raw stresslets are not a complete suspension stress.
StokesianResponse solve_stokesian(const std::vector<double>& positions,
    const std::vector<double>& forces, const std::vector<double>& torques,
    const std::vector<double>& background_velocity,
    const std::vector<double>& background_angular_velocity,
    const std::vector<double>& strain, const StokesianSpec& = StokesianSpec{});

} // namespace bdmc
