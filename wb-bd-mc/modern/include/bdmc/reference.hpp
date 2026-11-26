#pragma once

#include "bdmc/core.hpp"
#include "bdmc/execution.hpp"
#include "bdmc/spatial.hpp"

#include <array>
#include <mutex>
#include <vector>

namespace bdmc {

// These are additive reference models. They do not change Simulation, its
// historical conventions, RNG policy defaults, or restart format.
struct RepulsiveSpec {
    std::size_t n_particles = 64;
    double box_length = 10.0;
    double dt = 0.001;
    double mobility = 1.0;
    double diffusion = 1.0;
    double epsilon = 1.0;
    double sigma = 1.0;
    std::uint64_t seed = 0;
    NoisePolicy noise_policy = NoisePolicy::Gaussian;
    SpatialBackend force_search = SpatialBackend::CellList;
    std::string execution_backend = "serial";
    int backend_threads = 1;
    void validate() const;
};

struct RepulsiveFrame {
    std::vector<double> positions, unwrapped_positions, forces;
    double potential_energy = 0.0;
    // Unnormalised positive pair virial sum r_ij (outer) F_on_j, counted once.
    std::array<double, 9> pair_virial{};
    std::uint64_t step = 0;
    double time = 0.0;
};

double repulsive_pair_energy(double r, double epsilon, double sigma);
double repulsive_radial_force(double r, double epsilon, double sigma);

class RepulsiveSimulation {
public:
    explicit RepulsiveSimulation(RepulsiveSpec spec, std::vector<double> initial_positions = {});
    RepulsiveSpec spec() const { return spec_; }
    bool explicit_initial_state() const { return explicit_initial_state_; }
    void advance(std::uint64_t steps);
    RepulsiveFrame observe() const;
    std::vector<RepulsiveFrame> sample(std::uint64_t steps, std::uint64_t stride = 1,
                                     bool include_initial = true);
private:
    RepulsiveFrame evaluate(std::vector<double> positions, std::vector<double> unwrapped,
                            std::uint64_t step) const;
    void advance_unlocked(std::uint64_t steps);
    const RepulsiveSpec spec_;
    const bool explicit_initial_state_;
    std::unique_ptr<ExecutionKernels> kernels_;
    OrthorhombicBox box_;
    SpatialSearch search_;
    NoiseSource noise_;
    RepulsiveFrame state_;
    mutable std::mutex mutex_;
};

enum class SpringLaw { Hookean, FENE };
enum class FlowProtocol { NoFlow, SimpleShear, StepShear, OscillatoryShear };
std::string spring_law_name(SpringLaw);
SpringLaw parse_spring_law(const std::string&);
std::string flow_protocol_name(FlowProtocol);
FlowProtocol parse_flow_protocol(const std::string&);

// This engine evolves q=r_tail-r_head directly in unbounded connector space.
// Each bead contributes mobility and diffusion, giving connector drift -2*mu*F
// and noise variance 4*D*dt. It is distinct from periodic bead-pair Simulation.
struct DumbbellSpec {
    std::size_t n_dumbbells = 64;
    double dt = 0.001;
    double mobility = 1.0;   // mobility of each bead, not connector mobility
    double diffusion = 1.0;  // diffusivity of each bead
    double spring_constant = 1.0;
    SpringLaw spring_law = SpringLaw::Hookean;
    double q0 = 10.0;
    FlowProtocol protocol = FlowProtocol::NoFlow;
    double shear_rate = 0.0;
    double step_strain = 0.0;
    double strain_amplitude = 0.0;
    double angular_frequency = 1.0;
    std::uint64_t seed = 0;
    NoisePolicy noise_policy = NoisePolicy::Gaussian;
    std::string execution_backend = "serial";
    int backend_threads = 1;
    void validate() const;
};

struct DumbbellFrame {
    std::vector<double> connectors, spring_forces;
    std::array<double, 9> conformation{}, polymer_stress{};
    double mean_energy = 0.0;
    double shear_rate = 0.0;
    std::uint64_t step = 0;
    double time = 0.0;
};

double dumbbell_spring_energy(const std::array<double, 3>& q, const DumbbellSpec& spec);
std::array<double, 3> dumbbell_spring_force(const std::array<double, 3>& q,
                                         const DumbbellSpec& spec);
double protocol_shear_rate(const DumbbellSpec& spec, double time);

class DumbbellSimulation {
public:
    explicit DumbbellSimulation(DumbbellSpec spec, std::vector<double> initial_connectors = {});
    DumbbellSpec spec() const { return spec_; }
    bool explicit_initial_state() const { return explicit_initial_state_; }
    void advance(std::uint64_t steps);
    DumbbellFrame observe() const;
    std::vector<DumbbellFrame> sample(std::uint64_t steps, std::uint64_t stride = 1,
                                    bool include_initial = true);
private:
    DumbbellFrame evaluate(std::vector<double> connectors, std::uint64_t step) const;
    void advance_unlocked(std::uint64_t steps);
    const DumbbellSpec spec_;
    const bool explicit_initial_state_;
    std::unique_ptr<ExecutionKernels> kernels_;
    NoiseSource noise_;
    DumbbellFrame state_;
    mutable std::mutex mutex_;
};

} // namespace bdmc
