#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace bdmc {

using ParticleId = std::uint64_t;

enum class Model { Brownian, Dumbbell };
enum class NoisePolicy { LegacyUniform, Gaussian };

std::string model_name(Model model);
std::string noise_name(NoisePolicy policy);
Model parse_model(const std::string& name);
NoisePolicy parse_noise_policy(const std::string& name);

// Configuration is a value copied into Simulation: later edits to a caller's
// spec cannot change a running trajectory. In these reduced units, mobility*F*dt
// is drift and 2*diffusivity*dt is coordinate variance; changing either transport
// coefficient changes the model, not just its execution speed.
struct SimulationSpec {
    Model model = Model::Brownian;
    std::size_t n_particles = 64;
    std::size_t dimension = 3;
    std::vector<double> box_lengths = {10.0, 10.0, 10.0};
    double dt = 0.001;
    double mobility = 1.0;
    double diffusivity = 1.0;
    double spring_alpha = 1.0;
    std::uint64_t seed = 0;
    NoisePolicy noise_policy = NoisePolicy::LegacyUniform;

    void validate() const;
};

class OrthorhombicBox {
public:
    explicit OrthorhombicBox(std::vector<double> lengths);
    std::size_t dimension() const noexcept { return lengths_.size(); }
    const std::vector<double>& lengths() const noexcept { return lengths_; }
    double wrap(double value, std::size_t coordinate) const;
    // The tie convention is [-L/2, L/2); topology never depends on cell lists.
    double minimum_image(double displacement, std::size_t coordinate) const;

private:
    std::vector<double> lengths_;
};

// IDs identify particles across topology, observations and restart. Flat
// buffers use particle-major (N,dimension) order; wrapped coordinates serve
// periodic forces while unwrapped coordinates retain displacement history.
// Reordering storage must move IDs and both coordinate buffers together.
struct ParticleState {
    std::size_t dimension = 0;
    std::vector<ParticleId> ids;
    std::vector<double> positions;
    std::vector<double> unwrapped_positions;
    std::uint64_t step = 0;
    double time = 0.0;

    std::size_t size() const noexcept { return ids.size(); }
    std::size_t index_of(ParticleId id) const;
    void validate() const;
};

struct Bond {
    ParticleId first = 0;
    ParticleId second = 0;
    friend bool operator==(const Bond& a, const Bond& b) noexcept {
        return a.first == b.first && a.second == b.second;
    }
};

class FixedTopology {
public:
    FixedTopology() = default;
    FixedTopology(std::vector<Bond> bonds, const std::vector<ParticleId>& ids,
                  bool require_disjoint_pairs = false);
    static FixedTopology dumbbells(const std::vector<ParticleId>& ids);
    const std::vector<Bond>& bonds() const noexcept { return bonds_; }

private:
    std::vector<Bond> bonds_;
};

enum class BondAction { Add, Remove };
struct BondEdit {
    BondAction action = BondAction::Add;
    Bond bond;
};

// A simple-edge graph seam for the reference core: no loops or parallel edges.
// Exact physical chain/end identities live separately in NetworkTopology
// (network.hpp). Adding edges here does not implement transient-network MC.
class DynamicTopology {
public:
    explicit DynamicTopology(std::vector<ParticleId> ids);
    const std::vector<Bond>& bonds() const noexcept { return bonds_; }
    // Apply the complete proposed batch only if every edit is valid. This is a
    // graph-storage contract, not a chain-end model or an MC acceptance law.
    void apply(const std::vector<BondEdit>& edits);

private:
    std::vector<ParticleId> ids_;
    std::vector<Bond> bonds_;
};

struct ForceResult {
    std::vector<double> values;
    double potential_energy = 0.0;
};

// Terms add into a zero-initialized force/energy buffer. They read an immutable
// state and topology; they never own particle storage or perform integration.
using ForceTerm = std::function<void(const ParticleState&, const OrthorhombicBox&,
                                     const FixedTopology&, const DynamicTopology&,
                                     ForceResult&)>;
struct HarmonicBondForce {
    double coefficient = 1.0;
    bool use_dynamic_bonds = false;
    void operator()(const ParticleState&, const OrthorhombicBox&,
                    const FixedTopology&, const DynamicTopology&, ForceResult&) const;
};
struct ConstantForce {
    std::vector<double> per_coordinate;
    void operator()(const ParticleState&, const OrthorhombicBox&,
                    const FixedTopology&, const DynamicTopology&, ForceResult&) const;
};

// The Gaussian spare is part of the stochastic state, just like the engine
// cursor. Saving only the seed or engine would change an odd-draw continuation.
// Both policies have unit variance, but their finite-step fourth moments differ.
class NoiseSource {
public:
    struct Snapshot {
        std::string engine = "std::mt19937_64";
        std::string engine_state;
        NoisePolicy policy = NoisePolicy::LegacyUniform;
        bool has_spare = false;
        double spare = 0.0;
    };

    NoiseSource(std::uint64_t seed, NoisePolicy policy);
    double standardized();
    NoisePolicy policy() const noexcept { return policy_; }
    Snapshot snapshot() const;
    static NoiseSource from_snapshot(const Snapshot& snapshot);

private:
    double uniform();
    std::mt19937_64 engine_;
    NoisePolicy policy_;
    bool has_spare_ = false;
    double spare_ = 0.0;
};

struct ConstantMobility {
    double value = 1.0;
};

struct OverdampedEuler {
    // Forces have already been evaluated on the old state for all particles.
    void step(ParticleState&, const OrthorhombicBox&, const ForceResult&,
              ConstantMobility, double diffusivity, double dt, NoiseSource&) const;
};

struct TopologyUpdateContext {
    const ParticleState& state;
    const OrthorhombicBox& box;
    const FixedTopology& fixed;
    const DynamicTopology& dynamic;
    double dt;
};
// Called once before force evaluation at old-state time. The updater owns its
// RNG/acceptance law and returns edits; an empty callback means no updates.
using TopologyUpdater =
    std::function<std::vector<BondEdit>(const TopologyUpdateContext&)>;

struct ObservationFrame {
    ParticleState state;
    ForceResult force;
    std::vector<Bond> dynamic_bonds;
};
struct ObservationBatch {
    std::vector<ObservationFrame> frames;
    std::vector<Bond> fixed_bonds;
};

// Restart is a coupled snapshot of geometry, absolute clock, particle identity,
// graph and RNG. Observation arrays alone omit information needed for replay.
// This schema deliberately covers default reference components only.
struct SimulationCheckpoint {
    std::string schema = "bdmc.checkpoint";
    std::uint32_t version = 1;
    std::string core_version = "0.1.0";
    std::string backend = "modern_cpu_reference";
    SimulationSpec spec;
    std::string boundary = "orthorhombic_periodic";
    std::vector<double> box_lengths;
    ParticleState state;
    std::vector<Bond> fixed_bonds;
    std::vector<Bond> dynamic_bonds;
    NoiseSource::Snapshot noise;
    bool explicit_initial_positions = false;

    // Format/physical-state validation, not authentication or proof of history.
    void validate() const;
};

// This small owning orchestrator composes force, topology, noise and stepping
// responsibilities. Native force/update callbacks are the available runtime
// seams; replacing the concrete integrator/mobility needs explicit wiring and
// new physical checks. It is not a completed general network/HI controller.
class Simulation {
public:
    explicit Simulation(SimulationSpec spec,
                        std::vector<double> initial_positions = {},
                        std::vector<ParticleId> particle_ids = {});
    SimulationSpec spec() const;
    bool has_explicit_initial_positions() const noexcept { return explicit_initial_positions_; }
    ParticleState state() const;
    std::vector<Bond> fixed_bonds() const;
    std::vector<Bond> dynamic_bonds() const;
    void advance(std::uint64_t steps);
    ObservationFrame observe() const;
    // Sample at stride and always at the final requested step, with an optional
    // initial frame. Advancing zero steps produces only the optional initial.
    ObservationBatch sample(std::uint64_t steps, std::uint64_t stride,
                            bool include_initial = true);
    SimulationCheckpoint checkpoint() const;
    static std::unique_ptr<Simulation> from_checkpoint(const SimulationCheckpoint& snapshot);

    // Native-only component seams. Register before running; callbacks run
    // inside the simulation lock and must not call back into this simulation.
    void add_force_term(ForceTerm term);
    void set_topology_updater(TopologyUpdater updater);
    void apply_topology_edits(const std::vector<BondEdit>& edits);

private:
    struct RestoreTag {};
    Simulation(const SimulationCheckpoint& snapshot, RestoreTag);
    ForceResult force_unlocked() const;
    ObservationFrame observe_unlocked() const;
    void advance_unlocked(std::uint64_t steps);
    void validate_horizon(std::uint64_t steps) const;

    const SimulationSpec spec_;
    const bool explicit_initial_positions_;
    OrthorhombicBox box_;
    ParticleState state_;
    FixedTopology fixed_;
    DynamicTopology dynamic_;
    std::vector<ForceTerm> forces_;
    bool has_custom_force_terms_ = false;
    TopologyUpdater topology_updater_;
    NoiseSource noise_;
    OverdampedEuler integrator_;
    mutable std::mutex mutex_;
};

}  // namespace bdmc
