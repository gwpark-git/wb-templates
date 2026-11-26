#include "bdmc/core.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <limits>
#include <locale>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace bdmc {
namespace {

// All callers use literals. Construct the exception message only on failure.
void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}

bool finite_nonnegative(double value) {
    return std::isfinite(value) && value >= 0.0;
}

Bond canonical(Bond bond) {
    if (bond.second < bond.first) std::swap(bond.first, bond.second);
    return bond;
}

bool bond_less(const Bond& a, const Bond& b) {
    return std::make_pair(a.first, a.second) < std::make_pair(b.first, b.second);
}

void validate_ids(const std::vector<ParticleId>& ids) {
    std::unordered_set<ParticleId> seen;
    for (const auto id : ids) require(seen.insert(id).second, "duplicate particle ID");
}

void validate_bond(const Bond& bond, const std::unordered_set<ParticleId>& ids) {
    require(bond.first != bond.second, "self bond is not allowed");
    require(ids.count(bond.first) && ids.count(bond.second), "bond uses unknown particle ID");
}

SimulationSpec validated(SimulationSpec spec) {
    spec.validate();
    return spec;
}

ParticleState initial_state(const SimulationSpec& spec,
                            std::vector<double> positions,
                            std::vector<ParticleId> ids,
                            const OrthorhombicBox& box) {
    ParticleState state;
    state.dimension = spec.dimension;
    if (ids.empty()) {
        ids.resize(spec.n_particles);
        std::iota(ids.begin(), ids.end(), ParticleId{0});
    }
    require(ids.size() == spec.n_particles, "particle ID count differs from specification");
    state.ids = std::move(ids);
    if (positions.empty()) positions.assign(spec.n_particles * spec.dimension, 0.0);
    require(positions.size() == spec.n_particles * spec.dimension,
            "initial positions must have n_particles * dimension values");
    state.unwrapped_positions = positions;
    for (std::size_t k = 0; k < positions.size(); ++k) {
        positions[k] = box.wrap(positions[k], k % spec.dimension);
    }
    state.positions = std::move(positions);
    state.validate();
    return state;
}

}  // namespace

std::string model_name(Model model) {
    switch (model) {
        case Model::Brownian: return "brownian";
        case Model::Dumbbell: return "dumbbell";
    }
    throw std::invalid_argument("unknown model enum");
}

std::string noise_name(NoisePolicy policy) {
    switch (policy) {
        case NoisePolicy::LegacyUniform: return "legacy_uniform";
        case NoisePolicy::Gaussian: return "gaussian";
    }
    throw std::invalid_argument("unknown noise policy enum");
}

Model parse_model(const std::string& name) {
    if (name == "brownian") return Model::Brownian;
    if (name == "dumbbell") return Model::Dumbbell;
    throw std::invalid_argument("model must be 'brownian' or 'dumbbell'; transient-network MC is not implemented in the modern reference core");
}

NoisePolicy parse_noise_policy(const std::string& name) {
    if (name == "legacy_uniform") return NoisePolicy::LegacyUniform;
    if (name == "gaussian") return NoisePolicy::Gaussian;
    throw std::invalid_argument("noise_policy must be 'legacy_uniform' or 'gaussian'");
}

void SimulationSpec::validate() const {
    (void)model_name(model);
    (void)noise_name(noise_policy);
    require(n_particles > 0, "n_particles must be positive");
    require(dimension >= 1 && dimension <= 3, "dimension must be 1, 2 or 3");
    require(n_particles <= std::numeric_limits<std::size_t>::max() / dimension,
            "particle buffer size overflows size_t");
    require(box_lengths.size() == dimension, "box_lengths must match dimension");
    for (const auto length : box_lengths) {
        require(std::isfinite(length) && length > 0.0, "box lengths must be finite and positive");
    }
    require(std::isfinite(dt) && dt > 0.0, "dt must be finite and positive");
    require(finite_nonnegative(mobility), "mobility must be finite and nonnegative");
    require(finite_nonnegative(diffusivity), "diffusivity must be finite and nonnegative");
    require(std::isfinite(spring_alpha) && spring_alpha > 0.0,
            "spring_alpha must be finite and positive");
    require(std::isfinite(mobility * dt), "mobility * dt overflows");
    require(std::isfinite(2.0 * (diffusivity * dt)), "noise variance overflows");
    require(std::isfinite(static_cast<double>(dimension) * spring_alpha * spring_alpha),
            "spring coefficient overflows");
    if (model == Model::Dumbbell) {
        require(n_particles % 2 == 0, "dumbbell model requires an even particle count");
        // q_next=(1-2*mu*k*dt)q+noise. A stochastic invariant second
        // moment exists only for |1-2*mu*k*dt|<1 when mu is nonzero.
        const double coefficient = static_cast<double>(dimension) * spring_alpha * spring_alpha;
        require(mobility == 0.0 || mobility * coefficient * dt < 1.0,
                "dumbbell Euler timestep is outside the linear stability interval mu*k*dt < 1");
    }
}

OrthorhombicBox::OrthorhombicBox(std::vector<double> lengths)
    : lengths_(std::move(lengths)) {
    require(!lengths_.empty() && lengths_.size() <= 3, "box dimension must be 1, 2 or 3");
    for (const auto length : lengths_) {
        require(std::isfinite(length) && length > 0.0, "invalid box length");
    }
}

double OrthorhombicBox::wrap(double value, std::size_t coordinate) const {
    require(std::isfinite(value), "position must be finite");
    const double length = lengths_.at(coordinate);
    // For an in-box value fmod(v, L) returns v. Keep zero canonicalized.
    if (value >= 0.0 && value < length) return value == 0.0 ? 0.0 : value;
    double wrapped = std::fmod(value, length);
    if (wrapped < 0.0) wrapped += length;
    // Rounding an infinitesimal negative input can give exactly L.
    if (wrapped >= length) wrapped = 0.0;
    return wrapped == 0.0 ? 0.0 : wrapped;
}

double OrthorhombicBox::minimum_image(double displacement, std::size_t coordinate) const {
    require(std::isfinite(displacement), "displacement must be finite");
    const double length = lengths_.at(coordinate);
    double result = std::fmod(displacement, length);
    if (result >= 0.5 * length) result -= length;
    if (result < -0.5 * length) result += length;
    return result;
}

std::size_t ParticleState::index_of(ParticleId id) const {
    const auto found = std::find(ids.begin(), ids.end(), id);
    if (found == ids.end()) throw std::invalid_argument("unknown particle ID");
    return static_cast<std::size_t>(found - ids.begin());
}

void ParticleState::validate() const {
    require(dimension >= 1 && dimension <= 3, "invalid state dimension");
    require(!ids.empty(), "particle state is empty");
    validate_ids(ids);
    require(positions.size() == ids.size() * dimension, "invalid position buffer size");
    require(unwrapped_positions.size() == positions.size(), "invalid unwrapped buffer size");
    for (const auto value : positions) require(std::isfinite(value), "nonfinite position");
    for (const auto value : unwrapped_positions) require(std::isfinite(value), "nonfinite unwrapped position");
    require(finite_nonnegative(time), "invalid state time");
}

FixedTopology::FixedTopology(std::vector<Bond> bonds, const std::vector<ParticleId>& ids,
                             bool require_disjoint_pairs) {
    validate_ids(ids);
    const std::unordered_set<ParticleId> known(ids.begin(), ids.end());
    std::unordered_set<ParticleId> used;
    for (auto& bond : bonds) {
        bond = canonical(bond);
        validate_bond(bond, known);
        if (require_disjoint_pairs) {
            require(used.insert(bond.first).second && used.insert(bond.second).second,
                    "fixed dumbbell pairs must be disjoint");
        }
    }
    std::sort(bonds.begin(), bonds.end(), bond_less);
    require(std::adjacent_find(bonds.begin(), bonds.end()) == bonds.end(), "duplicate fixed bond");
    if (require_disjoint_pairs) require(used.size() == ids.size(), "every particle must have one dumbbell partner");
    bonds_ = std::move(bonds);
}

FixedTopology FixedTopology::dumbbells(const std::vector<ParticleId>& ids) {
    require(!ids.empty() && ids.size() % 2 == 0, "dumbbells require an even nonzero particle count");
    std::vector<Bond> pairs;
    pairs.reserve(ids.size() / 2);
    for (std::size_t i = 0; i < ids.size() / 2; ++i) pairs.push_back({ids[i], ids[i + ids.size() / 2]});
    return FixedTopology(std::move(pairs), ids, true);
}

DynamicTopology::DynamicTopology(std::vector<ParticleId> ids) : ids_(std::move(ids)) {
    validate_ids(ids_);
}

void DynamicTopology::apply(const std::vector<BondEdit>& edits) {
    // Validate on a candidate so a late invalid edit cannot leave half a graph.
    // This protects graph consistency only; proposal probabilities and chain-end
    // kinetics belong to a separate scientific controller.
    auto candidate = bonds_;
    const std::unordered_set<ParticleId> known(ids_.begin(), ids_.end());
    for (const auto& edit : edits) {
        const auto bond = canonical(edit.bond);
        validate_bond(bond, known);
        const auto found = std::lower_bound(candidate.begin(), candidate.end(), bond, bond_less);
        const bool present = found != candidate.end() && *found == bond;
        if (edit.action == BondAction::Add) {
            require(!present, "cannot add an existing dynamic bond");
            candidate.insert(found, bond);
        } else if (edit.action == BondAction::Remove) {
            require(present, "cannot remove a missing dynamic bond");
            candidate.erase(found);
        } else {
            throw std::invalid_argument("unknown bond action");
        }
    }
    bonds_.swap(candidate);
}

void HarmonicBondForce::operator()(const ParticleState& state, const OrthorhombicBox& box,
                                   const FixedTopology& fixed, const DynamicTopology& dynamic,
                                   ForceResult& result) const {
    require(finite_nonnegative(coefficient), "harmonic coefficient must be finite and nonnegative");
    require(result.values.size() == state.positions.size(), "force buffer shape mismatch");
    // Build one O(N) ID lookup, avoiding O(N^2) endpoint searches and retaining
    // stable IDs independently of storage order.
    std::unordered_map<ParticleId, std::size_t> storage;
    storage.reserve(state.size());
    for (std::size_t i = 0; i < state.size(); ++i) storage.emplace(state.ids[i], i);
    const auto& bonds = use_dynamic_bonds ? dynamic.bonds() : fixed.bonds();
    for (const auto& bond : bonds) {
        const auto i = storage.at(bond.first);
        const auto j = storage.at(bond.second);
        for (std::size_t c = 0; c < state.dimension; ++c) {
            const auto ic = i * state.dimension + c;
            const auto jc = j * state.dimension + c;
            const double q = box.minimum_image(state.positions[jc] - state.positions[ic], c);
            const double force = coefficient * q;
            result.values[ic] += force;
            result.values[jc] -= force;
            result.potential_energy += 0.5 * coefficient * q * q;
        }
    }
}

void ConstantForce::operator()(const ParticleState& state, const OrthorhombicBox&,
                               const FixedTopology&, const DynamicTopology&,
                               ForceResult& result) const {
    require(per_coordinate.size() == state.dimension, "constant force dimension mismatch");
    require(result.values.size() == state.positions.size(), "force buffer shape mismatch");
    for (const auto force : per_coordinate) require(std::isfinite(force), "constant force must be finite");
    for (std::size_t k = 0; k < result.values.size(); ++k) {
        result.values[k] += per_coordinate[k % state.dimension];
    }
    // A constant periodic drive has no globally single-valued periodic energy.
    // It contributes force only, so energy remains the conservative-term sum.
}

NoiseSource::NoiseSource(std::uint64_t seed, NoisePolicy policy)
    : engine_(seed), policy_(policy) {
    (void)noise_name(policy);
}

double NoiseSource::uniform() {
    // Explicit 53-bit conversion avoids implementation-specific uniform_real_distribution.
    return static_cast<double>(engine_() >> 11U) * 0x1.0p-53;
}

double NoiseSource::standardized() {
    // Variance matching preserves the chosen Brownian second moment, not the
    // entire finite-step law. Keep policy selection visible in run metadata.
    if (policy_ == NoisePolicy::LegacyUniform) return std::sqrt(3.0) * (2.0 * uniform() - 1.0);
    if (has_spare_) {
        has_spare_ = false;
        return spare_;
    }
    const double radius = std::sqrt(-2.0 * std::log(1.0 - uniform()));
    const double angle = 6.283185307179586476925286766559 * uniform();
    spare_ = radius * std::sin(angle);
    has_spare_ = true;
    return radius * std::cos(angle);
}

NoiseSource::Snapshot NoiseSource::snapshot() const {
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << engine_;
    require(static_cast<bool>(output), "could not serialize mt19937_64 state");
    return {"std::mt19937_64", output.str(), policy_, has_spare_, spare_};
}

NoiseSource NoiseSource::from_snapshot(const Snapshot& snapshot) {
    require(snapshot.engine == "std::mt19937_64", "unsupported checkpoint RNG engine");
    (void)noise_name(snapshot.policy);
    require(std::isfinite(snapshot.spare), "checkpoint Gaussian spare must be finite");
    if (snapshot.policy == NoisePolicy::LegacyUniform) {
        require(!snapshot.has_spare && snapshot.spare == 0.0,
                "uniform-noise checkpoint must not contain a Gaussian cache");
    }
    require(!snapshot.engine_state.empty() && snapshot.engine_state.size() <= 20000,
            "invalid checkpoint RNG text length");

    // Standard mt19937_64 stream format: 312 unsigned state words followed by
    // its cursor. Validate grammar and cursor independently before relying on
    // the standard-library extraction operator; no executable serialization.
    std::istringstream words(snapshot.engine_state);
    words.imbue(std::locale::classic());
    constexpr auto state_words = std::mt19937_64::state_size;
    std::size_t count = 0;
    bool any_nonzero = false;
    std::string token;
    while (words >> token) {
        require(count < state_words + 1, "checkpoint RNG has extra state words");
        std::uint64_t value = 0;
        const auto parsed = std::from_chars(token.data(), token.data() + token.size(), value);
        require(parsed.ec == std::errc{} && parsed.ptr == token.data() + token.size(),
                "checkpoint RNG state words must be unsigned decimal integers");
        if (count < state_words) any_nonzero = any_nonzero || value != 0;
        else require(value <= state_words, "checkpoint RNG cursor is out of range");
        ++count;
    }
    require(count == state_words + 1, "checkpoint RNG has missing state words");
    require(any_nonzero, "checkpoint RNG has a degenerate all-zero state");
    NoiseSource restored(0, snapshot.policy);
    std::istringstream input(snapshot.engine_state);
    input.imbue(std::locale::classic());
    input >> restored.engine_;
    require(!input.fail(), "could not parse checkpoint RNG state");
    restored.has_spare_ = snapshot.has_spare;
    restored.spare_ = snapshot.spare;
    return restored;
}

void OverdampedEuler::step(ParticleState& state, const OrthorhombicBox& box,
                           const ForceResult& forces, ConstantMobility mobility,
                           double diffusivity, double dt, NoiseSource& noise) const {
    require(forces.values.size() == state.positions.size(), "integrator force shape mismatch");
    require(std::isfinite(dt) && dt > 0.0 && finite_nonnegative(diffusivity) &&
            finite_nonnegative(mobility.value), "invalid integrator transport parameters");
    require(state.step < std::numeric_limits<std::uint64_t>::max(), "step counter overflow");
    const double scale = std::sqrt(2.0 * (diffusivity * dt));
    // Every drift uses the already-complete old-state force. Updating a bead
    // and then recomputing its partner's force would change the Euler scheme.
    // Separate candidate coordinate buffers also avoid exposing half a frame.
    std::vector<double> wrapped(state.positions.size());
    std::vector<double> unwrapped(state.positions.size());
    for (std::size_t k = 0; k < state.positions.size(); ++k) {
        const double random = diffusivity == 0.0 ? 0.0 : scale * noise.standardized();
        const double displacement = (mobility.value * dt) * forces.values[k] + random;
        unwrapped[k] = state.unwrapped_positions[k] + displacement;
        require(std::isfinite(unwrapped[k]), "integration produced nonfinite coordinates");
        wrapped[k] = box.wrap(state.positions[k] + displacement, k % state.dimension);
    }
    const auto next_step = state.step + 1;
    const double next_time = static_cast<double>(next_step) * dt;
    require(std::isfinite(next_time), "simulation time overflow");
    state.positions.swap(wrapped);
    state.unwrapped_positions.swap(unwrapped);
    state.step = next_step;
    state.time = next_time;
}

Simulation::Simulation(SimulationSpec spec, std::vector<double> initial_positions,
                         std::vector<ParticleId> particle_ids)
    : spec_(validated(std::move(spec))), explicit_initial_positions_(!initial_positions.empty()),
      box_(spec_.box_lengths),
      state_(initial_state(spec_, std::move(initial_positions), std::move(particle_ids), box_)),
      fixed_(spec_.model == Model::Dumbbell ? FixedTopology::dumbbells(state_.ids) : FixedTopology{}),
      dynamic_(state_.ids), noise_(spec_.seed, spec_.noise_policy) {
    if (spec_.model == Model::Dumbbell) {
        forces_.emplace_back(HarmonicBondForce{
            static_cast<double>(spec_.dimension) * spec_.spring_alpha * spec_.spring_alpha, false});
    }
}

void SimulationCheckpoint::validate() const {
    require(schema == "bdmc.checkpoint", "unsupported checkpoint schema");
    require(version == 1, "unsupported checkpoint schema version");
    require(core_version == "0.1.0", "unsupported checkpoint core version");
    require(backend == "modern_cpu_reference", "unsupported checkpoint backend");
    spec.validate();
    require(boundary == "orthorhombic_periodic", "unsupported checkpoint boundary");
    require(box_lengths == spec.box_lengths, "checkpoint geometry/spec box mismatch");
    state.validate();
    require(state.dimension == spec.dimension && state.size() == spec.n_particles,
            "checkpoint particle state/spec shape mismatch");
    const double expected_time = static_cast<double>(state.step) * spec.dt;
    require(std::isfinite(expected_time) && state.time == expected_time,
            "checkpoint time must equal absolute_step*dt");
    const OrthorhombicBox box(box_lengths);
    for (std::size_t k = 0; k < state.positions.size(); ++k) {
        const auto c = k % spec.dimension;
        const double length = box_lengths[c];
        require(state.positions[k] >= 0.0 && state.positions[k] < length,
                "checkpoint positions must be wrapped into [0,L)");
        const double residual = box.minimum_image(
            box.wrap(state.unwrapped_positions[k], c) - state.positions[k], c);
        // Wrapped and unwrapped coordinates are accumulated separately in the
        // existing integrator. Permit its small rounding drift, but never
        // normalize or otherwise change either checkpoint array on restore.
        const long double scale = std::max(static_cast<long double>(length),
            std::abs(static_cast<long double>(state.unwrapped_positions[k])));
        const long double budget = 128.0L * std::numeric_limits<double>::epsilon() *
            (static_cast<long double>(state.step) + 1.0L) * scale;
        const long double tolerance = std::min(1e-8L * static_cast<long double>(length),
            std::max(128.0L * std::numeric_limits<double>::epsilon() * length, budget));
        require(std::abs(static_cast<long double>(residual)) <= tolerance,
                "checkpoint wrapped/unwrapped coordinates disagree beyond rounding tolerance");
    }

    const auto expected_fixed = spec.model == Model::Dumbbell
        ? FixedTopology::dumbbells(state.ids).bonds() : std::vector<Bond>{};
    require(fixed_bonds == expected_fixed, "checkpoint fixed topology does not match default model pairing");
    DynamicTopology validated_dynamic(state.ids);
    std::vector<BondEdit> edits;
    edits.reserve(dynamic_bonds.size());
    for (const auto& bond : dynamic_bonds) edits.push_back({BondAction::Add, bond});
    validated_dynamic.apply(edits);
    require(validated_dynamic.bonds() == dynamic_bonds,
            "checkpoint dynamic topology must be canonical and sorted");
    require(noise.policy == spec.noise_policy, "checkpoint RNG/spec noise policy mismatch");
    (void)NoiseSource::from_snapshot(noise);
}

Simulation::Simulation(const SimulationCheckpoint& snapshot, RestoreTag)
    : spec_(snapshot.spec), explicit_initial_positions_(snapshot.explicit_initial_positions),
      box_(snapshot.box_lengths), state_(snapshot.state),
      fixed_(snapshot.fixed_bonds, state_.ids, spec_.model == Model::Dumbbell),
      dynamic_(state_.ids), noise_(NoiseSource::from_snapshot(snapshot.noise)) {
    if (spec_.model == Model::Dumbbell) {
        forces_.emplace_back(HarmonicBondForce{
            static_cast<double>(spec_.dimension) * spec_.spring_alpha * spec_.spring_alpha, false});
    }
    std::vector<BondEdit> edits;
    edits.reserve(snapshot.dynamic_bonds.size());
    for (const auto& bond : snapshot.dynamic_bonds) edits.push_back({BondAction::Add, bond});
    dynamic_.apply(edits);
}

SimulationCheckpoint Simulation::checkpoint() const {
    std::lock_guard<std::mutex> lock(mutex_);
    // Serializing positions while silently dropping a callback/extra force
    // would restore a different model. Refuse that unsupported state explicitly.
    require(!has_custom_force_terms_ && !topology_updater_,
            "checkpoint supports only default model forces and no active topology callback");
    SimulationCheckpoint snapshot;
    snapshot.spec = spec_;
    snapshot.box_lengths = box_.lengths();
    snapshot.state = state_;
    snapshot.fixed_bonds = fixed_.bonds();
    snapshot.dynamic_bonds = dynamic_.bonds();
    snapshot.noise = noise_.snapshot();
    snapshot.explicit_initial_positions = explicit_initial_positions_;
    snapshot.validate();
    return snapshot;
}

std::unique_ptr<Simulation> Simulation::from_checkpoint(const SimulationCheckpoint& snapshot) {
    snapshot.validate();
    return std::unique_ptr<Simulation>(new Simulation(snapshot, RestoreTag{}));
}

SimulationSpec Simulation::spec() const { return spec_; }

ParticleState Simulation::state() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return state_;
}

std::vector<Bond> Simulation::fixed_bonds() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return fixed_.bonds();
}

std::vector<Bond> Simulation::dynamic_bonds() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return dynamic_.bonds();
}

ForceResult Simulation::force_unlocked() const {
    // One shared additive buffer permits independent force laws. A new term
    // must preserve shape/finite values and declare whether it has an energy;
    // nonconservative drives cannot invent a periodic potential.
    ForceResult result{std::vector<double>(state_.positions.size(), 0.0), 0.0};
    for (const auto& term : forces_) term(state_, box_, fixed_, dynamic_, result);
    require(result.values.size() == state_.positions.size(), "force term changed buffer shape");
    for (const auto value : result.values) require(std::isfinite(value), "force term produced nonfinite force");
    require(std::isfinite(result.potential_energy), "force term produced nonfinite energy");
    return result;
}

void Simulation::validate_horizon(std::uint64_t steps) const {
    require(steps <= std::numeric_limits<std::uint64_t>::max() - state_.step, "step counter overflow");
    require(std::isfinite(static_cast<double>(state_.step + steps) * spec_.dt), "simulation time overflow");
}

void Simulation::advance_unlocked(std::uint64_t steps) {
    validate_horizon(steps);
    // Ordering is physical: topology edits at t_old, then all forces, then a
    // simultaneous position step. Observation stride must not move this clock
    // or consume extra random values. Custom callbacks own their acceptance law.
    for (std::uint64_t i = 0; i < steps; ++i) {
        if (topology_updater_) {
            dynamic_.apply(topology_updater_({state_, box_, fixed_, dynamic_, spec_.dt}));
        }
        const auto force = force_unlocked();
        integrator_.step(state_, box_, force, ConstantMobility{spec_.mobility},
                         spec_.diffusivity, spec_.dt, noise_);
    }
}

void Simulation::advance(std::uint64_t steps) {
    std::lock_guard<std::mutex> lock(mutex_);
    advance_unlocked(steps);
}

ObservationFrame Simulation::observe_unlocked() const {
    return {state_, force_unlocked(), dynamic_.bonds()};
}

ObservationFrame Simulation::observe() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return observe_unlocked();
}

ObservationBatch Simulation::sample(std::uint64_t steps, std::uint64_t stride,
                                     bool include_initial) {
    std::lock_guard<std::mutex> lock(mutex_);
    require(stride > 0, "sample stride must be positive");
    validate_horizon(steps);
    ObservationBatch batch;
    batch.fixed_bonds = fixed_.bonds();
    const auto chunks = steps / stride + (steps % stride != 0 ? 1 : 0);
    require(chunks < std::numeric_limits<std::size_t>::max(), "observation count overflows size_t");
    const auto count = static_cast<std::size_t>(chunks) + (include_initial ? 1 : 0);
    batch.frames.reserve(count);
    if (include_initial) batch.frames.push_back(observe_unlocked());
    std::uint64_t remaining = steps;
    while (remaining > 0) {
        const auto chunk = std::min(stride, remaining);
        advance_unlocked(chunk);
        batch.frames.push_back(observe_unlocked());
        remaining -= chunk;
    }
    return batch;
}

void Simulation::add_force_term(ForceTerm term) {
    std::lock_guard<std::mutex> lock(mutex_);
    require(static_cast<bool>(term), "force term cannot be empty");
    forces_.push_back(std::move(term));
    has_custom_force_terms_ = true;
}

void Simulation::set_topology_updater(TopologyUpdater updater) {
    std::lock_guard<std::mutex> lock(mutex_);
    topology_updater_ = std::move(updater);
}

void Simulation::apply_topology_edits(const std::vector<BondEdit>& edits) {
    std::lock_guard<std::mutex> lock(mutex_);
    dynamic_.apply(edits);
}

}  // namespace bdmc
