#include "bdmc/core.hpp"

#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Function>
void rejects(Function function, const std::string& message) {
    bool rejected = false;
    try { function(); } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, message);
}

bdmc::SimulationSpec spec_1d(std::size_t particles) {
    bdmc::SimulationSpec spec;
    spec.n_particles = particles;
    spec.dimension = 1;
    spec.box_lengths = {10.0};
    spec.seed = 19041;
    return spec;
}

void equal_state(const bdmc::Simulation& a, const bdmc::Simulation& b) {
    const auto x = a.observe();
    const auto y = b.observe();
    check(x.state.ids == y.state.ids, "restored IDs differ");
    check(x.state.positions == y.state.positions, "restored wrapped coordinates differ");
    check(x.state.unwrapped_positions == y.state.unwrapped_positions, "restored unwrapped coordinates differ");
    check(x.state.step == y.state.step && x.state.time == y.state.time, "restored clock differs");
    check(x.force.values == y.force.values && x.force.potential_energy == y.force.potential_energy,
          "restored current-state force/energy differs");
    check(a.fixed_bonds() == b.fixed_bonds() && a.dynamic_bonds() == b.dynamic_bonds(), "restored topology differs");
    check(a.has_explicit_initial_positions() == b.has_explicit_initial_positions(), "restored initialization metadata differs");
    const auto rng_a = a.checkpoint().noise;
    const auto rng_b = b.checkpoint().noise;
    check(rng_a.engine_state == rng_b.engine_state && rng_a.has_spare == rng_b.has_spare &&
          rng_a.spare == rng_b.spare && rng_a.policy == rng_b.policy, "restored RNG/cache differs");
}

void brownian_uniform_exact_continuation() {
    auto spec = spec_1d(7);
    spec.diffusivity = 0.7;
    const std::vector<double> initial{-20.1, -0.01, 0.0, 0.2, 2.0, 9.99, 15.7};
    const std::vector<bdmc::ParticleId> ids{100, 3, 90, 71, 4, 81, 11};
    bdmc::Simulation continuous(spec, initial, ids), split(spec, initial, ids);
    continuous.advance(93);
    split.advance(17);
    const auto snapshot = split.checkpoint();
    check(snapshot.schema == "bdmc.checkpoint" && snapshot.version == 1, "checkpoint schema identity");
    auto restored = bdmc::Simulation::from_checkpoint(snapshot);
    equal_state(split, *restored);
    split.advance(76);
    restored->advance(76);
    equal_state(continuous, split);  // checkpoint did not consume noise
    equal_state(continuous, *restored);
}

void gaussian_odd_and_even_cache_continuation() {
    for (const std::uint64_t prefix : {0, 1, 2, 7}) {
        auto spec = spec_1d(3);
        spec.noise_policy = bdmc::NoisePolicy::Gaussian;
        bdmc::Simulation continuous(spec), split(spec);
        continuous.advance(60);
        split.advance(prefix);
        const auto snapshot = split.checkpoint();
        check(snapshot.noise.has_spare == (prefix % 2 != 0), "odd Gaussian draw cache captured");
        check(!snapshot.explicit_initial_positions, "default initialization provenance captured");
        auto restored = bdmc::Simulation::from_checkpoint(snapshot);
        restored->advance(60 - prefix);
        equal_state(continuous, *restored);
    }
}

void dumbbell_both_policies_and_roundoff() {
    for (const auto policy : {bdmc::NoisePolicy::LegacyUniform, bdmc::NoisePolicy::Gaussian}) {
        auto spec = spec_1d(4);
        spec.model = bdmc::Model::Dumbbell;
        spec.spring_alpha = 1.2;
        spec.dt = 0.004;
        spec.mobility = 0.8;
        spec.diffusivity = 0.6;
        spec.noise_policy = policy;
        const std::vector<double> initial{9.9, -0.3, 0.2, 0.8};
        const std::vector<bdmc::ParticleId> ids{30, 10, 80, 7};
        bdmc::Simulation continuous(spec, initial, ids), split(spec, initial, ids);
        continuous.advance(1700);
        split.advance(1601);
        auto restored = bdmc::Simulation::from_checkpoint(split.checkpoint());
        restored->advance(99);
        equal_state(continuous, *restored);
    }
}

void dynamic_storage_roundtrip_and_later_edits() {
    auto spec = spec_1d(4);
    spec.model = bdmc::Model::Dumbbell;
    bdmc::Simulation original(spec, {0.1, 0.2, 0.4, 0.7}, {40, 10, 90, 5});
    original.apply_topology_edits({{bdmc::BondAction::Add, {90, 10}},
                                    {bdmc::BondAction::Add, {40, 5}}});
    original.advance(11);
    auto restored = bdmc::Simulation::from_checkpoint(original.checkpoint());
    equal_state(original, *restored);
    const std::vector<bdmc::BondEdit> edits{
        {bdmc::BondAction::Remove, {10, 90}}, {bdmc::BondAction::Add, {90, 5}}};
    original.apply_topology_edits(edits);
    restored->apply_topology_edits(edits);
    original.advance(19);
    restored->advance(19);
    equal_state(original, *restored);
}

void unsupported_components_rejected_explicitly() {
    auto spec = spec_1d(2);
    bdmc::Simulation custom_force(spec);
    custom_force.add_force_term(bdmc::ConstantForce{{0.0}});
    rejects([&] { custom_force.checkpoint(); }, "even a zero custom force is not silently serialized");
    bdmc::Simulation updater(spec);
    updater.set_topology_updater([](const bdmc::TopologyUpdateContext&) {
        return std::vector<bdmc::BondEdit>{{bdmc::BondAction::Add, {0, 1}}};
    });
    rejects([&] { updater.checkpoint(); }, "active topology callback rejected");
    updater.advance(1);
    updater.set_topology_updater({});
    auto restored = bdmc::Simulation::from_checkpoint(updater.checkpoint());
    equal_state(updater, *restored);  // static edited graph is serializable
}

void corrupted_checkpoint_rejected() {
    auto spec = spec_1d(4);
    spec.model = bdmc::Model::Dumbbell;
    bdmc::Simulation original(spec, {0.1, 0.2, 0.4, 0.8}, {10, 20, 30, 40});
    original.advance(7);
    const auto good = original.checkpoint();
    using Mutator = std::function<void(bdmc::SimulationCheckpoint&)>;
    const std::vector<std::pair<std::string, Mutator>> corruptions{
        {"schema", [](auto& c) { c.schema = "other"; }},
        {"version", [](auto& c) { c.version = 2; }},
        {"core version", [](auto& c) { c.core_version = "future"; }},
        {"backend", [](auto& c) { c.backend = "gpu"; }},
        {"boundary", [](auto& c) { c.boundary = "shear"; }},
        {"geometry/spec box mismatch", [](auto& c) { c.box_lengths[0] *= 2; }},
        {"particle count", [](auto& c) { ++c.spec.n_particles; }},
        {"state dimension", [](auto& c) { c.state.dimension = 2; }},
        {"duplicate IDs", [](auto& c) { c.state.ids[1] = c.state.ids[0]; }},
        {"position shape", [](auto& c) { c.state.positions.pop_back(); }},
        {"unwrapped shape", [](auto& c) { c.state.unwrapped_positions.pop_back(); }},
        {"nonfinite coordinates", [](auto& c) { c.state.positions[0] = std::numeric_limits<double>::quiet_NaN(); }},
        {"nonfinite unwrapped coordinates", [](auto& c) { c.state.unwrapped_positions[0] = std::numeric_limits<double>::infinity(); }},
        {"wrapped bounds", [](auto& c) { c.state.positions[0] = c.box_lengths[0]; }},
        {"inconsistent unwrapped coordinate", [](auto& c) { c.state.unwrapped_positions[0] += 0.1; }},
        {"clock", [](auto& c) { c.state.time += c.spec.dt; }},
        {"step", [](auto& c) { ++c.state.step; }},
        {"different fixed pairing", [](auto& c) { c.fixed_bonds = {{10, 20}, {30, 40}}; }},
        {"missing fixed bond", [](auto& c) { c.fixed_bonds.pop_back(); }},
        {"unknown dynamic endpoint", [](auto& c) { c.dynamic_bonds = {{10, 999}}; }},
        {"self dynamic bond", [](auto& c) { c.dynamic_bonds = {{10, 10}}; }},
        {"duplicate dynamic bond", [](auto& c) { c.dynamic_bonds = {{10, 20}, {10, 20}}; }},
        {"noncanonical dynamic bond", [](auto& c) { c.dynamic_bonds = {{20, 10}}; }},
        {"RNG engine", [](auto& c) { c.noise.engine = "unknown"; }},
        {"RNG/spec policy", [](auto& c) { c.noise.policy = bdmc::NoisePolicy::Gaussian; }},
        {"uniform cache", [](auto& c) { c.noise.has_spare = true; }},
        {"nonfinite spare", [](auto& c) { c.noise.spare = std::numeric_limits<double>::quiet_NaN(); }},
        {"malformed RNG", [](auto& c) { c.noise.engine_state = "not RNG state"; }},
        {"missing RNG words", [](auto& c) { c.noise.engine_state = "1 2 3"; }},
        {"extra RNG words", [](auto& c) { c.noise.engine_state += " 0"; }},
        {"signed RNG word", [](auto& c) { c.noise.engine_state = "-1" + c.noise.engine_state.substr(c.noise.engine_state.find(' ')); }},
        {"out-of-range RNG cursor", [](auto& c) { c.noise.engine_state.replace(c.noise.engine_state.rfind(' ') + 1, std::string::npos, "313"); }},
        {"zero RNG state", [](auto& c) {
            c.noise.engine_state.clear();
            for (std::size_t i = 0; i <= std::mt19937_64::state_size; ++i) c.noise.engine_state += "0 ";
        }},
    };
    for (const auto& corruption : corruptions) {
        auto malformed = good;
        corruption.second(malformed);
        rejects([&] { bdmc::Simulation::from_checkpoint(malformed); }, "corrupt checkpoint accepted: " + corruption.first);
    }
    auto restored = bdmc::Simulation::from_checkpoint(good);
    equal_state(original, *restored);
    check(original.state().step == 7, "failed restore attempts never mutate source simulation");
    std::cout << "  rejected " << corruptions.size() << " structural/state/RNG corruptions\n";
}

void rng_snapshot_cache_roundtrip() {
    bdmc::NoiseSource original(1701, bdmc::NoisePolicy::Gaussian);
    (void)original.standardized();
    const auto snapshot = original.snapshot();
    check(snapshot.has_spare, "standalone Gaussian spare captured");
    auto restored = bdmc::NoiseSource::from_snapshot(snapshot);
    for (std::size_t i = 0; i < 101; ++i) {
        check(original.standardized() == restored.standardized(), "RNG/cache continuation must be exact");
    }
    check(original.snapshot().engine_state == restored.snapshot().engine_state, "RNG state stays synchronized");
}

}  // namespace

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"brownian_uniform_exact_continuation", brownian_uniform_exact_continuation},
        {"gaussian_odd_and_even_cache_continuation", gaussian_odd_and_even_cache_continuation},
        {"dumbbell_both_policies_and_roundoff", dumbbell_both_policies_and_roundoff},
        {"dynamic_storage_roundtrip_and_later_edits", dynamic_storage_roundtrip_and_later_edits},
        {"unsupported_components_rejected_explicitly", unsupported_components_rejected_explicitly},
        {"corrupted_checkpoint_rejected", corrupted_checkpoint_rejected},
        {"rng_snapshot_cache_roundtrip", rng_snapshot_cache_roundtrip},
    };
    int failures = 0;
    for (const auto& test : tests) {
        try {
            test.second();
            std::cout << "PASS " << test.first << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cerr << "FAIL " << test.first << ": " << error.what() << '\n';
        }
    }
    std::cout << tests.size() - failures << '/' << tests.size() << " checkpoint check groups passed\n";
    return failures == 0 ? 0 : 1;
}
