#include "bdmc/core.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void close(double actual, double expected, double tolerance, const std::string& message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message + ": actual=" + std::to_string(actual) +
                                 ", expected=" + std::to_string(expected));
    }
}

template <typename Function>
void rejects(Function function, const std::string& message) {
    bool rejected = false;
    try { function(); } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, message);
}

bdmc::SimulationSpec spec_1d(std::size_t particles = 2) {
    bdmc::SimulationSpec spec;
    spec.n_particles = particles;
    spec.dimension = 1;
    spec.box_lengths = {1000.0};
    return spec;
}

void geometry() {
    bdmc::OrthorhombicBox box({10.0, 6.0});
    close(box.wrap(-0.2, 0), 9.8, 1e-14, "negative periodic wrapping");
    close(box.wrap(30.2, 0), 0.2, 1e-14, "multiple positive crossings");
    close(box.wrap(-20.0, 0), 0.0, 0.0, "exact boundary wrapping");
    close(box.minimum_image(9.6, 0), -0.4, 1e-14, "minimum image");
    close(box.minimum_image(5.0, 0), -5.0, 0.0, "positive half-box tie");
    close(box.minimum_image(-5.0, 0), -5.0, 0.0, "negative half-box tie");
    close(box.minimum_image(16.0, 1), -2.0, 0.0, "orthorhombic coordinate");
    rejects([&] { box.wrap(std::numeric_limits<double>::infinity(), 0); }, "reject nonfinite geometry");
    rejects([] { bdmc::OrthorhombicBox invalid({0.0}); }, "reject zero box");
}

void specifications() {
    auto spec = spec_1d();
    spec.validate();
    spec.dt = 0.0;
    rejects([&] { spec.validate(); }, "reject nonpositive dt");
    spec = spec_1d();
    spec.diffusivity = -1.0;
    rejects([&] { spec.validate(); }, "reject negative diffusivity");
    spec = spec_1d();
    spec.model = bdmc::Model::Dumbbell;
    spec.n_particles = 3;
    rejects([&] { spec.validate(); }, "reject odd dumbbell population");
    spec.n_particles = 2;
    spec.dt = 1.0;
    rejects([&] { spec.validate(); }, "reject unstable Euler spring step");
    spec = spec_1d();
    spec.dimension = 4;
    rejects([&] { spec.validate(); }, "reject unsupported dimensions");
    spec = spec_1d();
    spec.box_lengths = {1.0, 2.0};
    rejects([&] { spec.validate(); }, "reject box mismatch");
    spec = spec_1d();
    spec.dimension = 3;
    spec.box_lengths = {1.0, 1.0, 1.0};
    spec.n_particles = std::numeric_limits<std::size_t>::max();
    rejects([&] { spec.validate(); }, "reject buffer size overflow");
    rejects([] { bdmc::parse_model("transient_network"); }, "no invented network model");
    rejects([] { bdmc::parse_noise_policy("default"); }, "noise policy is explicit");
}

void identity_and_topology() {
    const std::vector<bdmc::ParticleId> ids{40, 10, 90, 5};
    const auto fixed = bdmc::FixedTopology::dumbbells(ids);
    check(fixed.bonds() == std::vector<bdmc::Bond>{{5, 10}, {40, 90}}, "legacy pairing uses initial storage order with stable IDs");
    rejects([&] { bdmc::FixedTopology bad({{40, 90}, {90, 40}}, ids); }, "duplicate bond rejected");
    rejects([&] { bdmc::FixedTopology bad({{40, 40}}, ids); }, "self bond rejected");
    rejects([&] { bdmc::FixedTopology bad({{40, 99}}, ids); }, "unknown endpoint rejected");
    rejects([&] { bdmc::FixedTopology bad({{40, 90}, {40, 10}}, ids, true); }, "non-disjoint pairs rejected");
    bdmc::DynamicTopology dynamic(ids);
    dynamic.apply({{bdmc::BondAction::Add, {90, 40}}});
    check(dynamic.bonds() == std::vector<bdmc::Bond>{{40, 90}}, "dynamic canonical storage");
    rejects([&] {
        dynamic.apply({{bdmc::BondAction::Add, {10, 5}}, {bdmc::BondAction::Remove, {40, 10}}});
    }, "invalid edit batch rejected");
    check(dynamic.bonds() == std::vector<bdmc::Bond>{{40, 90}}, "failed edit batch is atomic");
    dynamic.apply({{bdmc::BondAction::Remove, {90, 40}}});
    check(dynamic.bonds().empty(), "remove canonical dynamic edge");
    auto spec = spec_1d(4);
    rejects([&] { bdmc::Simulation bad(spec, {}, {1, 1, 2, 3}); }, "duplicate particle IDs rejected");
    rejects([&] { bdmc::Simulation bad(spec, {0.0}); }, "initial position shape rejected");
}

void harmonic_force_and_old_state_update() {
    auto spec = spec_1d();
    spec.model = bdmc::Model::Dumbbell;
    spec.box_lengths = {10.0};
    spec.spring_alpha = 2.0;
    spec.diffusivity = 0.0;
    spec.mobility = 0.75;
    spec.dt = 0.01;
    bdmc::Simulation simulation(spec, {9.8, 0.2}, {100, 7});
    const auto initial = simulation.observe();
    close(initial.force.values[0], 1.6, 1e-13, "minimum-image harmonic force");
    close(initial.force.values[1], -1.6, 1e-13, "opposite endpoint force");
    close(initial.force.potential_energy, 0.32, 1e-13, "harmonic pair energy counted once");
    const double eps = 1e-5;
    bdmc::Simulation plus(spec, {9.8 + eps, 0.2});
    bdmc::Simulation minus(spec, {9.8 - eps, 0.2});
    const double negative_derivative = -(plus.observe().force.potential_energy - minus.observe().force.potential_energy) / (2 * eps);
    close(initial.force.values[0], negative_derivative, 1e-8, "force matches independent energy finite difference");
    simulation.advance(1);
    const auto after = simulation.observe();
    close(after.state.positions[0], 9.812, 1e-13, "first bead Euler update");
    close(after.state.positions[1], 0.188, 1e-13, "second bead reads old-state force");
    close(after.state.unwrapped_positions[0] + after.state.unwrapped_positions[1], 10.0, 1e-13,
          "deterministic pair center conserved");
    close(after.force.values[0], 1.504, 1e-13, "observation force belongs to updated state");
}

void wrapping_and_component_composition() {
    auto spec = spec_1d(1);
    spec.box_lengths = {1.0};
    spec.dt = 0.25;
    spec.diffusivity = 0.0;
    bdmc::Simulation simulation(spec, {0.9}, {99});
    simulation.add_force_term(bdmc::ConstantForce{{3.0}});
    simulation.add_force_term(bdmc::ConstantForce{{1.0}});
    simulation.advance(2);
    const auto after = simulation.observe();
    close(after.state.positions[0], 0.9, 1e-14, "wrapped displacement with multiple crossings");
    close(after.state.unwrapped_positions[0], 2.9, 1e-14, "unwrapped displacement retained");
    close(after.force.values[0], 4.0, 0.0, "force terms add");
    check(after.state.ids[0] == 99, "stable ID survives stepping");
    close(after.state.time, 0.5, 0.0, "absolute step time convention");

    auto pair_spec = spec_1d();
    pair_spec.diffusivity = 0.0;
    pair_spec.dt = 0.1;
    bdmc::Simulation dynamic(pair_spec, {0.0, 1.0});
    dynamic.add_force_term(bdmc::HarmonicBondForce{2.0, true});
    dynamic.set_topology_updater([](const bdmc::TopologyUpdateContext& context) {
        check(context.dt == 0.1 && context.state.step == 0, "topology update time/phase");
        return std::vector<bdmc::BondEdit>{{bdmc::BondAction::Add, {0, 1}}};
    });
    dynamic.advance(1);
    close(dynamic.state().positions[0], 0.2, 1e-14, "topology edits precede force evaluation");
    check(dynamic.dynamic_bonds().size() == 1, "dynamic storage separate from fixed storage");
    check(dynamic.fixed_bonds().empty(), "dynamic edges do not alter permanent topology");
    dynamic.set_topology_updater({});
    dynamic.apply_topology_edits({{bdmc::BondAction::Remove, {1, 0}}});
    dynamic.advance(1);
    close(dynamic.state().positions[0], 0.2, 1e-14, "removed bond stops contributing force");
}

void seeded_stream_and_observation_schedule() {
    auto spec = spec_1d(32);
    spec.seed = 1234;
    bdmc::Simulation first(spec), split(spec), batched(spec);
    first.advance(17);
    split.advance(5);
    split.advance(12);
    check(first.state().unwrapped_positions == split.state().unwrapped_positions, "batch boundaries preserve RNG stream");
    auto batch = batched.sample(17, 6, true);
    check(batch.frames.size() == 4, "include initial and final partial stride");
    const std::vector<std::uint64_t> expected{0, 6, 12, 17};
    for (std::size_t i = 0; i < expected.size(); ++i) {
        check(batch.frames[i].state.step == expected[i], "sample absolute step");
        close(batch.frames[i].state.time, expected[i] * spec.dt, 0.0, "sample physical time");
    }
    check(first.state().unwrapped_positions == batch.frames.back().state.unwrapped_positions,
          "observers do not consume RNG");
    check(batched.sample(0, 1, false).frames.empty(), "zero-step empty sample");
    check(batched.sample(0, 1, true).frames.size() == 1, "zero-step initial sample");
    rejects([&] { batched.sample(1, 0, true); }, "zero sampling stride rejected");
    spec.seed = 1235;
    bdmc::Simulation other(spec);
    other.advance(17);
    check(first.state().unwrapped_positions != other.state().unwrapped_positions, "different seeds differ");
    auto copy = first.state();
    copy.positions[0] = -123.0;
    check(first.state().positions[0] >= 0.0, "state snapshots own storage");
}

void uniform_brownian_moments() {
    // One-step independent coordinates: no wrapped-coordinate MSD ambiguity.
    auto spec = spec_1d(120000);
    spec.dt = 0.02;
    spec.diffusivity = 0.7;
    spec.seed = 81023;
    bdmc::Simulation simulation(spec);
    simulation.advance(1);
    const auto values = simulation.state().unwrapped_positions;
    double mean = 0.0, second = 0.0, fourth = 0.0;
    const double variance = 2.0 * spec.diffusivity * spec.dt;
    for (const auto value : values) {
        check(std::abs(value) <= std::sqrt(3.0 * variance), "uniform finite-step support");
        mean += value;
        second += value * value;
        fourth += value * value * value * value;
    }
    mean /= values.size(); second /= values.size(); fourth /= values.size();
    close(mean, 0.0, 5.0 * std::sqrt(variance / values.size()), "Brownian increment mean");
    close(second, variance, 0.015 * variance, "Brownian increment variance=2Ddt");
    close(fourth / (variance * variance), 1.8, 0.035, "finite-step uniform fourth moment differs from Gaussian");
    std::cout << "  Brownian mean=" << mean << " second=" << second << " expected=" << variance
              << " fourth/variance^2=" << fourth / (variance * variance) << '\n';
}

void gaussian_brownian_moments() {
    auto spec = spec_1d(100000);
    spec.noise_policy = bdmc::NoisePolicy::Gaussian;
    spec.dt = 0.015;
    spec.seed = 9191;
    bdmc::Simulation simulation(spec);
    simulation.advance(1);
    const auto values = simulation.state().unwrapped_positions;
    const double variance = 2.0 * spec.dt;
    double mean = 0.0, second = 0.0, fourth = 0.0;
    for (const auto value : values) {
        mean += value; second += value * value; fourth += value * value * value * value;
    }
    mean /= values.size(); second /= values.size(); fourth /= values.size();
    close(mean, 0.0, 5.0 * std::sqrt(variance / values.size()), "Gaussian increment mean");
    close(second, variance, 0.025 * variance, "Gaussian increment variance");
    close(fourth / (variance * variance), 3.0, 0.12, "Gaussian fourth moment");
}

void dumbbell_discrete_time_moments() {
    // Independent pairs in a large periodic box: no connector approaches L/2.
    // q'=a*q+eta, a=1-2*mu*k*dt, Var(eta)=4*D*dt.
    constexpr std::size_t pairs = 20000;
    constexpr std::uint64_t steps = 80;
    auto spec = spec_1d(2 * pairs);
    spec.model = bdmc::Model::Dumbbell;
    spec.spring_alpha = 1.25;
    spec.mobility = 0.8;
    spec.diffusivity = 0.6;
    spec.dt = 0.01;
    spec.seed = 71937;
    std::vector<double> initial(2 * pairs, 0.0);
    std::fill(initial.begin() + pairs, initial.end(), 1.0);
    bdmc::Simulation simulation(spec, initial);
    simulation.advance(steps);
    const auto state = simulation.state();
    const double k = spec.spring_alpha * spec.spring_alpha;
    const double a = 1.0 - 2.0 * spec.mobility * k * spec.dt;
    const double expected_mean = std::pow(a, static_cast<double>(steps));
    const double expected_variance = 4.0 * spec.diffusivity * spec.dt *
        (1.0 - std::pow(a, 2.0 * steps)) / (1.0 - a * a);
    double mean = 0.0, second = 0.0;
    for (std::size_t i = 0; i < pairs; ++i) {
        const double q = state.unwrapped_positions[i + pairs] - state.unwrapped_positions[i];
        mean += q;
        second += (q - expected_mean) * (q - expected_mean);
        check(std::abs(q) < 0.5 * spec.box_lengths[0], "large-box linear comparison stays in minimum-image branch");
    }
    mean /= pairs; second /= pairs;
    close(mean, expected_mean, 5.0 * std::sqrt(expected_variance / pairs), "two-moving-bead discrete relaxation");
    close(second, expected_variance, 0.045 * expected_variance, "discrete-time connector covariance");
    std::cout << "  Dumbbell mean=" << mean << " expected=" << expected_mean
              << " variance=" << second << " expected=" << expected_variance << '\n';
}

void failed_force_does_not_advance() {
    auto spec = spec_1d();
    bdmc::Simulation simulation(spec);
    simulation.add_force_term([](const bdmc::ParticleState&, const bdmc::OrthorhombicBox&,
                                 const bdmc::FixedTopology&, const bdmc::DynamicTopology&,
                                 bdmc::ForceResult& result) {
        result.values[0] = std::numeric_limits<double>::quiet_NaN();
    });
    rejects([&] { simulation.advance(1); }, "nonfinite component force rejected");
    check(simulation.state().step == 0, "failed force does not advance positions/time");
}

}  // namespace

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"geometry", geometry},
        {"specifications", specifications},
        {"identity_and_topology", identity_and_topology},
        {"harmonic_force_and_old_state_update", harmonic_force_and_old_state_update},
        {"wrapping_and_component_composition", wrapping_and_component_composition},
        {"seeded_stream_and_observation_schedule", seeded_stream_and_observation_schedule},
        {"uniform_brownian_moments", uniform_brownian_moments},
        {"gaussian_brownian_moments", gaussian_brownian_moments},
        {"dumbbell_discrete_time_moments", dumbbell_discrete_time_moments},
        {"failed_force_does_not_advance", failed_force_does_not_advance},
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
    std::cout << tests.size() - failures << '/' << tests.size() << " native check groups passed\n";
    return failures == 0 ? 0 : 1;
}
