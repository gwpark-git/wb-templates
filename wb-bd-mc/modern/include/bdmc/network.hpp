#pragma once

#include "bdmc/core.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace bdmc {

// Chain IDs are row indices. End IDs 0..Nc-1 are HEAD and Nc..2Nc-1
// are TAIL, matching CHAIN_HANDLE. Particle labels never refer to cell slots.
using ChainEndpoints = std::array<std::size_t, 2>;
struct NetworkAdjacency {
    std::size_t particle = 0;
    // Self: twice the loop count. Other: bridge multiplicity.
    std::size_t weight = 0;
};
struct EndpointMove {
    std::size_t end = 0, expected_particle = 0, destination = 0;
};

// Labeled endpoints are authoritative; adjacency is a derived multiplicity
// view. Two bridges joining the same particles remain different chains, and a
// loop retains two distinct ends. This correctness reference rebuilds aggregates
// on transaction; an incremental index must preserve the same labeled events.
class NetworkTopology {
public:
    NetworkTopology(std::size_t n_particles, std::vector<ChainEndpoints> endpoints);
    static NetworkTopology loops(std::size_t n_particles, std::size_t chains_per_particle);
    std::size_t n_particles() const noexcept { return n_particles_; }
    std::size_t n_chains() const noexcept { return endpoints_.size(); }
    std::uint64_t epoch() const noexcept { return epoch_; }
    const std::vector<ChainEndpoints>& endpoints() const noexcept { return endpoints_; }
    const std::vector<std::size_t>& attached_ends() const noexcept { return attached_ends_; }
    const std::vector<std::vector<NetworkAdjacency>>& adjacency() const noexcept { return adjacency_; }
    std::size_t particle_of_end(std::size_t end) const;
    std::size_t opposite_end(std::size_t end) const;
    std::size_t weight(std::size_t anchor, std::size_t other) const;
    std::size_t n_bridges() const noexcept { return n_bridges_; }
    std::size_t n_distinct_pairs() const noexcept { return n_distinct_pairs_; }
    // Ends attached at 'other' whose opposite end stays at 'anchor'. Loops
    // supply both labeled ends; a bridge supplies one end per physical chain.
    std::vector<std::size_t> compatible_ends(std::size_t anchor, std::size_t other) const;
    // Atomic batch, interpreted in order; any invalid expectation/ID rolls
    // back the whole batch. This primitive imposes no scientific count policy.
    void transact(const std::vector<EndpointMove>& moves);
    void validate() const;
private:
    void rebuild();
    void require_particle(std::size_t particle) const;
    std::size_t n_particles_;
    std::vector<ChainEndpoints> endpoints_;
    std::vector<std::size_t> attached_ends_;
    std::vector<std::vector<NetworkAdjacency>> adjacency_;
    std::size_t n_bridges_ = 0, n_distinct_pairs_ = 0;
    std::uint64_t epoch_ = 0;
};

enum class NetworkAction { Cancel = 0, Add = 1, OppositeDelete = 2, Move = 3 };
enum class NetworkRejection { None, SameDestination, CountOrMultiplicity, Acceptance };
std::string network_action_name(NetworkAction);
std::string network_rejection_name(NetworkRejection);

// Exact bounded subset of TRUTH_MAP: MULTIPLE with endpoint-count bounds,
// or SINGLE (ignores count bounds, as the original SINGLE checks do).
struct NetworkConstraints {
    bool multiple_connections = true;
    std::size_t min_attached_ends = 0;
    std::size_t max_attached_ends = std::numeric_limits<std::size_t>::max();
    void validate() const;
};
struct NetworkClassification {
    NetworkAction action = NetworkAction::Cancel;
    NetworkRejection rejection = NetworkRejection::None;
};
NetworkClassification classify_network_action(const NetworkTopology&, std::size_t anchor,
    std::size_t other, std::size_t destination, const NetworkConstraints&);

enum class NetworkTransitionLaw { Uniform, Metropolis, Weighted, Dissociation, FirstOrder };
enum class NetworkDetachmentLaw { Uniform, Weighted };
struct NetworkKinetics {
    NetworkTransitionLaw law = NetworkTransitionLaw::Uniform;
    double l_cap = 0.0, energy_barrier = 0.0, delta_t0 = 0.0;
    void validate() const;
};
// Inputs are the tension of the OLD connector, before the proposed change.
// For positive delta_t0, dissociation/first-order require exp(tension*l_cap)
// to be positive finite. Overflow/underflow intermediates fail explicitly;
// no extension of the legacy law to those numerical extremes is claimed.
double network_transition_probability(double old_tension, const NetworkKinetics&);
double network_detachment_weight(double old_tension, double l_cap, NetworkDetachmentLaw);

// N_steps_block is part of the rate conversion and not just output cadence.
// Keep the kinetic clock explicit when changing batching: the same number of
// BD steps with a different block schedule is not automatically the same MC law.
struct NetworkSchedule {
    double dt_over_tauR = 0.001, repulsion_coefficient = 1.0, Rt = 1.0;
    std::uint64_t steps_per_block = 1;
    void validate() const;
    double delta_t0() const;
    std::size_t attempts_per_block(const NetworkTopology& topology) const noexcept {
        return topology.n_chains();
    }
};

struct NetworkWeight { std::size_t particle = 0; double weight = 0.0; };
struct NetworkProbability { std::size_t particle = 0; double probability = 0.0; };
// PROBABILITY_AREA=false only. Unique particle IDs and finite nonnegative
// raw values required; all-zero maps fail explicitly. No radial Jacobian.
std::vector<NetworkProbability> network_discrete_distribution(
    std::size_t n_particles, const std::vector<NetworkWeight>& raw);
std::vector<NetworkProbability> network_suggestion_distribution(
    const NetworkTopology&, std::size_t anchor, const std::vector<NetworkWeight>& detachment);

struct NetworkProposal {
    std::size_t anchor = 0, other = 0, destination = 0, end = 0;
    NetworkAction action = NetworkAction::Cancel;
    NetworkRejection rejection = NetworkRejection::None;
    double proposal_probability = 0.0, acceptance_probability = 0.0;
    double mutation_probability = 0.0;
    std::uint64_t topology_epoch = 0;
};
// Exact finite measure for ordinary SERIAL attempts at frozen geometry:
// uniform nonempty anchor -> weighted old neighbour -> discrete destination
// -> uniform compatible opposite end. No RNG is owned or consumed here.
// Each outer row corresponds to a stable particle ID. Empty anchors' rows
// may be empty. Old tensions rows cover all adjacency entries including self.
// Rejected proposal mass remains present; it is never resampled away.
std::vector<NetworkProposal> enumerate_network_proposals(const NetworkTopology&,
    const std::vector<std::vector<NetworkWeight>>& association,
    const std::vector<std::vector<NetworkWeight>>& detachment,
    const std::vector<std::vector<NetworkWeight>>& old_tensions,
    const NetworkConstraints&, const NetworkKinetics&);

struct NetworkAttemptResult {
    NetworkAction action = NetworkAction::Cancel;
    NetworkRejection rejection = NetworkRejection::None;
    double acceptance_probability = 0.0;
};
// Caller supplies a proposal's choices and a uniform roll in [0,1). The
// compatible end is explicit; neither constraint nor rate rejection retries.
// Stale proposals fail before mutation. No claim of legacy RNG identity.
NetworkAttemptResult attempt_network_transition(NetworkTopology&, const NetworkProposal&,
    double old_tension, double acceptance_roll, const NetworkConstraints&, const NetworkKinetics&);

enum class NetworkConnectorLaw { Hookean, FENE };
struct NetworkConnectorSpec {
    NetworkConnectorLaw law = NetworkConnectorLaw::Hookean;
    double spring_constant = 3.0, q0 = 10.0, scale = 1.0;
    void validate() const;
};
struct NetworkConnectorFrame {
    std::vector<double> forces;
    double potential_energy = 0.0;
    // Once per physical bridge, positive tension q outer F_tension.
    // No volume normalisation or random/repulsive/solvent contribution.
    std::array<double, 9> connector_virial{};
};
// 3D orthorhombic PBC, independent of all spatial-neighbour caches. Loops
// contribute no interparticle force/energy. FENE requires every bridge r<q0.
NetworkConnectorFrame evaluate_network_connectors(const NetworkTopology&,
    const std::vector<double>& positions, const std::vector<double>& box_lengths,
    const NetworkConnectorSpec&);

} // namespace bdmc
