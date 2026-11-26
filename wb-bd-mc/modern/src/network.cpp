#include "bdmc/network.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace bdmc {
namespace {
void finite(double x, const char* name) {
    if (!std::isfinite(x)) throw std::invalid_argument(std::string(name) + " must be finite");
}
void nonnegative(double x, const char* name) {
    if (!std::isfinite(x) || x < 0) throw std::invalid_argument(std::string(name) + " must be finite and nonnegative");
}
void positive(double x, const char* name) {
    if (!std::isfinite(x) || x <= 0) throw std::invalid_argument(std::string(name) + " must be finite and positive");
}
std::map<std::size_t, double> checked_map(std::size_t n, const std::vector<NetworkWeight>& row,
                                       bool require_nonnegative) {
    std::map<std::size_t, double> result;
    for (const auto& item : row) {
        if (item.particle >= n) throw std::invalid_argument("network weight particle ID out of range");
        if (require_nonnegative) nonnegative(item.weight, "network weight");
        else finite(item.weight, "network tension");
        if (!result.emplace(item.particle, item.weight).second)
            throw std::invalid_argument("duplicate particle in network weight map");
    }
    return result;
}
double lookup(const std::map<std::size_t, double>& values, std::size_t p) {
    const auto found = values.find(p);
    if (found == values.end()) throw std::invalid_argument("network map omits an existing adjacency entry");
    return found->second;
}
} // namespace

NetworkTopology::NetworkTopology(std::size_t n_particles, std::vector<ChainEndpoints> endpoints)
    : n_particles_(n_particles), endpoints_(std::move(endpoints)) {
    if (!n_particles_) throw std::invalid_argument("network requires at least one particle");
    if (n_chains() > std::numeric_limits<std::size_t>::max() / 2)
        throw std::invalid_argument("network end-ID range overflows");
    rebuild();
}
NetworkTopology NetworkTopology::loops(std::size_t n, std::size_t per_particle) {
    if (!n || per_particle > std::numeric_limits<std::size_t>::max() / n / 2)
        throw std::invalid_argument("invalid loop initialization size");
    std::vector<ChainEndpoints> endpoints;
    endpoints.reserve(n * per_particle);
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t c = 0; c < per_particle; ++c) endpoints.push_back({p, p});
    return NetworkTopology(n, std::move(endpoints));
}
void NetworkTopology::require_particle(std::size_t p) const {
    if (p >= n_particles_) throw std::invalid_argument("network particle ID out of range");
}
std::size_t NetworkTopology::particle_of_end(std::size_t end) const {
    if (end >= 2 * n_chains()) throw std::invalid_argument("network end ID out of range");
    return endpoints_[end % n_chains()][end / n_chains()];
}
std::size_t NetworkTopology::opposite_end(std::size_t end) const {
    if (end >= 2 * n_chains()) throw std::invalid_argument("network end ID out of range");
    return end < n_chains() ? end + n_chains() : end - n_chains();
}
std::size_t NetworkTopology::weight(std::size_t anchor, std::size_t other) const {
    require_particle(anchor); require_particle(other);
    for (const auto& item : adjacency_[anchor]) if (item.particle == other) return item.weight;
    return 0;
}
std::vector<std::size_t> NetworkTopology::compatible_ends(std::size_t anchor, std::size_t other) const {
    require_particle(anchor); require_particle(other);
    // This full scan is the transparent reference for a future incident-end
    // index. Preserve exact labels and both loop ends when optimizing it; bridge
    // multiplicity alone cannot recover which physical end moved.
    std::vector<std::size_t> ends;
    for (std::size_t end = 0; end < 2 * n_chains(); ++end)
        if (particle_of_end(end) == other && particle_of_end(opposite_end(end)) == anchor)
            ends.push_back(end);
    return ends;
}
void NetworkTopology::rebuild() {
    // Endpoint rows are the source of truth; regenerate all aggregate views
    // rather than letting graph counts drift away from labeled-chain history.
    std::vector<std::map<std::size_t, std::size_t>> counts(n_particles_);
    attached_ends_.assign(n_particles_, 0);
    n_bridges_ = 0;
    for (std::size_t p = 0; p < n_particles_; ++p) counts[p][p] = 0;
    for (const auto& ends : endpoints_) {
        const auto a = ends[0], b = ends[1]; require_particle(a); require_particle(b);
        ++attached_ends_[a]; ++attached_ends_[b];
        if (a == b) counts[a][a] += 2;
        else { ++counts[a][b]; ++counts[b][a]; ++n_bridges_; }
    }
    adjacency_.assign(n_particles_, {}); n_distinct_pairs_ = 0;
    for (std::size_t p = 0; p < n_particles_; ++p) {
        // Permanent self first, then target-ID order. This preserves the
        // discrete measure, not historical insertion order or RNG trajectories.
        adjacency_[p].push_back({p, counts[p][p]});
        for (const auto& item : counts[p]) if (item.first != p) {
            adjacency_[p].push_back({item.first, item.second});
            if (p < item.first) ++n_distinct_pairs_;
        }
    }
}
void NetworkTopology::transact(const std::vector<EndpointMove>& moves) {
    if (moves.empty()) return;
    if (epoch_ == std::numeric_limits<std::uint64_t>::max()) throw std::overflow_error("topology epoch overflow");
    auto next = endpoints_;
    bool changed = false;
    for (const auto& move : moves) {
        require_particle(move.expected_particle); require_particle(move.destination);
        if (move.end >= 2 * n_chains()) throw std::invalid_argument("network transaction end ID out of range");
        auto& attached = next[move.end % n_chains()][move.end / n_chains()];
        if (attached != move.expected_particle) throw std::invalid_argument("network transaction expectation mismatch");
        changed = changed || attached != move.destination;
        attached = move.destination;
    }
    if (!changed || next == endpoints_) return;
    NetworkTopology candidate(n_particles_, std::move(next));
    // No throwing operation follows the first swap; failed allocation above
    // therefore also leaves all persistent representations unchanged.
    endpoints_.swap(candidate.endpoints_); attached_ends_.swap(candidate.attached_ends_);
    adjacency_.swap(candidate.adjacency_);
    n_bridges_ = candidate.n_bridges_; n_distinct_pairs_ = candidate.n_distinct_pairs_; ++epoch_;
}
void NetworkTopology::validate() const {
    NetworkTopology expected(n_particles_, endpoints_);
    if (attached_ends_ != expected.attached_ends_ || n_bridges_ != expected.n_bridges_ ||
        n_distinct_pairs_ != expected.n_distinct_pairs_ || adjacency_.size() != expected.adjacency_.size())
        throw std::logic_error("network endpoint/aggregate invariant failed");
    std::size_t total = 0;
    for (std::size_t p = 0; p < n_particles_; ++p) {
        if (adjacency_[p].size() != expected.adjacency_[p].size()) throw std::logic_error("network adjacency size mismatch");
        for (std::size_t k = 0; k < adjacency_[p].size(); ++k) {
            const auto& a = adjacency_[p][k]; const auto& b = expected.adjacency_[p][k];
            if (a.particle != b.particle || a.weight != b.weight) throw std::logic_error("network adjacency mismatch");
            total += a.weight;
        }
    }
    if (total != 2 * n_chains()) throw std::logic_error("network global endpoint conservation failed");
}

std::string network_action_name(NetworkAction action) {
    switch (action) {
    case NetworkAction::Cancel: return "cancel";
    case NetworkAction::Add: return "add";
    case NetworkAction::OppositeDelete: return "opposite_delete";
    case NetworkAction::Move: return "move";
    }
    throw std::invalid_argument("invalid network action");
}
std::string network_rejection_name(NetworkRejection reason) {
    switch (reason) {
    case NetworkRejection::None: return "none";
    case NetworkRejection::SameDestination: return "same_destination";
    case NetworkRejection::CountOrMultiplicity: return "count_or_multiplicity";
    case NetworkRejection::Acceptance: return "acceptance";
    }
    throw std::invalid_argument("invalid network rejection");
}
void NetworkConstraints::validate() const {
    if (min_attached_ends > max_attached_ends) throw std::invalid_argument("network count bounds reversed");
}
NetworkClassification classify_network_action(const NetworkTopology& topology, std::size_t anchor,
    std::size_t other, std::size_t destination, const NetworkConstraints& constraints) {
    constraints.validate();
    if (anchor >= topology.n_particles() || other >= topology.n_particles() || destination >= topology.n_particles())
        throw std::invalid_argument("network proposal particle ID out of range");
    if (!topology.weight(anchor, other)) throw std::invalid_argument("network proposal selects an absent old connector");
    if (other == destination) return {NetworkAction::Cancel, NetworkRejection::SameDestination};
    NetworkAction action = anchor == other ? NetworkAction::Add :
        (anchor == destination ? NetworkAction::OppositeDelete : NetworkAction::Move);
    bool allowed = true;
    if (constraints.multiple_connections) {
        // The opposite endpoint, never the fixed anchor end, is moved.
        const auto& ends = topology.attached_ends();
        allowed = ends[other] > constraints.min_attached_ends && ends[destination] < constraints.max_attached_ends;
    } else if (action != NetworkAction::OppositeDelete) {
        allowed = topology.weight(anchor, destination) == 0;
    }
    if (!allowed) return {NetworkAction::Cancel, NetworkRejection::CountOrMultiplicity};
    return {action, NetworkRejection::None};
}

void NetworkKinetics::validate() const {
    nonnegative(l_cap, "l_cap"); finite(energy_barrier, "energy_barrier"); nonnegative(delta_t0, "delta_t0");
    switch (law) {
    case NetworkTransitionLaw::Uniform: case NetworkTransitionLaw::Metropolis:
    case NetworkTransitionLaw::Weighted: case NetworkTransitionLaw::Dissociation:
    case NetworkTransitionLaw::FirstOrder: return;
    }
    throw std::invalid_argument("unsupported network transition law");
}
double network_transition_probability(double tension, const NetworkKinetics& kinetics) {
    kinetics.validate(); finite(tension, "old tension");
    if (kinetics.law == NetworkTransitionLaw::Uniform) return 1;
    const double exponent = tension * kinetics.l_cap;
    finite(exponent, "tension*l_cap");
    if (kinetics.law == NetworkTransitionLaw::Metropolis) {
        const double barrier_exponent = exponent - kinetics.energy_barrier;
        return barrier_exponent >= 0 ? 1 : std::exp(barrier_exponent);
    }
    if (kinetics.law == NetworkTransitionLaw::Weighted) return exponent >= 0 ? 1 : std::exp(exponent);
    // Keep the legacy rate arithmetic within its supported finite domain.
    // A huge rate times a tiny delta can have a finite hazard; allowing an
    // overflowing intermediate to silently saturate would misrepresent it.
    if (kinetics.delta_t0 == 0) return 0;
    const double rate = std::exp(exponent);
    if (!std::isfinite(rate) || rate == 0)
        throw std::invalid_argument("kinetic rate exponential is outside the positive finite supported domain");
    const double hazard = kinetics.delta_t0 * rate;
    if (kinetics.law == NetworkTransitionLaw::Dissociation) return std::min(1.0, hazard);
    // expm1 retains tiny hazards; same mathematical FIRST_ORDER law.
    return -std::expm1(-hazard);
}
double network_detachment_weight(double tension, double l_cap, NetworkDetachmentLaw law) {
    finite(tension, "old tension"); nonnegative(l_cap, "l_cap");
    if (law == NetworkDetachmentLaw::Uniform) return 1;
    if (law != NetworkDetachmentLaw::Weighted) throw std::invalid_argument("unsupported detachment law");
    const double value = std::exp(tension * l_cap);
    if (!std::isfinite(value)) throw std::invalid_argument("detachment weight overflows; no implicit reweighting");
    return value;
}
void NetworkSchedule::validate() const {
    positive(dt_over_tauR, "dt/tauR"); positive(repulsion_coefficient, "repulsion coefficient"); positive(Rt, "Rt");
    if (!steps_per_block) throw std::invalid_argument("steps_per_block must be positive");
}
double NetworkSchedule::delta_t0() const {
    validate();
    // Block length enters the attempt hazard. Treating it only as a recording
    // stride would change the physical rate even with unchanged BD dt.
    const double value = dt_over_tauR * static_cast<double>(steps_per_block) / (repulsion_coefficient * Rt);
    if (!std::isfinite(value) || value <= 0) throw std::invalid_argument("network rate conversion is not finite positive");
    return value;
}

std::vector<NetworkProbability> network_discrete_distribution(std::size_t n, const std::vector<NetworkWeight>& raw) {
    if (!n) throw std::invalid_argument("distribution requires at least one particle");
    checked_map(n, raw, true);
    double total = 0;
    for (const auto& item : raw) total += item.weight;
    positive(total, "network normalization");
    std::vector<NetworkProbability> result; result.reserve(raw.size());
    for (const auto& item : raw) result.push_back({item.particle, item.weight / total});
    return result;
}
std::vector<NetworkProbability> network_suggestion_distribution(const NetworkTopology& topology, std::size_t anchor,
    const std::vector<NetworkWeight>& detachment) {
    if (anchor >= topology.n_particles()) throw std::invalid_argument("suggestion anchor ID out of range");
    const auto weights = checked_map(topology.n_particles(), detachment, true);
    std::vector<NetworkWeight> raw;
    for (const auto& item : topology.adjacency()[anchor])
        raw.push_back({item.particle, static_cast<double>(item.weight) * lookup(weights, item.particle)});
    return network_discrete_distribution(topology.n_particles(), raw);
}
std::vector<NetworkProposal> enumerate_network_proposals(const NetworkTopology& topology,
    const std::vector<std::vector<NetworkWeight>>& association,
    const std::vector<std::vector<NetworkWeight>>& detachment,
    const std::vector<std::vector<NetworkWeight>>& tensions,
    const NetworkConstraints& constraints, const NetworkKinetics& kinetics) {
    const auto n = topology.n_particles(); constraints.validate(); kinetics.validate();
    if (association.size() != n || detachment.size() != n || tensions.size() != n)
        throw std::invalid_argument("proposal maps need one row per stable particle");
    const auto nonempty = static_cast<std::size_t>(std::count_if(topology.attached_ends().begin(),
        topology.attached_ends().end(), [](auto count) { return count != 0; }));
    if (!nonempty) throw std::invalid_argument("no nonempty network anchor exists");
    // Enumerate the complete finite proposal measure for inspection/testing,
    // including cancellations. Resampling until mutation would condition away
    // rejected mass and change the stochastic clock; no RNG is drawn here.
    std::vector<NetworkProposal> result;
    for (std::size_t anchor = 0; anchor < n; ++anchor) {
        if (!topology.attached_ends()[anchor]) continue;
        const auto destinations = network_discrete_distribution(n, association[anchor]);
        const auto suggestions = network_suggestion_distribution(topology, anchor, detachment[anchor]);
        const auto old_tensions = checked_map(n, tensions[anchor], false);
        for (const auto& old : suggestions) {
            if (old.probability == 0) continue;
            const double acceptance = network_transition_probability(lookup(old_tensions, old.particle), kinetics);
            const auto ends = topology.compatible_ends(anchor, old.particle);
            if (ends.empty()) throw std::logic_error("suggestion without compatible labeled endpoint");
            for (const auto& dest : destinations) {
                if (dest.probability == 0) continue;
                const auto classified = classify_network_action(topology, anchor, old.particle, dest.particle, constraints);
                const double mass = old.probability * dest.probability /
                    static_cast<double>(nonempty) / static_cast<double>(ends.size());
                for (const auto end : ends) result.push_back({anchor, old.particle, dest.particle, end,
                    classified.action, classified.rejection, mass, acceptance,
                    classified.action == NetworkAction::Cancel ? 0 : mass * acceptance, topology.epoch()});
            }
        }
    }
    return result;
}
NetworkAttemptResult attempt_network_transition(NetworkTopology& topology, const NetworkProposal& proposal,
    double tension, double roll, const NetworkConstraints& constraints, const NetworkKinetics& kinetics) {
    if (proposal.topology_epoch != topology.epoch()) throw std::invalid_argument("stale network proposal epoch");
    if (!std::isfinite(roll) || roll < 0 || roll >= 1) throw std::invalid_argument("acceptance roll must be in [0,1)");
    const auto classified = classify_network_action(topology, proposal.anchor, proposal.other, proposal.destination, constraints);
    const auto ends = topology.compatible_ends(proposal.anchor, proposal.other);
    if (std::find(ends.begin(), ends.end(), proposal.end) == ends.end())
        throw std::invalid_argument("selected labeled end is incompatible with old connector");
    const double p = network_transition_probability(tension, kinetics);
    // Production tests acceptance before action classification. This determines
    // rejection reporting when the same proposed destination also violates a rule.
    if (!(roll < p)) return {NetworkAction::Cancel, NetworkRejection::Acceptance, p};
    if (classified.action == NetworkAction::Cancel) return {classified.action, classified.rejection, p};
    topology.transact({EndpointMove{proposal.end, proposal.other, proposal.destination}});
    return {classified.action, NetworkRejection::None, p};
}

void NetworkConnectorSpec::validate() const {
    nonnegative(spring_constant, "network spring constant"); positive(q0, "network q0"); nonnegative(scale, "network force/energy scale");
    if (law != NetworkConnectorLaw::Hookean && law != NetworkConnectorLaw::FENE)
        throw std::invalid_argument("unsupported network connector law");
}
NetworkConnectorFrame evaluate_network_connectors(const NetworkTopology& topology,
    const std::vector<double>& positions, const std::vector<double>& lengths, const NetworkConnectorSpec& spec) {
    spec.validate();
    if (topology.n_particles() > std::numeric_limits<std::size_t>::max() / 3 ||
        positions.size() != topology.n_particles() * 3 || lengths.size() != 3)
        throw std::invalid_argument("network connector geometry requires N by 3 positions and three box lengths");
    for (const auto value : positions) finite(value, "network position");
    OrthorhombicBox box(lengths);
    NetworkConnectorFrame result; result.forces.assign(positions.size(), 0);
    for (std::size_t i = 0; i < topology.n_particles(); ++i) for (const auto& edge : topology.adjacency()[i]) {
        const auto j = edge.particle; if (j <= i) continue;
        std::array<double, 3> q{}; double r2 = 0;
        for (std::size_t c = 0; c < 3; ++c) { q[c] = box.minimum_image(positions[j*3+c] - positions[i*3+c], c); r2 += q[c]*q[c]; }
        if (!std::isfinite(r2)) throw std::invalid_argument("network bonded displacement overflows");
        double factor = spec.spring_constant;
        double energy = 0.5 * spec.spring_constant * r2;
        if (spec.law == NetworkConnectorLaw::FENE) {
            const double ratio = std::sqrt(r2) / spec.q0;
            if (!(ratio < 1)) throw std::invalid_argument("network FENE bridge must satisfy r<q0");
            const double x = ratio * ratio;
            factor /= 1 - x;
            energy = -0.5 * spec.spring_constant * spec.q0 * spec.q0 * std::log1p(-x);
        }
        const double multiplicity = static_cast<double>(edge.weight) * spec.scale;
        factor *= multiplicity; energy *= multiplicity;
        if (!std::isfinite(factor) || !std::isfinite(energy)) throw std::invalid_argument("network connector output overflows");
        result.potential_energy += energy;
        for (std::size_t a = 0; a < 3; ++a) {
            const double f = factor * q[a]; result.forces[i*3+a] += f; result.forces[j*3+a] -= f;
            for (std::size_t b = 0; b < 3; ++b) result.connector_virial[a*3+b] += q[a]*factor*q[b];
        }
    }
    if (!std::isfinite(result.potential_energy)) throw std::invalid_argument("network energy accumulation overflows");
    for (const auto value : result.forces) finite(value, "network force output");
    for (const auto value : result.connector_virial) finite(value, "network connector virial output");
    return result;
}

} // namespace bdmc
