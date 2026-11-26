#include "bdmc/network.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <limits>

namespace py = pybind11;
namespace {
py::ssize_t size(std::size_t n) {
    if (n > static_cast<std::size_t>(std::numeric_limits<py::ssize_t>::max())) throw std::overflow_error("network array size overflows ssize_t");
    return static_cast<py::ssize_t>(n);
}
template<class T> py::array_t<T> copied(const std::vector<T>& values, const std::vector<py::ssize_t>& shape) {
    py::array_t<T> output(shape); std::copy(values.begin(), values.end(), output.mutable_data()); return output;
}
std::vector<bdmc::ChainEndpoints> endpoints(py::object object) {
    const auto a = py::array::ensure(object);
    if (!a || a.ndim() != 2 || a.shape(1) != 2 || (a.dtype().kind() != 'i' && a.dtype().kind() != 'u'))
        throw std::invalid_argument("endpoints require an integer array of shape (n_chains,2)");
    const auto values = py::array_t<std::uint64_t, py::array::c_style | py::array::forcecast>::ensure(a);
    if (a.dtype().kind() == 'i') {
        const auto signed_values = py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>::ensure(a);
        for (py::ssize_t i = 0; i < a.size(); ++i)
            if (signed_values.data()[i] < 0) throw std::invalid_argument("endpoint particle IDs must be nonnegative");
    }
    std::vector<bdmc::ChainEndpoints> output; output.reserve(static_cast<std::size_t>(a.shape(0)));
    for (py::ssize_t i = 0; i < a.shape(0); ++i) {
        const auto first = values.data()[2*i], second = values.data()[2*i+1];
        if (first > std::numeric_limits<std::size_t>::max() || second > std::numeric_limits<std::size_t>::max())
            throw std::invalid_argument("endpoint particle IDs overflow size_t");
        output.push_back({static_cast<std::size_t>(first), static_cast<std::size_t>(second)});
    }
    return output;
}
std::vector<double> positions(py::object object, std::size_t count) {
    const auto a = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(object);
    if (!a || a.ndim() != 2 || a.shape(0) != size(count) || a.shape(1) != 3)
        throw std::invalid_argument("positions require an array of shape (n_particles,3)");
    return {a.data(), a.data() + a.size()};
}
py::dict metadata() {
    py::dict d;
    d["schema"] = "bdmc.network.components.v1";
    d["chain_end_ids"] = "0..Nc-1 HEAD; Nc..2Nc-1 TAIL; stable chain row IDs";
    d["adjacency"] = "CSR; permanent self first; self weight=2*loops; other weight=bridge multiplicity; sorted other IDs";
    d["ownership"] = "all NumPy results are independent C-contiguous owned copies";
    d["proposal_ownership"] = "Python proposal choices and cached masses are readonly; explicit manual choices carry no evaluated proposal mass";
    d["proposal_scope"] = "ordinary serial frozen-geometry attempts; PROBABILITY_AREA=false; caller supplies raw maps and old tensions";
    d["proposal_measure"] = "uniform nonempty anchor; multiplicity*detachment old-neighbor weight; discrete raw association weight; uniform compatible opposite labeled end";
    d["rng"] = "none; caller supplies explicit choices and acceptance roll; no legacy RNG sequence identity";
    d["kinetic_timing"] = "delta_t0=(dt/tauR)*steps_per_block/(C*Rt); N_chains ordinary attempts per legacy block";
    d["supported_constraints"] = "MULTIPLE endpoint-count checks or SINGLE bridge checks; SINGLE ignores count bounds as legacy does";
    d["unsupported"] = "full MC/BD controller; cutoff-destruction; PROBABILITY_AREA; renewal; parallel locks; checkpoint/restart; network stationarity or rheology validation";
    return d;
}
py::dict snapshot(const bdmc::NetworkTopology& topology) {
    // CSR is an exchange view of aggregate adjacency, not the authoritative
    // labeled-chain state. Retain endpoints alongside it so loop/bridge identity
    // survives a Python round trip or an alternative graph implementation.
    std::vector<std::uint64_t> ends, attached, ids, weights, offsets{0};
    ends.reserve(topology.n_chains()*2);
    for (const auto& row : topology.endpoints()) { ends.push_back(row[0]); ends.push_back(row[1]); }
    for (const auto value : topology.attached_ends()) attached.push_back(value);
    for (const auto& row : topology.adjacency()) {
        for (const auto& item : row) { ids.push_back(item.particle); weights.push_back(item.weight); }
        offsets.push_back(ids.size());
    }
    py::dict d;
    d["endpoints"] = copied(ends, {size(topology.n_chains()), 2});
    d["attached_ends"] = copied(attached, {size(topology.n_particles())});
    d["adjacency_particles"] = copied(ids, {size(ids.size())});
    d["adjacency_weights"] = copied(weights, {size(weights.size())});
    d["adjacency_offsets"] = copied(offsets, {size(offsets.size())});
    d["n_bridges"] = topology.n_bridges(); d["n_distinct_pairs"] = topology.n_distinct_pairs();
    d["epoch"] = topology.epoch(); d["metadata"] = metadata(); return d;
}
bdmc::NetworkTransitionLaw transition_law(const std::string& name) {
    if (name == "uniform") return bdmc::NetworkTransitionLaw::Uniform;
    if (name == "metropolis") return bdmc::NetworkTransitionLaw::Metropolis;
    if (name == "weighted") return bdmc::NetworkTransitionLaw::Weighted;
    if (name == "dissociation") return bdmc::NetworkTransitionLaw::Dissociation;
    if (name == "first_order") return bdmc::NetworkTransitionLaw::FirstOrder;
    throw std::invalid_argument("unsupported network transition law");
}
bdmc::NetworkDetachmentLaw detachment_law(const std::string& name) {
    if (name == "uniform") return bdmc::NetworkDetachmentLaw::Uniform;
    if (name == "weighted") return bdmc::NetworkDetachmentLaw::Weighted;
    throw std::invalid_argument("unsupported network detachment law");
}
using MapRow = std::vector<std::pair<std::size_t, double>>;
using MapRows = std::vector<MapRow>;
std::vector<bdmc::NetworkWeight> weights(const MapRow& raw) {
    std::vector<bdmc::NetworkWeight> out; out.reserve(raw.size());
    for (const auto& item : raw) out.push_back({item.first, item.second}); return out;
}
std::vector<std::vector<bdmc::NetworkWeight>> weights(const MapRows& raw) {
    std::vector<std::vector<bdmc::NetworkWeight>> out; out.reserve(raw.size());
    for (const auto& row : raw) out.push_back(weights(row)); return out;
}
py::dict distribution(const std::vector<bdmc::NetworkProbability>& raw) {
    std::vector<std::uint64_t> ids; std::vector<double> values;
    for (const auto& item : raw) { ids.push_back(item.particle); values.push_back(item.probability); }
    py::dict d; d["particles"] = copied(ids, {size(ids.size())}); d["probabilities"] = copied(values, {size(values.size())}); return d;
}
} // namespace

PYBIND11_MODULE(_network, m) {
    m.doc() = "Owned labeled transient-network topology and legacy-mapped serial proposal primitives";
    py::class_<bdmc::NetworkTopology>(m, "Topology")
        .def(py::init([](std::size_t n, py::object e) { return bdmc::NetworkTopology(n, endpoints(e)); }), py::arg("n_particles"), py::arg("endpoints"))
        .def_static("loops", &bdmc::NetworkTopology::loops, py::arg("n_particles"), py::arg("chains_per_particle"))
        .def_property_readonly("n_particles", &bdmc::NetworkTopology::n_particles)
        .def_property_readonly("n_chains", &bdmc::NetworkTopology::n_chains)
        .def_property_readonly("epoch", &bdmc::NetworkTopology::epoch)
        .def("snapshot", &snapshot).def("validate", &bdmc::NetworkTopology::validate)
        .def("particle_of_end", &bdmc::NetworkTopology::particle_of_end, py::arg("end"))
        .def("opposite_end", &bdmc::NetworkTopology::opposite_end, py::arg("end"))
        .def("compatible_ends", [](const bdmc::NetworkTopology& t, std::size_t a, std::size_t o) {
            const auto ends = t.compatible_ends(a, o); std::vector<std::uint64_t> values(ends.begin(), ends.end()); return copied(values, {size(values.size())});
        }, py::arg("anchor"), py::arg("other"))
        .def("transact", [](bdmc::NetworkTopology& t, const std::vector<std::array<std::size_t,3>>& input) {
            std::vector<bdmc::EndpointMove> moves; moves.reserve(input.size());
            for (const auto& v : input) moves.push_back({v[0],v[1],v[2]}); t.transact(moves);
        }, py::arg("moves"), "Atomic list of (end_id, expected_particle, destination), interpreted in order");
    py::class_<bdmc::NetworkConstraints>(m, "Constraints")
        .def(py::init([](bool multiple, std::size_t lower, std::size_t upper) {
            bdmc::NetworkConstraints c; c.multiple_connections=multiple; c.min_attached_ends=lower; c.max_attached_ends=upper; c.validate(); return c;
        }), py::arg("multiple_connections")=true, py::arg("min_attached_ends")=0,
            py::arg("max_attached_ends")=std::numeric_limits<std::size_t>::max())
        .def_readonly("multiple_connections", &bdmc::NetworkConstraints::multiple_connections)
        .def_readonly("min_attached_ends", &bdmc::NetworkConstraints::min_attached_ends)
        .def_readonly("max_attached_ends", &bdmc::NetworkConstraints::max_attached_ends);
    py::class_<bdmc::NetworkKinetics>(m, "Kinetics")
        .def(py::init([](const std::string& law, double lcap, double barrier, double delta) {
            bdmc::NetworkKinetics k; k.law=transition_law(law); k.l_cap=lcap; k.energy_barrier=barrier; k.delta_t0=delta; k.validate(); return k;
        }), py::arg("law")="uniform", py::arg("l_cap")=0., py::arg("energy_barrier")=0., py::arg("delta_t0")=0.)
        .def_readonly("l_cap", &bdmc::NetworkKinetics::l_cap)
        .def_readonly("energy_barrier", &bdmc::NetworkKinetics::energy_barrier)
        .def_readonly("delta_t0", &bdmc::NetworkKinetics::delta_t0);
    py::class_<bdmc::NetworkSchedule>(m, "Schedule")
        .def(py::init([](double dt, std::uint64_t steps, double C, double Rt) {
            bdmc::NetworkSchedule s; s.dt_over_tauR=dt; s.steps_per_block=steps; s.repulsion_coefficient=C; s.Rt=Rt; s.validate(); return s;
        }), py::arg("dt_over_tauR"), py::arg("steps_per_block"), py::arg("repulsion_coefficient"), py::arg("Rt"))
        .def_property_readonly("delta_t0", &bdmc::NetworkSchedule::delta_t0)
        .def("attempts_per_block", &bdmc::NetworkSchedule::attempts_per_block, py::arg("topology"));
    py::class_<bdmc::NetworkProposal>(m, "Proposal")
        .def(py::init([](std::size_t anchor, std::size_t other, std::size_t destination,
                        std::size_t end, std::uint64_t epoch) {
            bdmc::NetworkProposal p; p.anchor=anchor; p.other=other; p.destination=destination;
            p.end=end; p.topology_epoch=epoch; return p;
        }), py::arg("anchor"), py::arg("other"), py::arg("destination"), py::arg("end"), py::arg("topology_epoch"),
            "Explicit manual choices; no evaluated action/proposal mass; attempt_transition recomputes the decision")
        .def_readonly("anchor", &bdmc::NetworkProposal::anchor).def_readonly("other", &bdmc::NetworkProposal::other)
        .def_readonly("destination", &bdmc::NetworkProposal::destination).def_readonly("end", &bdmc::NetworkProposal::end)
        .def_readonly("topology_epoch", &bdmc::NetworkProposal::topology_epoch)
        .def_property_readonly("action", [](const bdmc::NetworkProposal& p) { return bdmc::network_action_name(p.action); })
        .def_property_readonly("rejection", [](const bdmc::NetworkProposal& p) { return bdmc::network_rejection_name(p.rejection); })
        .def_readonly("proposal_probability", &bdmc::NetworkProposal::proposal_probability)
        .def_readonly("acceptance_probability", &bdmc::NetworkProposal::acceptance_probability)
        .def_readonly("mutation_probability", &bdmc::NetworkProposal::mutation_probability);
    m.def("metadata", &metadata);
    m.def("transition_probability", &bdmc::network_transition_probability, py::arg("old_tension"), py::arg("kinetics"));
    m.def("detachment_weight", [](double tension, double lcap, const std::string& law) {
        return bdmc::network_detachment_weight(tension,lcap,detachment_law(law));
    }, py::arg("old_tension"), py::arg("l_cap"), py::arg("law")="uniform");
    m.def("discrete_distribution", [](std::size_t n, const MapRow& raw) {
        return distribution(bdmc::network_discrete_distribution(n,weights(raw)));
    }, py::arg("n_particles"), py::arg("raw_weights"));
    m.def("suggestion_distribution", [](const bdmc::NetworkTopology& t, std::size_t a, const MapRow& raw) {
        return distribution(bdmc::network_suggestion_distribution(t,a,weights(raw)));
    }, py::arg("topology"), py::arg("anchor"), py::arg("detachment_weights"));
    m.def("classify_action", [](const bdmc::NetworkTopology& t, std::size_t a, std::size_t o, std::size_t d, const bdmc::NetworkConstraints& c) {
        const auto r=bdmc::classify_network_action(t,a,o,d,c); py::dict out;
        out["action"]=bdmc::network_action_name(r.action); out["rejection"]=bdmc::network_rejection_name(r.rejection); return out;
    }, py::arg("topology"), py::arg("anchor"), py::arg("other"), py::arg("destination"), py::arg("constraints")=bdmc::NetworkConstraints{});
    m.def("enumerate_proposals", [](const bdmc::NetworkTopology& t, const MapRows& a, const MapRows& d,
        const MapRows& tensions, const bdmc::NetworkConstraints& c, const bdmc::NetworkKinetics& k) {
        return bdmc::enumerate_network_proposals(t,weights(a),weights(d),weights(tensions),c,k);
    }, py::arg("topology"), py::arg("association"), py::arg("detachment"), py::arg("old_tensions"),
        py::arg("constraints")=bdmc::NetworkConstraints{}, py::arg("kinetics")=bdmc::NetworkKinetics{});
    m.def("attempt_transition", [](bdmc::NetworkTopology& t, const bdmc::NetworkProposal& p, double tension,
        double roll, const bdmc::NetworkConstraints& c, const bdmc::NetworkKinetics& k) {
        const auto r=bdmc::attempt_network_transition(t,p,tension,roll,c,k); py::dict out;
        out["action"]=bdmc::network_action_name(r.action); out["rejection"]=bdmc::network_rejection_name(r.rejection);
        out["acceptance_probability"]=r.acceptance_probability; return out;
    }, py::arg("topology"), py::arg("proposal"), py::arg("old_tension"), py::arg("acceptance_roll"),
        py::arg("constraints")=bdmc::NetworkConstraints{}, py::arg("kinetics")=bdmc::NetworkKinetics{});
    py::class_<bdmc::NetworkConnectorSpec>(m, "ConnectorSpec")
        .def(py::init([](const std::string& law, double h, double q0, double scale) {
            bdmc::NetworkConnectorSpec s;
            if (law=="hookean") s.law=bdmc::NetworkConnectorLaw::Hookean;
            else if (law=="fene") s.law=bdmc::NetworkConnectorLaw::FENE;
            else throw std::invalid_argument("unsupported network connector law");
            s.spring_constant=h; s.q0=q0; s.scale=scale; s.validate(); return s;
        }), py::arg("law")="hookean", py::arg("spring_constant")=3., py::arg("q0")=10., py::arg("scale")=1.)
        .def_readonly("spring_constant", &bdmc::NetworkConnectorSpec::spring_constant)
        .def_readonly("q0", &bdmc::NetworkConnectorSpec::q0).def_readonly("scale", &bdmc::NetworkConnectorSpec::scale);
    m.def("connector_forces", [](const bdmc::NetworkTopology& t, py::object input, const std::vector<double>& lengths,
        const bdmc::NetworkConnectorSpec& spec) {
        const auto f=bdmc::evaluate_network_connectors(t,positions(input,t.n_particles()),lengths,spec); py::dict out;
        out["forces"]=copied(f.forces,{size(t.n_particles()),3}); out["potential_energy"]=f.potential_energy;
        std::vector<double> virial(f.connector_virial.begin(),f.connector_virial.end());
        out["connector_virial"]=copied(virial,{3,3});
        out["convention"]="positive connector tension q outer F_tension; once per physical bridge; no volume normalization; scale applies to force, energy and virial";
        return out;
    }, py::arg("topology"), py::arg("positions"), py::arg("box_lengths"), py::arg("spec")=bdmc::NetworkConnectorSpec{});
}
