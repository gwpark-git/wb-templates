"""Labeled network topology, finite legacy proposal measure and bonded geometry."""

from itertools import product

import numpy as np
import pytest

from bdmc import _network as net


def proposal(topology, anchor, other, destination, end):
    return net.Proposal(anchor, other, destination, end, topology.epoch)


def dense_snapshot(topology):
    snapshot = topology.snapshot()
    dense = np.zeros((topology.n_particles, topology.n_particles), dtype=np.uint64)
    for i, (start, stop) in enumerate(zip(snapshot["adjacency_offsets"][:-1], snapshot["adjacency_offsets"][1:])):
        dense[i, snapshot["adjacency_particles"][start:stop]] = snapshot["adjacency_weights"][start:stop]
    return snapshot, dense


def check_independent_counts(topology):
    snapshot, dense = dense_snapshot(topology)
    expected = np.zeros_like(dense)
    for a, b in snapshot["endpoints"]:
        expected[a, b] += 1
        expected[b, a] += 1
    np.testing.assert_array_equal(dense, expected)
    np.testing.assert_array_equal(snapshot["attached_ends"], np.bincount(snapshot["endpoints"].ravel().astype(np.int64), minlength=topology.n_particles))
    assert snapshot["attached_ends"].sum() == 2 * topology.n_chains
    assert snapshot["n_bridges"] == np.count_nonzero(snapshot["endpoints"][:, 0] != snapshot["endpoints"][:, 1])
    assert snapshot["n_distinct_pairs"] == np.count_nonzero(np.triu(dense, 1))
    topology.validate()


def test_owned_endpoint_state_and_csr_counts():
    initial = np.array([[0, 0], [0, 2], [2, 0], [1, 1]])
    topology = net.Topology(3, initial)
    initial[:] = 1
    snapshot, dense = dense_snapshot(topology)
    np.testing.assert_array_equal(snapshot["endpoints"], [[0, 0], [0, 2], [2, 0], [1, 1]])
    np.testing.assert_array_equal(dense, [[2, 0, 2], [0, 2, 0], [2, 0, 0]])
    assert snapshot["n_bridges"] == 2
    assert snapshot["n_distinct_pairs"] == 1
    for key in ("endpoints", "attached_ends", "adjacency_particles", "adjacency_weights", "adjacency_offsets"):
        assert snapshot[key].dtype == np.uint64
        assert snapshot[key].flags.owndata and snapshot[key].flags.c_contiguous
        snapshot[key][:] = 99
    check_independent_counts(topology)
    assert topology.particle_of_end(1) == 0
    assert topology.particle_of_end(5) == 2
    assert topology.opposite_end(1) == 5
    np.testing.assert_array_equal(topology.compatible_ends(0, 0), [0, 4])
    np.testing.assert_array_equal(topology.compatible_ends(0, 2), [2, 5])
    topology.compatible_ends(0, 2)[:] = 99
    np.testing.assert_array_equal(topology.compatible_ends(0, 2), [2, 5])


def test_atomic_batch_rollback_and_identity():
    topology = net.Topology.loops(3, 2)
    topology.transact([(0, 0, 2), (1, 0, 2)])
    before = topology.snapshot()
    with pytest.raises(ValueError, match="expectation"):
        topology.transact([(0, 2, 1), (1, 1, 0)])
    for key in ("endpoints", "attached_ends", "adjacency_weights"):
        np.testing.assert_array_equal(topology.snapshot()[key], before[key])
    assert topology.epoch == before["epoch"]
    with pytest.raises(ValueError, match="end ID"):
        topology.transact([(0, 2, 1), (99, 0, 1)])
    assert topology.epoch == before["epoch"]
    topology.transact([(0, 2, 1), (0, 1, 2)])
    assert topology.epoch == before["epoch"]
    topology.transact([(0, 2, 0), (1, 2, 0)])
    np.testing.assert_array_equal(topology.snapshot()["endpoints"], net.Topology.loops(3, 2).snapshot()["endpoints"])
    check_independent_counts(topology)


def test_serial_four_actions_constraints_rejection_and_staleness():
    topology = net.Topology(3, [[0, 0], [0, 0]])
    cases = [(0, 0, 1, 0, "add"), (0, 0, 1, 1, "add"),
             (0, 1, 2, 0, "move"), (0, 1, 1, 1, "cancel"),
             (0, 1, 0, 1, "opposite_delete"), (0, 2, 0, 0, "opposite_delete")]
    for anchor, other, destination, end, action in cases:
        result = net.attempt_transition(topology, proposal(topology, anchor, other, destination, end), 0, 0)
        assert result["action"] == action
        check_independent_counts(topology)
    fixed = net.Constraints(min_attached_ends=4, max_attached_ends=4)
    assert net.classify_action(topology, 0, 0, 1, fixed)["rejection"] == "count_or_multiplicity"
    p = proposal(topology, 0, 0, 1, 0)
    kinetics = net.Kinetics("dissociation", delta_t0=0.25)
    assert net.attempt_transition(topology, p, 0, 0.25, kinetics=kinetics)["rejection"] == "acceptance"
    assert topology.epoch == p.topology_epoch
    assert net.attempt_transition(topology, p, 0, 0.249, kinetics=kinetics)["action"] == "add"
    with pytest.raises(ValueError, match="stale"):
        net.attempt_transition(topology, p, 0, 0)
    single = net.Constraints(multiple_connections=False, min_attached_ends=100, max_attached_ends=100)
    assert net.classify_action(topology, 0, 0, 1, single)["rejection"] == "count_or_multiplicity"
    assert net.classify_action(topology, 0, 1, 0, single)["action"] == "opposite_delete"


@pytest.mark.parametrize("law", ["uniform", "metropolis", "weighted", "dissociation", "first_order"])
def test_mapped_kinetics_and_finite_small_hazard(law):
    for tension in (-2, 0, 0.7, 9):
        for delta in (0, 0.001, 0.1, 2):
            k = net.Kinetics(law, l_cap=0.4, energy_barrier=2, delta_t0=delta)
            rate = np.exp(tension * 0.4)
            expected = {"uniform": 1, "metropolis": min(1, np.exp(tension * 0.4 - 2)),
                        "weighted": min(1, rate), "dissociation": min(1, delta * rate),
                        "first_order": -np.expm1(-delta * rate)}[law]
            assert net.transition_probability(tension, k) == pytest.approx(expected, abs=2e-15)
    assert net.transition_probability(1, net.Kinetics("first_order", delta_t0=1e-20)) == pytest.approx(1e-20, abs=1e-35)
    assert net.detachment_weight(9, 0.4, "weighted") == pytest.approx(np.exp(3.6))
    assert net.detachment_weight(9, 0.4) == 1


def test_finite_discrete_measure_and_schedule():
    candidate = net.discrete_distribution(3, [(0, 1), (2, 3), (1, 0)])
    np.testing.assert_array_equal(candidate["particles"], [0, 2, 1])
    np.testing.assert_allclose(candidate["probabilities"], [0.25, 0.75, 0])
    topology = net.Topology(3, [[0, 0], [0, 1], [1, 0], [0, 2]])
    old = net.suggestion_distribution(topology, 0, [(0, 2), (1, 3), (2, 5)])
    np.testing.assert_array_equal(old["particles"], [0, 1, 2])
    np.testing.assert_allclose(old["probabilities"], np.array([4, 6, 5]) / 15)
    schedule = net.Schedule(0.003, 7, 2, 5)
    assert schedule.delta_t0 == pytest.approx(0.0021)
    assert schedule.attempts_per_block(topology) == 4


@pytest.mark.parametrize("n,chains", [(2, 1), (2, 2), (3, 1), (3, 2)])
def test_exhaustive_labeled_end_proposal_measure(n, chains):
    association = [[(j, float(j + 1)) for j in range(n)] for i in range(n)]
    detachment = [[(j, float(i + j + 1)) for j in range(n)] for i in range(n)]
    tensions = [[(j, 0.2 * (i + j)) for j in range(n)] for i in range(n)]
    for state in product(range(n), repeat=2 * chains):
        endpoints = np.array(state).reshape(chains, 2)
        topology = net.Topology(n, endpoints)
        actual = net.enumerate_proposals(topology, association, detachment, tensions)
        counts = np.bincount(endpoints.ravel(), minlength=n)
        occupied = np.count_nonzero(counts)
        expected = {}
        for anchor in range(n):
            if counts[anchor] == 0:
                continue
            denominator = sum(anchor + endpoints[c, 1 - side] + 1
                              for c in range(chains) for side in range(2) if endpoints[c, side] == anchor)
            for c in range(chains):
                for fixed_side in range(2):
                    if endpoints[c, fixed_side] != anchor:
                        continue
                    other = endpoints[c, 1 - fixed_side]
                    end = c + chains * (1 - fixed_side)
                    for destination in range(n):
                        expected[(anchor, other, destination, end)] = ((anchor + other + 1) / denominator /
                            occupied * (destination + 1) / sum(range(1, n + 1)))
        assert len(actual) == len(expected)
        assert sum(p.proposal_probability for p in actual) == pytest.approx(1)
        for p in actual:
            assert p.proposal_probability == pytest.approx(expected[(p.anchor, p.other, p.destination, p.end)])
            assert p.mutation_probability == pytest.approx(p.proposal_probability if p.other != p.destination else 0)


def test_constraint_and_kinetic_rejected_mass_is_retained():
    topology = net.Topology(2, [[0, 0], [0, 1]])
    maps = [[(0, 1), (1, 1)], [(0, 1), (1, 1)]]
    tensions = [[(0, 0), (1, 0)], [(0, 0), (1, 0)]]
    proposals = net.enumerate_proposals(topology, maps, maps, tensions,
        net.Constraints(multiple_connections=False), net.Kinetics("dissociation", delta_t0=0.3))
    assert sum(p.proposal_probability for p in proposals) == pytest.approx(1)
    assert any(p.rejection == "count_or_multiplicity" for p in proposals)
    assert any(p.rejection == "same_destination" for p in proposals)
    accepted_mass = sum(p.mutation_probability for p in proposals)
    assert 0 < accepted_mass < 0.3
    assert sum(p.proposal_probability * (1 - p.acceptance_probability) for p in proposals) == pytest.approx(0.7)


def test_proposal_choices_and_cached_mass_cannot_diverge_by_python_mutation():
    topology = net.Topology(2, [[0, 0]])
    maps = [[(0, 1), (1, 1)], [(0, 1), (1, 1)]]
    tensions = [[(0, 0), (1, 0)], [(0, 0), (1, 0)]]
    generated = next(p for p in net.enumerate_proposals(topology, maps, maps, tensions) if p.action == "add")
    manual = net.Proposal(0, 0, 1, 0, topology.epoch)
    for p in (generated, manual):
        before = (p.anchor, p.other, p.destination, p.end, p.topology_epoch, p.proposal_probability)
        for field in ("anchor", "other", "destination", "end", "topology_epoch", "proposal_probability"):
            with pytest.raises(AttributeError):
                setattr(p, field, 99)
        assert (p.anchor, p.other, p.destination, p.end, p.topology_epoch, p.proposal_probability) == before
    assert manual.proposal_probability == 0
    assert net.attempt_transition(topology, manual, 0, 0)["action"] == "add"
    np.testing.assert_array_equal(topology.snapshot()["endpoints"], [[1, 0]])


@pytest.mark.parametrize("law", ["hookean", "fene"])
def test_bonded_force_numpy_chain_oracle_and_copy(law):
    endpoints = np.array([[0, 3], [3, 0], [1, 1], [1, 2]])
    topology = net.Topology(4, endpoints)
    x = np.array([[0.1, 0.2, 0.3], [1, 2, 3], [1.3, 2.4, 3.2], [6.1, 0.2, 0.3]])
    box, h, q0, scale = np.array([20, 20, 20]), 2.4, 8, 0.25
    forces, virial, energy = np.zeros_like(x), np.zeros((3, 3)), 0.0
    for a, b in endpoints:
        if a == b:
            continue
        q = x[b] - x[a]
        q -= box * np.floor(q / box + 0.5)
        r2 = q @ q
        tension = scale * h * q / (1 - r2 / q0**2 if law == "fene" else 1)
        forces[a] += tension
        forces[b] -= tension
        virial += np.outer(q, tension)
        energy += -scale * h * q0**2 / 2 * np.log1p(-r2 / q0**2) if law == "fene" else scale * h * r2 / 2
    result = net.connector_forces(topology, x, box, net.ConnectorSpec(law, h, q0, scale))
    np.testing.assert_allclose(result["forces"], forces, atol=2e-14)
    np.testing.assert_allclose(result["connector_virial"], virial, atol=2e-14)
    assert result["potential_energy"] == pytest.approx(energy)
    np.testing.assert_allclose(result["forces"].sum(axis=0), 0, atol=2e-15)
    for key in ("forces", "connector_virial"):
        assert result[key].flags.owndata and result[key].flags.c_contiguous
        result[key][:] = 99
    np.testing.assert_array_equal(topology.snapshot()["endpoints"], endpoints)
    pair = net.Topology(2, [[0, 1]])
    periodic = net.connector_forces(pair, [[9.8, 0, 0], [0.3, 0, 0]], [10, 10, 10])
    np.testing.assert_allclose(periodic["forces"], [[1.5, 0, 0], [-1.5, 0, 0]], atol=1e-14)


def test_invalid_domains_fail_without_state_mutation():
    for endpoints in ([[0, 2]], [[0, -1]], [[0.5, 1]], [[0, 1, 2]], [0, 1]):
        with pytest.raises((ValueError, TypeError)):
            net.Topology(2, endpoints)
    topology = net.Topology(2, np.empty((0, 2), dtype=np.int64))
    assert topology.snapshot()["endpoints"].shape == (0, 2)
    with pytest.raises(ValueError):
        topology.particle_of_end(0)
    with pytest.raises(ValueError, match="no nonempty"):
        net.enumerate_proposals(topology, [[], []], [[], []], [[], []])
    for raw in ([(0, 0), (1, 0)], [(0, 1), (0, 1)], [(0, -1)], [(2, 1)], [(0, np.inf)]):
        with pytest.raises(ValueError):
            net.discrete_distribution(2, raw)
    pair = net.Topology(2, [[0, 1]])
    with pytest.raises(ValueError, match="r<q0"):
        net.connector_forces(pair, [[0, 0, 0], [0.5, 0, 0]], [10, 10, 10], net.ConnectorSpec("fene", q0=0.5))
    with pytest.raises(ValueError, match="shape"):
        net.connector_forces(pair, [0, 0, 0, 1, 0, 0], [10, 10, 10])
    np.testing.assert_array_equal(pair.snapshot()["endpoints"], [[0, 1]])
    with pytest.raises(ValueError):
        net.Kinetics("glauber")
    with pytest.raises(ValueError):
        net.Kinetics("dissociation", delta_t0=-1)
    for law in ("dissociation", "first_order"):
        with pytest.raises(ValueError, match="kinetic rate exponential"):
            net.transition_probability(720, net.Kinetics(law, l_cap=1, delta_t0=1e-320))
        with pytest.raises(ValueError, match="kinetic rate exponential"):
            net.transition_probability(-750, net.Kinetics(law, l_cap=1, delta_t0=1e308))
    with pytest.raises(ValueError):
        net.Schedule(0.001, 0, 1, 1)
    assert "full MC/BD controller" in net.metadata()["unsupported"]
