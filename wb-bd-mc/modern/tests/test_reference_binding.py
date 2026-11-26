"""Owned NumPy interface and physical conventions of the additive references."""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from bdmc import _reference as ref


def test_repulsive_sampling_ownership_and_initialization():
    spec = ref.RepulsiveSpec(n_particles=2, diffusion=0)
    initial = np.array([[-0.2, 0, 0], [0.3, 0, 0]])
    sim = ref.RepulsiveSimulation(spec, initial_positions=initial)
    initial[:] = 99
    spec.dt = 0.7
    before = sim.observe()
    np.testing.assert_allclose(before["positions"], [[9.8, 0, 0], [0.3, 0, 0]])
    np.testing.assert_allclose(before["unwrapped_positions"], [[-0.2, 0, 0], [0.3, 0, 0]])
    assert sim.spec.dt == 0.001
    before["positions"][:] = -9
    assert np.min(sim.observe()["positions"]) >= 0
    batch = sim.sample(7, stride=3)
    np.testing.assert_array_equal(batch["steps"], [0, 3, 6, 7])
    np.testing.assert_allclose(batch["times"], np.array([0, 3, 6, 7]) * 0.001)
    for key in ("positions", "unwrapped_positions", "forces"):
        assert batch[key].shape == (4, 2, 3)
        assert batch[key].dtype == np.float64
        assert batch[key].flags.c_contiguous
        assert batch[key].flags.owndata
    assert batch["pair_virial"].shape == (4, 3, 3)
    assert batch["potential_energy"].shape == (4,)
    assert batch["steps"].dtype == np.uint64
    assert sim.sample(0, include_initial=False)["positions"].shape == (0, 2, 3)
    batch["unwrapped_positions"][:] = 100
    assert np.max(sim.observe()["unwrapped_positions"]) < 100


def test_repulsive_numpy_allpairs_oracle():
    rng = np.random.default_rng(932)
    n, box, sigma, epsilon = 70, 5.0, 1.3, 2.1
    x = rng.uniform(-box, 2 * box, (n, 3))
    expected_force = np.zeros_like(x)
    expected_virial = np.zeros((3, 3))
    expected_energy = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            q = x[j] - x[i]
            q -= box * np.floor(q / box + 0.5)
            r = np.linalg.norm(q)
            if r < sigma:
                f = epsilon / sigma * (1 - (r / sigma) ** 2) * q / r
                expected_force[i] -= f
                expected_force[j] += f
                expected_energy += epsilon / 3 * (1 - r / sigma) ** 2 * (2 + r / sigma)
                expected_virial += np.outer(q, f)
    for search in ("cell_list", "brute_force"):
        s = ref.RepulsiveSpec(n_particles=n, box_length=box, sigma=sigma, epsilon=epsilon, force_search=search)
        f = ref.RepulsiveSimulation(s, x).observe()
        np.testing.assert_allclose(f["forces"], expected_force, atol=2e-14)
        np.testing.assert_allclose(f["pair_virial"], expected_virial, atol=5e-14)
        assert f["potential_energy"] == pytest.approx(expected_energy, abs=2e-14)
        np.testing.assert_allclose(f["forces"].sum(axis=0), 0, atol=1e-14)


@pytest.mark.parametrize("model", ["repulsive", "dumbbell"])
def test_batch_boundaries_and_thread_serialization(model):
    if model == "repulsive":
        spec = ref.RepulsiveSpec(n_particles=13, seed=42)
        factory, key = ref.RepulsiveSimulation, "unwrapped_positions"
    else:
        spec = ref.DumbbellSpec(n_dumbbells=13, seed=42)
        factory, key = ref.DumbbellSimulation, "connectors"
    whole, split, shared = [factory(spec) for _ in range(3)]
    whole.advance(21)
    split.sample(8, stride=3)
    split.advance(13)
    np.testing.assert_array_equal(whole.observe()[key], split.observe()[key])
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(shared.advance, [10, 11]))
    np.testing.assert_array_equal(whole.observe()[key], shared.observe()[key])


def test_dumbbell_observables_step_shear_and_owned_arrays():
    spec = ref.DumbbellSpec(n_dumbbells=2, spring_constant=2, mobility=0.5,
                            diffusion=0.7, protocol="step_shear", step_strain=-0.5)
    initial = np.array([[0.3, 0.8, 0.1], [-0.3, -0.8, -0.1]])
    expected = initial.copy()
    expected[:, 0] -= 0.5 * expected[:, 1]
    sim = ref.DumbbellSimulation(spec, initial_connectors=initial)
    initial[:] = 8
    f = sim.observe()
    np.testing.assert_allclose(f["connectors"], expected)
    np.testing.assert_allclose(f["spring_forces"], 2 * expected)
    moment = expected.T @ expected / 2
    np.testing.assert_allclose(f["conformation"], moment)
    np.testing.assert_allclose(f["polymer_stress"], 2 * moment - 1.4 * np.eye(3))
    assert f["mean_energy"] == pytest.approx(np.mean(np.sum(expected**2, axis=1)))
    assert f["time"] == f["step"] == f["shear_rate"] == 0
    f["connectors"][:] = 88
    np.testing.assert_allclose(sim.observe()["connectors"], expected)
    batch = sim.sample(4, stride=3)
    np.testing.assert_array_equal(batch["steps"], [0, 3, 4])
    assert batch["connectors"].shape == batch["spring_forces"].shape == (3, 2, 3)
    assert batch["conformation"].shape == batch["polymer_stress"].shape == (3, 3, 3)
    assert batch["mean_energy"].shape == batch["shear_rates"].shape == (3,)
    for key in ("connectors", "spring_forces", "conformation", "polymer_stress"):
        assert batch[key].flags.owndata and batch[key].flags.c_contiguous
    assert sim.sample(0, include_initial=False)["connectors"].shape == (0, 2, 3)


@pytest.mark.parametrize("rate", [-120.0, -1.7, 1.7, 120.0])
def test_signed_hookean_shear_without_periodic_box(rate):
    s = ref.DumbbellSpec(n_dumbbells=1, diffusion=0, dt=0.01,
                        protocol="simple_shear", shear_rate=rate)
    q = np.array([[0.2, 0.7, -0.1]])
    sim = ref.DumbbellSimulation(s, q)
    sim.advance(1)
    expected = 0.98 * q
    expected[0, 0] += rate * 0.01 * q[0, 1]
    np.testing.assert_allclose(sim.observe()["connectors"], expected, atol=1e-15)
    assert sim.observe()["shear_rate"] == rate
    assert "unbounded connector" in sim.metadata["boundary"]


def test_fene_force_energy_domain_and_transaction():
    s = ref.DumbbellSpec(n_dumbbells=1, spring_law="fene", q0=2,
                        spring_constant=1.8, diffusion=0, dt=0.001)
    q = np.array([0.3, -0.4, 0.2])
    ratio = q @ q / 4
    np.testing.assert_allclose(ref.dumbbell_spring_force(q, s), 1.8 * q / (1 - ratio))
    assert ref.dumbbell_spring_energy(q, s) == pytest.approx(-1.8 * 4 / 2 * np.log1p(-ratio))
    with pytest.raises(ValueError, match="FENE"):
        ref.DumbbellSimulation(s, [[2.0, 0, 0]])
    s.protocol, s.step_strain = "step_shear", 5
    with pytest.raises(ValueError, match="FENE"):
        ref.DumbbellSimulation(s, [[0, 0.5, 0]])
    s.protocol, s.step_strain, s.q0, s.dt = "no_flow", 0, 1, 0.1
    sim = ref.DumbbellSimulation(s, [[0.99, 0, 0]])
    before = sim.observe()
    for _ in range(2):
        with pytest.raises(ValueError, match="FENE"):
            sim.advance(1)
        after = sim.observe()
        np.testing.assert_array_equal(after["connectors"], before["connectors"])
        assert after["step"] == after["time"] == 0


@pytest.mark.parametrize("constructor,kwargs", [
    (ref.RepulsiveSpec, {"n_particles": 0}),
    (ref.RepulsiveSpec, {"sigma": 6}),
    (ref.RepulsiveSpec, {"epsilon": -1}),
    (ref.RepulsiveSpec, {"force_search": "unknown"}),
    (ref.RepulsiveSpec, {"mobility": 0}),
    (ref.RepulsiveSpec, {"diffusion": -0.1}),
    (ref.DumbbellSpec, {"n_dumbbells": 0}),
    (ref.DumbbellSpec, {"dt": 1}),
    (ref.DumbbellSpec, {"spring_constant": 0}),
    (ref.DumbbellSpec, {"spring_law": "unknown"}),
    (ref.DumbbellSpec, {"protocol": "periodic_shear"}),
    (ref.DumbbellSpec, {"q0": 0}),
    (ref.DumbbellSpec, {"noise_policy": "unknown"}),
    (ref.DumbbellSpec, {"shear_rate": np.inf}),
])
def test_invalid_specs(constructor, kwargs):
    with pytest.raises((ValueError, TypeError)):
        constructor(**kwargs)


def test_invalid_arrays_coincidence_mutations_and_sampling():
    r = ref.RepulsiveSpec(n_particles=2)
    for x in (np.zeros(6), np.zeros((2, 2)), [[0, 0, 0], [np.nan, 0, 0]]):
        with pytest.raises(ValueError):
            ref.RepulsiveSimulation(r, x)
    with pytest.raises(ValueError, match="coincident"):
        ref.RepulsiveSimulation(r, [[0, 0, 0], [10, 0, 0]])
    r.dt = -1
    with pytest.raises(ValueError):
        ref.RepulsiveSimulation(r)
    d = ref.DumbbellSpec(n_dumbbells=2)
    with pytest.raises(ValueError):
        ref.DumbbellSimulation(d, [[0, 0, 0]])
    sim = ref.DumbbellSimulation(d)
    with pytest.raises(ValueError):
        sim.sample(1, stride=0)
    with pytest.raises(TypeError):
        sim.advance(-1)


def test_declared_metadata_and_restart_scope():
    r = ref.RepulsiveSimulation(ref.RepulsiveSpec(diffusion=0.7, mobility=0.35))
    assert r.spec.to_dict()["noise_policy"] == "gaussian"
    assert r.spec.dimension == 3
    assert r.metadata["thermal_energy"] == pytest.approx(2)
    assert "epsilon/3" in r.metadata["potential"]
    assert "explicit mapping" in r.metadata["time_unit"]
    assert "regular grid" in r.metadata["initialization"]
    d = ref.DumbbellSimulation(ref.DumbbellSpec(noise_policy="legacy_uniform"))
    assert d.spec.to_dict()["spring_law"] == "hookean"
    assert "not Gaussian" in d.metadata["noise_finite_step_law"]
    assert "positive-Cauchy" in d.metadata["polymer_stress_convention"]
    assert "not domain-safe" in d.metadata["fene_domain_policy"]
    for sim in (r, d):
        assert sim.metadata["schema"] == "bdmc.reference.observation.v1"
        assert "unsupported" in sim.metadata["restart_contract"]
        assert not hasattr(sim, "checkpoint")
        assert not hasattr(sim, "from_checkpoint")
