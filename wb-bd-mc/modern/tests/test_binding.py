"""Direct NumPy/native interface contracts, independent of the ML package."""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from bdmc import _core


def test_batched_arrays_times_and_empty_batch():
    spec = _core.SimulationSpec(n_particles=4, dimension=2, dt=0.02, seed=45)
    simulation = _core.Simulation(spec)
    batch = simulation.sample(7, stride=3)
    assert batch["positions"].shape == (4, 4, 2)
    assert batch["forces"].shape == batch["positions"].shape
    assert batch["unwrapped_positions"].shape == batch["positions"].shape
    np.testing.assert_array_equal(batch["steps"], [0, 3, 6, 7])
    np.testing.assert_allclose(batch["times"], np.array([0, 3, 6, 7]) * 0.02)
    assert batch["steps"].dtype == np.dtype("uint64")
    assert batch["positions"].dtype == np.dtype("float64")
    assert batch["positions"].flags.c_contiguous
    assert np.all(batch["forces"] == 0)
    assert np.all(batch["potential_energy"] == 0)
    assert batch["fixed_bonds"].shape == (0, 2)
    empty = simulation.sample(0, include_initial=False)
    assert empty["positions"].shape == (0, 4, 2)
    assert empty["particle_ids"].shape == (4,)
    subsequent = simulation.sample(2, stride=5, include_initial=False)
    np.testing.assert_array_equal(subsequent["steps"], [9])


def test_array_ownership_and_specification_snapshot():
    spec = _core.SimulationSpec(n_particles=2, dimension=1, diffusivity=0)
    original = np.array([[-0.2], [10.3]])
    simulation = _core.Simulation(spec, initial_positions=original)
    original[:] = 77
    spec.dt = 4
    snapshot = simulation.observe()
    np.testing.assert_allclose(snapshot["positions"], [[9.8], [0.3]])
    np.testing.assert_allclose(snapshot["unwrapped_positions"], [[-0.2], [10.3]])
    assert simulation.spec.dt == 0.001
    snapshot["positions"][:] = -99
    assert np.all(simulation.observe()["positions"] >= 0)
    batch = simulation.sample(1)
    batch["unwrapped_positions"][:] = 42
    np.testing.assert_allclose(simulation.observe()["unwrapped_positions"], [[-0.2], [10.3]])
    assert simulation.metadata["initialization"] == "explicitly supplied positions"


def test_seed_stream_and_concurrent_advance_serialization():
    spec = _core.SimulationSpec(n_particles=100, dimension=3, seed=7)
    reference, split, shared = (_core.Simulation(spec) for _ in range(3))
    reference.advance(60)
    split.sample(23, stride=5)
    split.advance(37)
    np.testing.assert_array_equal(
        reference.observe()["unwrapped_positions"], split.observe()["unwrapped_positions"]
    )
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(shared.advance, [30, 30]))
    np.testing.assert_array_equal(
        reference.observe()["unwrapped_positions"], shared.observe()["unwrapped_positions"]
    )
    assert shared.observe()["step"] == 60


def test_dumbbell_exact_force_and_simultaneous_update():
    spec = _core.SimulationSpec(
        model="dumbbell", n_particles=2, dimension=1, box_lengths=[10],
        spring_alpha=2, mobility=0.75, diffusivity=0, dt=0.01,
    )
    simulation = _core.Simulation(
        spec, initial_positions=np.array([[9.8], [0.2]]),
        particle_ids=np.array([100, 7], dtype=np.uint64),
    )
    before = simulation.observe()
    np.testing.assert_array_equal(before["fixed_bonds"], [[7, 100]])
    np.testing.assert_allclose(before["forces"], [[1.6], [-1.6]], atol=1e-13)
    assert before["potential_energy"] == pytest.approx(0.32)
    simulation.advance(1)
    np.testing.assert_allclose(simulation.observe()["positions"], [[9.812], [0.188]])
    assert simulation.metadata["spring_coefficient"] == 4


def test_uniform_noise_contract_and_declared_gaussian_alternative():
    dt, diffusion = 0.02, 0.7
    common = dict(n_particles=40000, dimension=1, dt=dt, diffusivity=diffusion, seed=234)
    uniform = _core.Simulation(_core.SimulationSpec(**common)).sample(1, include_initial=False)
    x = uniform["unwrapped_positions"].ravel()
    variance = 2 * diffusion * dt
    assert np.max(np.abs(x)) <= np.sqrt(3 * variance)
    assert np.mean(x*x) == pytest.approx(variance, rel=0.025)
    assert np.mean(x**4) / variance**2 == pytest.approx(1.8, abs=0.07)
    gaussian = _core.Simulation(
        _core.SimulationSpec(**common, noise_policy="gaussian")
    ).sample(1, include_initial=False)
    y = gaussian["unwrapped_positions"].ravel()
    assert np.mean(y*y) == pytest.approx(variance, rel=0.035)
    assert np.mean(y**4) / variance**2 == pytest.approx(3, abs=0.18)
    assert "uniform" in uniform["metadata"]["noise_finite_step_law"]
    assert "Gaussian" in gaussian["metadata"]["noise_finite_step_law"]


@pytest.mark.parametrize("kwargs", [
    {"dt": 0}, {"diffusivity": -1}, {"mobility": float("nan")},
    {"dimension": 4}, {"dimension": 10**12}, {"n_particles": 0},
    {"box_lengths": [0, 1, 1]}, {"model": "transient_network"},
    {"noise_policy": "arbitrary"}, {"model": "dumbbell", "n_particles": 3},
    {"model": "dumbbell", "dt": 1},
])
def test_invalid_physical_specifications_rejected(kwargs):
    with pytest.raises((ValueError, TypeError)):
        _core.SimulationSpec(**kwargs)


def test_invalid_inputs_and_mutated_specification_rejected():
    spec = _core.SimulationSpec(n_particles=2, dimension=1)
    with pytest.raises(ValueError):
        _core.Simulation(spec, initial_positions=np.zeros(2))
    with pytest.raises(ValueError):
        _core.Simulation(spec, initial_positions=np.array([[0], [np.nan]]))
    with pytest.raises(ValueError):
        _core.Simulation(spec, particle_ids=[4, 4])
    with pytest.raises((TypeError, ValueError)):
        _core.Simulation(spec, particle_ids=[-1, 3])
    with pytest.raises(ValueError):
        _core.Simulation(spec, particle_ids=[])
    simulation = _core.Simulation(spec)
    with pytest.raises(ValueError):
        simulation.sample(1, stride=0)
    with pytest.raises(TypeError):
        simulation.advance(-1)
    spec.dt = -1
    with pytest.raises(ValueError):
        _core.Simulation(spec)


def test_metadata_records_physics_and_backend_limits():
    spec = _core.SimulationSpec(n_particles=2, dimension=1)
    metadata = _core.Simulation(spec).metadata
    assert metadata["schema"] == "bdmc.observation.v1"
    assert metadata["noise_policy"] == "legacy_uniform"
    assert metadata["initialization"] == "zero positions"
    assert metadata["time_unit"].startswith("tau_B")
    assert "not identical" in metadata["rng_compatibility"]
    assert "deferred" in metadata["scientific_validation"]
    assert spec.to_dict()["box_lengths"] == [10]
