"""Sparse spatial candidates remain separate from the simulation physics."""

import itertools

import numpy as np
import pytest

from bdmc import _core


def assert_backends_equal(positions, lengths, cutoff, ids=None):
    cell = _core.neighbor_pairs(positions, lengths, cutoff, particle_ids=ids)
    brute = _core.neighbor_pairs(positions, lengths, cutoff, "brute_force", particle_ids=ids)
    for name in ("pair_indices", "pair_ids", "distances", "vectors"):
        np.testing.assert_array_equal(cell[name], brute[name])
    n = len(positions)
    assert cell["stats"]["candidate_pairs"] <= n * (n - 1) // 2
    assert brute["stats"]["candidate_pairs"] == n * (n - 1) // 2
    assert cell["stats"]["occupied_cells"] <= n
    assert cell["stats"]["stored_particle_indices"] == n
    indices = cell["pair_indices"]
    assert np.all(indices[:, 0] < indices[:, 1])
    assert len(np.unique(indices, axis=0)) == len(indices)
    return cell


def test_analytic_geometry_ownership_and_readonly_outputs():
    positions = np.array([[0.0], [0.5], [1.0], [5.0], [9.0]])
    result = assert_backends_equal(positions, [10], 1.0, [99, 4, 72, 6, 1000])
    np.testing.assert_array_equal(result["pair_indices"], [[0, 1], [0, 2], [0, 4], [1, 2]])
    np.testing.assert_array_equal(result["pair_ids"], [[99, 4], [99, 72], [99, 1000], [4, 72]])
    np.testing.assert_array_equal(result["vectors"], [[0.5], [1.0], [-1.0], [0.5]])
    np.testing.assert_array_equal(result["distances"], [0.5, 1.0, 1.0, 0.5])
    for name in ("pair_indices", "pair_ids", "distances", "vectors"):
        assert result[name].flags.c_contiguous
        assert result[name].flags.owndata
        assert not result[name].flags.writeable
        with pytest.raises(ValueError):
            result[name].flat[0] = 123
    positions[:] = -777
    assert result["vectors"][0, 0] == 0.5


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_randomized_arbitrary_wrapping_matches_brute_and_numpy(dimension):
    generator = np.random.default_rng(8401 + dimension)
    lengths = np.arange(3, 3 + dimension, dtype=float)
    for cutoff in [0.03, 0.4, 1.5]:
        positions = generator.uniform(-3, 4, size=(90, dimension)) * lengths
        result = assert_backends_equal(positions, lengths, cutoff)
        wrapped = positions % lengths
        vectors = wrapped[None, :, :] - wrapped[:, None, :]
        vectors -= lengths * np.floor(vectors / lengths + 0.5)
        distances = np.linalg.norm(vectors, axis=-1)
        expected_i, expected_j = np.where(np.triu(distances <= cutoff, k=1))
        np.testing.assert_array_equal(result["pair_indices"], np.column_stack([expected_i, expected_j]))
        np.testing.assert_allclose(result["vectors"], vectors[expected_i, expected_j], atol=2e-14)


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_two_cell_axes_particles_on_boundaries_and_half_cutoff(dimension):
    coordinates = [-2.0, -np.finfo(float).eps, 0, np.nextafter(0.0, 1.0),
                   0.5, 1, np.nextafter(2.0, 0.0), 2, 4]
    positions = np.array(list(itertools.product(coordinates, repeat=dimension)))
    assert_backends_equal(positions, [2] * dimension, 1.0)
    assert_backends_equal(positions, [2] * dimension, np.nextafter(1.0, 0.0))


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_empty_single_and_coincident_particles(dimension):
    empty = assert_backends_equal(np.empty((0, dimension)), [1] * dimension, 0.5)
    assert empty["pair_indices"].shape == (0, 2)
    assert empty["vectors"].shape == (0, dimension)
    assert empty["distances"].shape == (0,)
    assert_backends_equal(np.zeros((1, dimension)), [1] * dimension, 0.5)
    overlap = assert_backends_equal(np.zeros((4, dimension)), [1] * dimension, np.nextafter(0., 1.))
    assert len(overlap["distances"]) == 6
    assert np.all(overlap["distances"] == 0)


def test_sparse_huge_grid_and_large_finite_geometry():
    positions = np.array([
        [0, 0, 0], [0.004, 0.002, 0], [1e12 - 0.003, 0, 0],
        [5e11, 1e12, 1.5e12], [-1e12, 4e12, -3e12], [0, 0, 0],
    ])
    result = assert_backends_equal(positions, [1e12, 2e12, 3e12], 0.01)
    assert result["stats"]["dense_grid_cells"] > np.iinfo(np.uint64).max
    assert result["stats"]["occupied_cells"] <= len(positions)
    assert result["stats"]["stored_particle_indices"] == len(positions)
    large = assert_backends_equal(
        np.array([[0, 0, 0], [0.5e307, 0.5e307, 0.5e307]]), [1e308] * 3, 1e307
    )
    assert len(large["distances"]) == 1
    assert np.isfinite(large["distances"][0])


def test_stable_ids_and_reordered_storage():
    positions = np.array([[0.1], [0.2], [1.0], [4.9], [9.9]])
    ids = np.array([91, 4, 72, 6, np.iinfo(np.uint64).max], dtype=np.uint64)
    original = assert_backends_equal(positions, [10], 1, ids)
    permutation = [3, 0, 4, 1, 2]
    changed = assert_backends_equal(positions[permutation], [10], 1, ids[permutation])

    def keyed(result):
        return {tuple(sorted(pair)): distance
                for pair, distance in zip(result["pair_ids"], result["distances"])}

    assert keyed(original) == keyed(changed)


@pytest.mark.parametrize("positions,lengths,cutoff,backend,ids", [
    (np.zeros((2, 1)), [1], 0, "cell_list", None),
    (np.zeros((2, 1)), [1], -1, "cell_list", None),
    (np.zeros((2, 1)), [1], 0.500001, "cell_list", None),
    (np.zeros((2, 1)), [1], np.inf, "cell_list", None),
    (np.zeros((2, 1)), [1], np.nan, "cell_list", None),
    (np.zeros((2, 1)), [1], 0.1, "unknown", None),
    (np.zeros((2, 1)), [0], 0.1, "cell_list", None),
    (np.zeros((2, 1)), [np.inf], 0.1, "cell_list", None),
    (np.zeros((2, 4)), [1] * 4, 0.1, "cell_list", None),
    (np.zeros((2, 2)), [1], 0.1, "cell_list", None),
    (np.array([[np.nan]]), [1], 0.1, "cell_list", None),
    (np.zeros((2, 1)), [1], 0.1, "cell_list", [4, 4]),
    (np.zeros((2, 1)), [1], 0.1, "cell_list", []),
    (np.zeros((2, 1)), [1], 0.1, "cell_list", [-1, 3]),
])
def test_invalid_search_inputs_rejected(positions, lengths, cutoff, backend, ids):
    with pytest.raises(ValueError):
        _core.neighbor_pairs(positions, lengths, cutoff, backend, ids)


def test_search_is_independent_of_fixed_bonds_and_rng():
    spec = _core.SimulationSpec(model="dumbbell", n_particles=2, dimension=1, seed=93)
    initial = np.array([[1.0], [4.0]])
    first, reference = (_core.Simulation(spec, initial) for _ in range(2))
    assert np.linalg.norm(first.observe()["forces"]) > 0
    result = _core.neighbor_pairs(first.observe()["positions"], [10], 0.1)
    assert result["pair_indices"].shape == (0, 2)
    np.testing.assert_array_equal(first.fixed_bonds, [[0, 1]])
    first.advance(15)
    reference.advance(15)
    np.testing.assert_array_equal(first.observe()["unwrapped_positions"], reference.observe()["unwrapped_positions"])
