"""Independent geometric oracles for radial/3D-directional/sheared observers."""

import itertools

import numpy as np
import pytest

from bdmc import _analysis


def radial_oracle(p, lengths, edges):
    p = np.asarray(p, dtype=float)
    if p.ndim == 2:
        p = p[None]
    d = p.shape[-1]
    counts = np.zeros(len(edges) - 1, dtype=np.uint64)
    for frame in p:
        wrapped = frame % lengths
        for i, j in itertools.combinations(range(len(frame)), 2):
            dr = wrapped[j] - wrapped[i]
            dr -= np.asarray(lengths) * np.floor(dr / lengths + .5)
            counts += (2 * np.histogram([np.linalg.norm(dr)], bins=edges)[0]).astype(np.uint64)
    unit_ball = (2., np.pi, 4 * np.pi / 3)[d - 1]
    volumes = unit_ball * np.diff(np.asarray(edges) ** d)
    expected = len(p) * p.shape[1] * (p.shape[1] - 1) * volumes / np.prod(lengths)
    return counts, volumes, counts / expected


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_radial_matches_independent_pair_loop_and_exact_shells(dimension):
    rng = np.random.default_rng(20261111 + dimension)
    lengths = np.arange(8, 8 + dimension, dtype=float)
    p = rng.uniform(-3, 4, (3, 31, dimension)) * lengths
    edges = [0, .4, 1., 2., 4.]
    counts, volumes, g = radial_oracle(p, lengths, edges)
    result = _analysis.pair_correlation(p, lengths, edges)
    brute = _analysis.pair_correlation(p, lengths, edges, backend="brute_force")
    np.testing.assert_array_equal(result["counts"], counts)
    np.testing.assert_array_equal(brute["counts"], counts)
    np.testing.assert_allclose(result["g"], g)
    np.testing.assert_allclose(result["bin_volumes"], volumes)
    np.testing.assert_allclose(result["cumulative_coordination"], np.cumsum(counts) / 3 / 31)
    for key in ["counts", "g", "bin_volumes", "expected_counts", "radial_edges", "box_lengths"]:
        assert result[key].flags.owndata
        assert result[key].flags.c_contiguous
        assert not result[key].flags.writeable
    before = result["g"].copy()
    p[:] = 0
    np.testing.assert_array_equal(result["g"], before)
    assert brute["candidate_pairs"] == 3 * 31 * 30 // 2


def test_finite_population_normalization_and_batch_counting():
    p = [[0, 0, 0], [1, 0, 0]]
    single = _analysis.pair_correlation(p, [10] * 3, [0, 1.5, 2])
    bulk = _analysis.pair_correlation(p, [10] * 3, [0, 1.5, 2], normalization="bulk")
    twice = _analysis.pair_correlation([p, p], [10] * 3, [0, 1.5, 2])
    np.testing.assert_array_equal(single["counts"], [2, 0])
    np.testing.assert_array_equal(twice["counts"], 2 * single["counts"])
    np.testing.assert_array_equal(twice["g"], single["g"])
    np.testing.assert_array_equal(bulk["g"], .5 * single["g"])
    np.testing.assert_array_equal(bulk["cumulative_coordination"], [1, 1])


def test_directional_sectors_azimuth_and_exact_angular_integration():
    p = [[0, 0, 0], [1, 1, 1]]
    r = _analysis.directional_pair_correlation(p, [10] * 3, [0, 2], [-1, 0, 1], [0, np.pi, 2 * np.pi])
    np.testing.assert_array_equal(r["counts"], [[[0, 1], [1, 0]]])
    np.testing.assert_allclose(r["bin_volumes"], 8 * np.pi / 3)
    radial = _analysis.pair_correlation(p, [10] * 3, [0, 2])
    np.testing.assert_array_equal(r["counts"].sum(axis=(1, 2)), radial["counts"])
    np.testing.assert_allclose((r["g"] * r["bin_volumes"]).sum(axis=(1, 2)), radial["g"] * radial["bin_volumes"])
    changed = _analysis.directional_pair_correlation(
        np.asarray(p)[:, [2, 0, 1]], [10] * 3, [0, 2], [-1, 0, 1], [0, np.pi, 2 * np.pi],
        polar_axis=[1, 0, 0], azimuth_reference=[0, 1, 0],
    )
    np.testing.assert_array_equal(changed["counts"], r["counts"])
    overlap = _analysis.directional_pair_correlation(np.zeros((2, 3)), [10] * 3, [0, 2])
    assert overlap["zero_distance_pairs"] == 2
    assert overlap["counts"].sum() == 0


def test_full_3d_planar_projection_and_finite_slab_normalization():
    p = [[0, 0, 0], [1, .5, .2]]
    args = (p, [10] * 3, [-2, 0, 2], [-1, 0, 1])
    full = _analysis.planar_pair_correlation(*args)
    np.testing.assert_array_equal(full["counts"], [[1, 0], [0, 1]])
    np.testing.assert_array_equal(full["bin_volumes"], np.full((2, 2), 20))
    slab = _analysis.planar_pair_correlation(*args, slab_half_width=.2)
    np.testing.assert_array_equal(slab["counts"], full["counts"])
    np.testing.assert_allclose(slab["g"], 25 * full["g"])
    omitted = _analysis.planar_pair_correlation(*args, slab_half_width=.1)
    assert not omitted["counts"].any()


def test_sheared_elongated_lattice_multiple_crossings_against_image_enumeration():
    lengths, shift = np.array([40., 8., 6.]), 127.
    p = np.array([[1., 0, 1], [8., 7.8, 1]])
    moved = p.copy()
    moved[1] += 5 * np.array([shift, lengths[1], 0]) + [-3 * lengths[0], 0, -4 * lengths[2]]
    for frame in [p, moved]:
        images = np.array([
            [mx * lengths[0] + my * shift, my * lengths[1], mz * lengths[2]]
            for mx in range(-64, 65) for my in range(-12, 13) for mz in range(-8, 9)
        ])
        nearest = np.linalg.norm(frame[1] - frame[0] - images, axis=1).min()
        assert nearest == pytest.approx(.2)
        r = _analysis.pair_correlation(frame, lengths, [0, .3, 1, 3], shear_displacement=shift)
        assert r["backend"] == "brute_force"
        np.testing.assert_array_equal(r["counts"], [2, 0, 0])
    args = (lengths, [0, .3, 1, 3], [-1, 1], [0, np.pi, 2 * np.pi])
    a = _analysis.directional_pair_correlation(p, *args, shear_displacement=shift)
    b = _analysis.directional_pair_correlation(moved, *args, shear_displacement=shift)
    np.testing.assert_array_equal(a["counts"], b["counts"])
    x, y = [-1, .1, 1], [-1, 0, 1]
    a = _analysis.planar_pair_correlation(p, lengths, x, y, shear_displacement=shift)
    b = _analysis.planar_pair_correlation(moved, lengths, x, y, shear_displacement=shift)
    np.testing.assert_array_equal(a["counts"], b["counts"])


def test_fixed_seed_ideal_gas_radial_and_solid_angle_normalization():
    p = np.random.default_rng(2026111101).uniform(0, 10, (40, 64, 3))
    r = _analysis.directional_pair_correlation(p, [10] * 3, [0, 1.5, 2.5, 4], [-1, 0, 1], [0, np.pi, 2 * np.pi])
    np.testing.assert_allclose(r["g"], 1., atol=.12, rtol=0)
    projected = _analysis.planar_pair_correlation(p, [10] * 3, [-5, 0, 5], [-5, 0, 5])
    np.testing.assert_allclose(projected["g"], 1., atol=.04, rtol=0)


@pytest.mark.parametrize("kwargs", [
    {"radial_edges": [0, 0]}, {"radial_edges": [-1, 1]}, {"radial_edges": [0, 6]},
    {"radial_edges": [0, np.nan]}, {"normalization": "unknown"}, {"backend": "unknown"},
    {"shear_displacement": 1, "backend": "cell_list"}, {"shear_displacement": np.inf},
])
def test_invalid_radial_contracts(kwargs):
    args = {"positions": [[0, 0, 0], [1, 0, 0]], "box_lengths": [10] * 3, "radial_edges": [0, 1]}
    args.update(kwargs)
    with pytest.raises(ValueError):
        _analysis.pair_correlation(**args)


def test_invalid_positions_angles_and_projection():
    for p in [np.zeros((1, 3)), np.zeros((0, 3)), np.zeros((0, 2, 3)), np.zeros(6), [[np.nan, 0, 0], [0, 0, 0]]]:
        with pytest.raises(ValueError):
            _analysis.pair_correlation(p, [10] * 3, [0, 1])
    p = [[0, 0, 0], [1, 1, 1]]
    for kwargs in [{"mu_edges": [-2, 1]}, {"phi_edges": [0, 7]}, {"polar_axis": [0, 0, 0]},
                   {"azimuth_reference": [0, 0, 1]}]:
        with pytest.raises(ValueError):
            _analysis.directional_pair_correlation(p, [10] * 3, [0, 2], **kwargs)
    for kwargs in [{"slab_half_width": 0}, {"slab_half_width": 6}, {"flow_axis": 1, "gradient_axis": 1}]:
        with pytest.raises(ValueError):
            _analysis.planar_pair_correlation(p, [10] * 3, [-1, 1], [-1, 1], **kwargs)


def test_native_autocorrelation_matches_independent_product_sums():
    values = np.array([[1., 5], [2, 4], [3, 1], [4, 0], [7, -1]])
    for center in [False, True]:
        data = values - values.mean(axis=0) if center else values
        expected = np.array([(data[:len(data)-lag] * data[lag:]).mean(axis=0) for lag in range(5)])
        result = _analysis.autocorrelation(values, 4, center=center)
        np.testing.assert_allclose(result["values"], expected)
        np.testing.assert_array_equal(result["origin_counts"], [5, 4, 3, 2, 1])
    old = _analysis.autocorrelation(values[:, 0], 1, window="legacy_fixed_window")
    np.testing.assert_array_equal(old["values"], [1, 2])
    for kwargs in [{"max_lag": -1}, {"max_lag": 5}, {"window": "unknown"}]:
        with pytest.raises(ValueError):
            _analysis.autocorrelation(values, **kwargs)
