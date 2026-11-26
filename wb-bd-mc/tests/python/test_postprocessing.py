"""Independent physical identities for corrected historical processing purposes."""
from itertools import product
from pathlib import Path
import numpy as np
import pytest

from bdmc import postprocessing as pp
from bdmc.analysis import pair_correlation, coordination_number, autocorrelation


def adjacency(n, edges=(), loops=()):
    rows = [{i: 0} for i in range(n)]
    for i, count in loops:
        rows[i][i] = 2 * count
    for i, j, weight in edges:
        rows[i][j] = weight; rows[j][i] = weight
    return tuple((np.array(list(row), dtype=int), np.array(list(row.values()), dtype=int)) for row in rows)


def test_topology_chains_distinct_edges_loops_and_components():
    a = adjacency(4, [(0, 3, 5)], [(0, 2), (1, 1)])
    stats = pp.topology_statistics(a)
    assert stats['bridge_chains'] == 5 and stats['loop_chains'] == 3
    assert stats['total_chains'] == 8 and stats['distinct_edges'] == 1
    assert stats['bridge_fraction'] == 5/8
    np.testing.assert_array_equal(stats['attached_ends'], [9, 2, 0, 5])
    np.testing.assert_array_equal(stats['cluster_sizes'], [2, 1, 1])
    pos = [[0., 0, 0], [3., 0, 0], [4., 0, 0], [1., 0, 0]]
    chain = pp.bond_observables(pos, a, [10]*3)
    distinct = pp.bond_observables(pos, a, [10]*3, measure='distinct_edges')
    assert chain['RR_sum'][0, 0] == 5 and distinct['RR_sum'][0, 0] == 1
    assert chain['RR_per_particle'][0, 0] == 1.25
    assert chain['RR_per_bridge'][0, 0] == 1
    np.testing.assert_allclose(chain['RR_sum'], np.diag([5, 0, 0]))


@pytest.mark.parametrize('shift', [-13., -3., 0., 4., 14.])
def test_exact_shear_images_against_independent_lattice_enumeration(shift):
    box = np.array([12., 5., 7.]); rng = np.random.default_rng(517)
    vectors = rng.uniform(-25, 25, (8, 3))
    actual = pp.minimum_image(vectors, box, shear_displacement=shift)
    lattice = np.array([[i*box[0]+j*shift, j*box[1], k*box[2]]
                        for i, j, k in product(range(-12, 13), repeat=3)])
    for v, a in zip(vectors, actual):
        expected = np.min(np.linalg.norm(v - lattice, axis=1))
        assert np.linalg.norm(a) == pytest.approx(expected, abs=1e-12)
    moved = vectors + np.array([2*box[0] + 3*shift, 3*box[1], -box[2]])
    np.testing.assert_allclose(pp.minimum_image(moved, box, shear_displacement=shift), actual, atol=1e-12)


def test_shear_short_radial_counts_and_angular_marginal():
    p = np.array([[1., .1, .2], [3.9, 9.9, .2], [6., 5., 5.]])
    g = pair_correlation(p, [10]*3, [0., .5, 1., 2.], shear_displacement=3.)
    d = pp.minimum_image(p[1]-p[0], [10]*3, shear_displacement=3.)
    assert np.linalg.norm(d) == pytest.approx(np.sqrt(.05))
    assert g['counts'][0] == 2
    assert coordination_number(g, .5) == pytest.approx(2/3)


def test_weighted_pdf_and_spherical_exact_sector_integral():
    radial = pp.histogram_density([0., .2, 1., 3.], [0., .5, 2.], weights=[2, 3, 5, 7])
    np.testing.assert_array_equal(radial['counts'], [5, 5])
    assert radial['outside_weight'] == 7
    assert np.sum(radial['pdf'] * radial['bin_measure']) == pytest.approx(10/17)
    q = np.array([[1, 0, 0], [-1, 0, 0], [0, 0, 1], [0, 0, 0.]])
    s = pp.spherical_bond_density(q, [0, .5, 2], [-1, 0, 1], [0, np.pi, 2*np.pi], weights=[1, 2, 3, 4])
    assert s['zero_vector_weight'] == 4
    assert s['bin_measure'].sum() == pytest.approx(4*np.pi*2**3/3)
    assert np.sum(s['pdf'] * s['bin_measure']) == pytest.approx(.6)
    assert s['counts'].sum() == 6
    empty = pp.bond_observables([[0, 0, 0]], adjacency(1, loops=[(0, 1)]), [4]*3)
    assert empty['bridge_weight'] == 0 and empty['RR_per_bridge'] is None


@pytest.mark.parametrize('n', [5, 8, 9])
@pytest.mark.parametrize('window', ['all_origins', 'legacy_fixed_window'])
@pytest.mark.parametrize('center', [False, True])
def test_linear_fft_direct_explicit_sums_and_current_native(n, window, center):
    rng = np.random.default_rng(n); x = rng.normal(size=(n, 2)); y = rng.normal(size=(n, 2))
    lag = min(2, n//2-1)
    a = pp.correlation(x, y, times=10 + np.arange(n)*.25, max_lag=lag, window=window, center=center)
    xc = x - x.mean(axis=0) if center else x
    yc = y - y.mean(axis=0) if center else y
    expected = []
    for k in range(lag + 1):
        m = n-k if window == 'all_origins' else n//2-1
        expected.append(sum(xc[j]*yc[j+k] for j in range(m)) / m)
    np.testing.assert_allclose(a['values'], expected, atol=2e-15)
    np.testing.assert_allclose(a['lag_times'], np.arange(lag+1)*.25)
    b = pp.correlation(x, y, max_lag=lag, window=window, center=center, method='direct')
    np.testing.assert_allclose(a['values'], b['values'], atol=2e-15)
    c = pp.correlation(x, max_lag=lag, window=window, center=center)
    native = autocorrelation(x, max_lag=lag, window=window, center=center)
    np.testing.assert_allclose(c['values'], native['values'], atol=2e-15)
    np.testing.assert_array_equal(c['origin_counts'], native['origin_counts'])


def test_isotropic_stress_projection_rotations_cross_terms_and_pressure_offset():
    rng = np.random.default_rng(515); a = rng.normal(size=(10, 6)); b = rng.normal(size=(10, 6))
    def matrices(s):
        return np.array([[[x[0], x[3], x[4]], [x[3], x[1], x[5]], [x[4], x[5], x[2]]] for x in s])
    def columns(m):
        return m[:, [0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]
    am = matrices(a); bm = matrices(b)
    am -= np.trace(am, axis1=1, axis2=2)[:, None, None]*np.eye(3)/3
    bm -= np.trace(bm, axis1=1, axis2=2)[:, None, None]*np.eye(3)/3
    expected = [np.sum(am[:10-k]*bm[k:]) / (10*(10-k)) for k in range(4)]
    actual = pp.stress_correlation(a, b, center=False, max_lag=3)
    np.testing.assert_allclose(actual['values'], expected, atol=1e-15)
    angle = .87; rot = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1.]])
    ar = columns(rot @ matrices(a) @ rot.T); br = columns(rot @ matrices(b) @ rot.T)
    np.testing.assert_allclose(pp.stress_correlation(ar, br, center=False, max_lag=3)['values'], expected, atol=1e-15)
    shifted = a.copy(); shifted[:, :3] += np.arange(10)[:, None]*100
    np.testing.assert_allclose(pp.stress_correlation(shifted, b, center=False, max_lag=3)['values'], expected, atol=1e-13)
    total = pp.stress_correlation(a+b, max_lag=3)['values']
    parts = sum(pp.stress_correlation(x, y, max_lag=3)['values'] for x, y in [(a,a), (b,b), (a,b), (b,a)])
    np.testing.assert_allclose(total, parts, atol=1e-15)


def test_green_kubo_units_and_finite_running_integral():
    corr = {'values': np.array([2., 1., 0.]), 'lag_times': np.array([0., .5, 1.]), 'centered': True}
    g = pp.green_kubo(corr, volume=4, thermal_energy=2, equilibrium=True)
    np.testing.assert_allclose(g['modulus_contribution'], [4, 2, 0])
    np.testing.assert_allclose(g['running_viscosity_contribution'], [0, 1.5, 2])
    with pytest.raises(ValueError, match='equilibrium'):
        pp.green_kubo(corr, volume=4, thermal_energy=2, equilibrium=False)


def test_unwrapping_msd_origins_image_counters_and_affine_subtraction():
    p = np.array([[[9., 1]], [[.5, 2]], [[2., 3]], [[3.5, 4]]])
    with pytest.raises(ValueError, match='image counters'):
        pp.unwrap_positions(p, [10, 10])
    u = pp.unwrap_positions(p, [10, 10], assume_small_steps=True)
    np.testing.assert_allclose(u[:, 0, 0], [9, 10.5, 12, 13.5])
    for origin in ('all', 'first'):
        result = pp.mean_square_displacement(u, times=np.arange(4)*.5 + 10, origin=origin)
        np.testing.assert_allclose(result['msd'], np.arange(4)**2*3.25)
        np.testing.assert_allclose(result['lag_times'], np.arange(4)*.5)
    flow = np.zeros_like(u); flow[:, :, 0] = np.arange(4)[:, None]*1.5
    np.testing.assert_allclose(pp.mean_square_displacement(u, affine_displacement=flow)['msd'], np.arange(4)**2)
    im = np.zeros_like(p, dtype=int); im[1:, :, 0] = 1
    np.testing.assert_array_equal(pp.unwrap_positions(p, [10, 10], images=im), u)
    with pytest.raises(ValueError, match='sheared'):
        pp.unwrap_positions(p, [10, 10], assume_small_steps=True, shear_displacements=[0]*4)


def test_structure_factor_self_forward_peak_translation_and_2d_transform():
    p = np.array([[0., 0], [1., 0]])
    q = np.array([[0., 0], [np.pi, 0], [2*np.pi, 0]])
    np.testing.assert_allclose(pp.structure_factor(p, q), [2, 0, 2], atol=1e-30)
    np.testing.assert_allclose(pp.structure_factor(p+7.3, q), [2, 0, 2], atol=1e-28)
    r = np.linspace(0, 1, 1001)
    assert pp.rdf_to_structure_factor(r, np.ones_like(r), [0, 1], number_density=.5)['values'].tolist() == [1., 1.]
    threed = pp.rdf_to_structure_factor(r, np.ones_like(r)*2, [0], number_density=.5)
    assert threed['values'][0] == pytest.approx(1+2*np.pi/3, rel=1e-6)
    pytest.importorskip('scipy.special')
    twod = pp.rdf_to_structure_factor(r, np.ones_like(r)*2, [0], number_density=.5, dimension=2)
    assert twod['values'][0] == pytest.approx(1+np.pi/2, abs=1e-12)


@pytest.mark.parametrize('law', ['gaussian', 'modified_gaussian', 'fene'])
def test_force_energy_derivative_virial_sign_and_all_unlinked_repulsive_pairs(law):
    # Three pairs, only 0--1 is bonded; particle 2 still repels both particles.
    pos = np.array([[0., 0, 0], [.5, 0, 0], [0, .5, 0]])
    adj = adjacency(3, [(0, 1, 2)], [(2, 1)])
    kwargs = dict(connector=law, alpha=1.5, repulsion_coefficient=25, fene_maximum=2.)
    result = pp.conservative_virial(pos, adj, [10]*3, **kwargs)
    # Isotropic scaling derivative: trace(RF/V)*V = dU/d(log length).
    h = 1e-6
    plus = pp.conservative_virial(pos*np.exp(h), adj, [10]*3, **kwargs)
    minus = pp.conservative_virial(pos*np.exp(-h), adj, [10]*3, **kwargs)
    for part in ('repulsive', 'connector'):
        derivative = (plus[part+'_energy'] - minus[part+'_energy'])/(2*h)
        assert 1000*np.trace(result[part]) == pytest.approx(derivative, rel=1e-9)
    assert np.trace(result['repulsive']) < 0 and np.trace(result['connector']) > 0
    assert result['repulsive_energy'] > (2-3*.5+.5**3)/3


def write_frames(prefix, frames):
    for extension in ('.traj', '.hash', '.weight', '.chain'):
        Path(str(prefix)+extension).write_text('')
    for t, positions, adj, ends in frames:
        aux = np.zeros_like(positions)
        row = np.concatenate([[t], np.concatenate([positions, aux], axis=1).ravel()])
        with Path(str(prefix)+'.traj').open('a') as f:
            f.write(' '.join(map(str, row))+'\n')
        for suffix, column in [('.hash', 0), ('.weight', 1)]:
            with Path(str(prefix)+suffix).open('a') as f:
                for pair in adj:
                    f.write(' '.join(map(str, pair[column]))+'\n')
        with Path(str(prefix)+'.chain').open('a') as f:
            f.write(' '.join(map(str, np.asarray(ends).ravel()))+'\n')


def test_synchronized_ragged_io_stride_labeled_chain_midpoint_and_rejection(tmp_path):
    prefix = tmp_path/'frame'; a = adjacency(2, [(0, 1, 1)])
    frames = [(float(t), np.array([[9., 0, 0], [1., 0, 0]]), a, [[0], [1]]) for t in range(5)]
    write_frames(prefix, frames)
    read = list(pp.iter_legacy_frames(prefix, 2, start=1, stride=2, with_chains=True))
    assert [r.index for r in read] == [1, 3]
    np.testing.assert_allclose(pp.chain_midpoints(read[0], [10]*3)['midpoints'], [[10., 0, 0]])
    Path(str(prefix)+'.weight').write_text('0 1\n0 2\n')
    with pytest.raises(ValueError, match='reciprocal'):
        list(pp.iter_legacy_frames(prefix, 2))
    write_frames(prefix, frames)
    with Path(str(prefix)+'.hash').open('a') as f:
        f.write('0\n')
    with pytest.raises(ValueError, match='extra frames'):
        list(pp.iter_legacy_frames(prefix, 2))
    write_frames(prefix, [frames[0], frames[1], frames[0]])
    with pytest.raises(ValueError, match='clock'):
        list(pp.iter_legacy_frames(prefix, 2))
    split = list(pp.iter_legacy_frames(prefix, 2, allow_time_resets=True))
    assert [f.segment for f in split] == [0, 0, 1]


def test_explicit_energy_schemas_and_distance_stream(tmp_path):
    tags = np.tile(np.arange(31), (4, 1)).astype(float)
    r = pp.stress_components(tags, schema='heur31')
    np.testing.assert_array_equal(r['configurational'][0], np.arange(18, 24)+np.arange(24, 30))
    with pytest.raises(ValueError, match='schema'):
        pp.stress_components(tags, schema='dumbbell24')
    path = tmp_path/'dist'; path.write_text('0 1\n2\n')
    assert list(pp.distance_stream(path)) == [0, 1, 2]
    path.write_text('7.829628e-018.581755e-01')
    np.testing.assert_array_equal(list(pp.distance_stream(path, legacy_scientific_printer=True)), [.7829628, .8581755])
    with pytest.raises(ValueError):
        list(pp.distance_stream(path))


def test_labeled_connectors_empirical_cdf_bell_factor_and_einstein_fit():
    q = np.arange(30.).reshape(5, 2, 3)/10
    c = pp.connector_correlation(q, max_lag=4)
    expected = np.array([np.mean(np.sum(q[:5-k]*q[k:], axis=2), axis=0) for k in range(5)])
    np.testing.assert_allclose(c['per_chain'], expected, atol=1e-15)
    np.testing.assert_allclose(c['values'], expected.mean(axis=1))
    e = pp.empirical_cdf([2, 1, 2, 0], weights=[1, 2, 3, 4])
    np.testing.assert_allclose(e['x'], [0, 1, 2])
    np.testing.assert_allclose(e['cdf'], [.4, .6, 1])
    assert pp.bell_factor([0, 1], capture_length=0) == 1
    assert pp.bell_factor([0, 1], alpha=1, capture_length=.1, weights=[1, 3]) == pytest.approx((1+3*np.exp(.3))/4)
    t = np.arange(10.)
    fit = pp.diffusion_fit(6*.2*t + .1, t, dimension=3, fit_range=[2, 8])
    assert fit['diffusion'] == pytest.approx(.2)
    assert fit['intercept'] == pytest.approx(.1)


def test_image_counters_for_both_shear_signs_and_geometric_repulsion():
    pos = np.zeros((2, 1, 3)); images = np.zeros_like(pos, dtype=int)
    images[:, :, 1] = 2
    unwrapped = pp.unwrap_positions(pos, [12, 10, 10], images=images, shear_displacements=[-3., 4.])
    np.testing.assert_allclose(unwrapped[:, 0, 0], [-6, 8])
    np.testing.assert_allclose(unwrapped[:, 0, 1], [20, 20])
    pos = [[0., 0, 0], [.5, 0, 0]]; a = adjacency(2, loops=[(0, 2), (1, 3)])
    geometric = pp.conservative_virial(pos, a, [4]*3, repulsion_coupling='geometric_mean',
                                       repulsion_base=.25, repulsion_coefficient=2.)
    average = pp.conservative_virial(pos, a, [4]*3, repulsion_coefficient=2.)
    np.testing.assert_allclose(geometric['repulsive'], 3*average['repulsive'])


def test_legacy_geometry_actual_clock_rate_and_rectangular_lengths():
    p = dict(N_dimension='3', box_dimension='10', LBx='12', LBy='8', LBz='6',
             repulsion_coefficient='25', MECHANICAL_PERTURBATION='FALSE')
    base = pp.legacy_geometry(p, 2.)
    np.testing.assert_array_equal(base['box_lengths'], [12, 8, 6])
    assert base['shear_displacement'] == 0
    for rate in [-.2, .2]:
        p.update(MECHANICAL_PERTURBATION='SIMPLE_SHEAR', Wi_tau_C=str(rate))
        assert pp.legacy_geometry(p, 2.)['shear_displacement'] == pytest.approx(rate*8*2)
    for strain in [-.3, .3]:
        p.update(MECHANICAL_PERTURBATION='STEP_SHEAR', gamma_0=str(strain))
        assert pp.legacy_geometry(p, 99.)['shear_displacement'] == pytest.approx(strain*8)
    p['shear_axis'] = '2'
    with pytest.raises(ValueError, match='xy only'):
        pp.legacy_geometry(p, 0.)


def test_absent_chain_fraction_and_negative_diffusion_are_undefined():
    assert pp.topology_statistics(adjacency(1))['bridge_fraction'] is None
    with pytest.raises(ValueError, match='negative fitted slope'):
        pp.diffusion_fit([3., 2., 1.], [0., 1., 2.], dimension=3, fit_range=[0, 2])


@pytest.mark.parametrize('call', [
    lambda: pp.correlation([1,2,3], times=[0,0,1]),
    lambda: pp.correlation([1,2,3], max_lag=True),
    lambda: pp.minimum_image([[1,2]], [1,-1]),
    lambda: pp.histogram_density([1], [0,0,2]),
    lambda: pp.validate_adjacency(((np.array([0]), np.array([1])),)),
    lambda: pp.conservative_virial([[0.,0,0]], adjacency(1), [2]*3, effective_distance=2),
    lambda: pp.bell_factor([1e6], alpha=1, capture_length=1),
])
def test_invalid_physical_domains_are_explicit(call):
    with pytest.raises(ValueError):
        call()
