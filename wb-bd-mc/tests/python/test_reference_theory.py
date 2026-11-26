"""Independent identities/domain controls for the reconstructed theory oracles.

These tests do not use a native force, integrator or native theory function as
expected data. Native-to-theory comparisons belong to the integration tests.
"""

import math

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from bdmc import reference_theory as theory


HOOKEAN = dict(diffusion=.75, mobility=.5, spring_constant=2.)
FENE = dict(spring_constant=2., q0=3., thermal_energy=.8)
PAIR = dict(epsilon=4., sigma=1.3, thermal_energy=.7, r_max=2.)


def test_single_bead_and_relative_brownian_noise_are_distinct():
    time = np.array([0., .5, 3.])
    np.testing.assert_allclose(theory.brownian_msd(time, diffusion=2), 12 * time)
    np.testing.assert_allclose(theory.brownian_msd(time, diffusion=2, relative=True),
                               24 * time)
    assert theory.thermal_energy(diffusion=.75, mobility=.5) == 1.5


def test_hookean_equilibrium_units_and_moment_relaxation_factor():
    m = theory.hookean_equilibrium_moment(**HOOKEAN)
    np.testing.assert_array_equal(m, .75 * np.eye(3))
    assert theory.hookean_relaxation_time(mobility=.5, spring_constant=2) == .25
    np.testing.assert_array_equal(theory.hookean_moment_rhs(m, **HOOKEAN), np.zeros((3, 3)))
    np.testing.assert_array_equal(theory.hookean_polymer_stress(m, **HOOKEAN), np.zeros((3, 3)))


def test_hookean_rhs_and_euler_step_against_hand_calculation():
    m = np.array([[2., .5], [.5, 3.]])
    k = np.array([[0., 2.], [0., 0.]])
    np.testing.assert_allclose(theory.hookean_moment_rhs(m, velocity_gradient=k, **HOOKEAN),
                               [[-3., 4.], [4., -9.]])
    actual = theory.hookean_em_moment_step(m, .1, velocity_gradient=k, **HOOKEAN)
    np.testing.assert_allclose(actual, [[1.86, .8], [.8, 2.22]])
    # Negative controls: moment-forward-Euler omits dt^2 B M B.T, while using
    # one-bead noise instead of relative noise misses .15 on both diagonals.
    assert not np.allclose(actual, m + .1 * np.array([[-3., 4.], [4., -9.]]))
    assert not np.allclose(actual, np.array([[1.86, .8], [.8, 2.22]]) - .15 * np.eye(2))


@pytest.mark.parametrize("steps", [0, 1, 7, 500, 10**15])
def test_euler_lifted_power_against_scalar_geometric_sum(steps):
    # q'=a*q+eta, a=.8 and E[eta^2]=.3: stationary moment=.3/(1-.64).
    a_squared = .64
    power = a_squared**steps
    expected = 3 * power + .3 / (1 - a_squared) * (1 - power)
    actual = theory.hookean_em_moment([[3.]], steps, .1, **HOOKEAN)
    np.testing.assert_allclose(actual, [[expected]], rtol=2e-13)


@pytest.mark.parametrize("dt", [.001, .1, .9])
def test_euler_equilibrium_has_matched_time_step_bias(dt):
    equilibrium = theory.hookean_em_equilibrium_moment(dt, **HOOKEAN, dimension=2)
    np.testing.assert_allclose(equilibrium, .75 / (1 - dt) * np.eye(2))
    np.testing.assert_allclose(theory.hookean_em_moment_step(equilibrium, dt, **HOOKEAN),
                               equilibrium)
    assert equilibrium[0, 0] > .75


@pytest.mark.parametrize("dt", [0., -1., 1., 1.1, float("inf")])
def test_euler_equilibrium_rejects_unstable_or_invalid_step(dt):
    with pytest.raises(ValueError):
        theory.hookean_em_equilibrium_moment(dt, **HOOKEAN)


def test_unstable_euler_endpoint_is_an_algebraic_oracle_not_stationarity():
    actual = theory.hookean_em_moment([[3.]], 2, 2., **HOOKEAN)
    # a=-3, noise variance=6, giving 3 -> 33 -> 303.
    np.testing.assert_allclose(actual, [[303.]])


def test_general_constant_flow_against_independent_ode_integration():
    m = np.array([[2., .3, -.1], [.3, 1., .2], [-.1, .2, .7]])
    k = np.array([[.1, .8, 0.], [-.2, -.1, .3], [0., -.2, .05]])
    # Spell out the scalar matrix components, rather than calling the oracle's
    # RHS or its matrix exponential to generate expected values.
    def rhs(_, values):
        matrix = values.reshape(3, 3)
        result = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                result[i, j] = sum(k[i, a] * matrix[a, j] + matrix[i, a] * k[j, a]
                                   for a in range(3)) - 4 * matrix[i, j] + 3 * (i == j)
        return result.ravel()

    time = np.array([[0., .07], [.3, 1.2]])
    solution = solve_ivp(rhs, (0, 1.2), m.ravel(), t_eval=time.ravel(),
                         rtol=2e-12, atol=2e-13, method="DOP853")
    assert solution.success
    expected = solution.y.T.reshape(time.shape + (3, 3))
    actual = theory.hookean_moment_evolution(m, time, velocity_gradient=k, **HOOKEAN)
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-12)


def test_continuum_singular_generator_and_zero_time():
    m = np.array([[2., .3], [.3, 1.]])
    # K=2*mu*H*I cancels connector drift, leaving free relative diffusion.
    actual = theory.hookean_moment_evolution(m, .4, velocity_gradient=2*np.eye(2), **HOOKEAN)
    np.testing.assert_allclose(actual, m + 1.2*np.eye(2))
    np.testing.assert_array_equal(theory.hookean_moment_evolution(m, 0., **HOOKEAN), m)


def test_continuum_noflow_long_time_and_closed_relaxation_agree():
    m = np.array([[2., .3], [.3, 1.]])
    times = [0., .3, 1000.]
    expected = .75*np.eye(2) + np.exp(-4*np.asarray(times))[:, None, None] * (m-.75*np.eye(2))
    np.testing.assert_allclose(theory.hookean_relaxation_moment(m, times, **HOOKEAN), expected)
    np.testing.assert_allclose(theory.hookean_moment_evolution(m, times, **HOOKEAN), expected,
                               atol=2e-13)


def test_finite_step_oracle_converges_to_continuum_at_first_order():
    m = np.array([[2., .3], [.3, 1.]])
    k = np.array([[0., 2.], [0., 0.]])
    exact = theory.hookean_moment_evolution(m, .4, velocity_gradient=k, **HOOKEAN)
    errors = [np.linalg.norm(theory.hookean_em_moment(m, n, .4/n,
               velocity_gradient=k, **HOOKEAN) - exact) for n in (20, 40, 80)]
    assert 1.9 < errors[0] / errors[1] < 2.1
    assert 1.9 < errors[1] / errors[2] < 2.1


@pytest.mark.parametrize("rate", [-8., 0., 8.])
@pytest.mark.parametrize("dimension", [2, 3])
def test_steady_shear_signed_stress_and_stationarity(rate, dimension):
    m = theory.hookean_steady_shear_moment(rate, **HOOKEAN, dimension=dimension)
    wi = rate / 4
    expected = .75*np.eye(dimension)
    expected[0, 0] = .75*(1+2*wi**2)
    expected[0, 1] = expected[1, 0] = .75*wi
    np.testing.assert_allclose(m, expected)
    k = np.zeros((dimension, dimension))
    k[0, 1] = rate
    np.testing.assert_allclose(theory.hookean_moment_rhs(m, velocity_gradient=k, **HOOKEAN),
                               np.zeros_like(m), atol=1e-13)
    stress = theory.hookean_polymer_stress(m, number_density=2, **HOOKEAN)
    assert stress[0, 1] == pytest.approx(3*wi)
    assert stress[0, 0]-stress[1, 1] == pytest.approx(6*wi**2)


def test_affine_step_shear_then_continuum_relaxation():
    f = np.eye(3)
    f[0, 1] = -.6
    initial = .75*np.eye(3)
    stepped = theory.hookean_step_deformation(initial, f)
    expected = .75*np.array([[1.36, -.6, 0.], [-.6, 1., 0.], [0., 0., 1.]])
    np.testing.assert_allclose(stepped, expected)
    relaxed = theory.hookean_relaxation_moment(stepped, .25, **HOOKEAN)
    stress = theory.hookean_polymer_stress(relaxed, number_density=2, **HOOKEAN)
    assert stress[0, 1] == pytest.approx(-1.8/math.e)
    assert stress[0, 0]-stress[1, 1] == pytest.approx(1.08/math.e)


def test_saos_moduli_against_relaxation_kernel_integrals():
    frequency = 1.7
    # G(t)=n*kT*exp(-t/lambda)=3*exp(-4t). Fourier constitutive
    # integrals yield G'=omega*int G(t)sin(omega*t)dt, and G'' with cos.
    storage = frequency * quad(lambda t: 3*math.exp(-4*t)*math.sin(frequency*t), 0, np.inf)[0]
    loss = frequency * quad(lambda t: 3*math.exp(-4*t)*math.cos(frequency*t), 0, np.inf)[0]
    actual = theory.hookean_saos_moduli(frequency, number_density=2, **HOOKEAN)
    np.testing.assert_allclose(actual, [storage, loss], rtol=2e-11)


def test_saos_zero_crossover_and_high_frequency_limits():
    storage, loss = theory.hookean_saos_moduli([0., 4., 4e12], number_density=2, **HOOKEAN)
    np.testing.assert_allclose(storage, [0., 1.5, 3.])
    np.testing.assert_allclose(loss, [0., 1.5, 3e-12], atol=1e-20)
    low_storage, low_loss = theory.hookean_saos_moduli(1e-7, number_density=2, **HOOKEAN)
    assert low_loss / 1e-7 == pytest.approx(.75)
    assert low_storage < 1e-13


@pytest.mark.parametrize("dimension", [1, 2, 3])
@pytest.mark.parametrize("b", [.5, 7., 100.])
def test_true_fene_radial_normalization_and_second_moment(dimension, b):
    # q0=2, kT=1 -> H=b/4. Integrals use adaptive quadrature on the
    # returned density; the expected second moment is the beta-law identity.
    params = dict(spring_constant=b/4, q0=2., thermal_energy=1., dimension=dimension)
    pdf = lambda q: theory.fene_radial_pdf(q, **params)
    normalization = quad(pdf, 0., 2., epsabs=2e-11, epsrel=2e-11)[0]
    second_moment = quad(lambda q: q*q*pdf(q), 0., 2., epsabs=2e-11, epsrel=2e-11)[0]
    assert normalization == pytest.approx(1., abs=3e-10)
    expected = dimension * 4 / (b + dimension + 2)
    assert second_moment == pytest.approx(expected, rel=3e-10)
    assert theory.fene_mean_square_extension(**params) == pytest.approx(expected)


def test_fene_pdf_includes_radial_jacobian_and_true_fene_denominator():
    q = np.array([.5, 1.])
    actual = theory.fene_radial_pdf(q, **FENE, dimension=3)
    b = 2*3**2/.8
    expected_ratio = (q[0]/q[1])**2 * ((1-q[0]**2/9)/(1-q[1]**2/9))**(b/2)
    assert actual[0]/actual[1] == pytest.approx(expected_ratio)
    true_value = theory.fene_mean_square_extension(**FENE)
    assert true_value == pytest.approx(27/(b+5))
    # The preaveraged FENE-P denominator b+d is a deliberate negative control.
    assert not math.isclose(true_value, 27/(b+3), rel_tol=.01)
    assert theory.fene_radial_pdf(0., **FENE, dimension=3) == 0.
    assert theory.fene_radial_pdf(0., **FENE, dimension=1) > 0.
    np.testing.assert_array_equal(theory.fene_radial_pdf([3., 4.], **FENE), [0., 0.])


def test_fene_force_derivative_and_equipartition():
    q = np.array([.1, .8, 2.5])
    delta = 1e-6
    force = theory.fene_force_magnitude(q, spring_constant=2., q0=3.)
    derivative = (theory.fene_potential(q+delta, spring_constant=2., q0=3.)
                  - theory.fene_potential(q-delta, spring_constant=2., q0=3.)) / (2*delta)
    np.testing.assert_allclose(force, derivative, rtol=2e-9)
    weighted = lambda radius: radius * theory.fene_force_magnitude(radius, spring_constant=2., q0=3.) \
        * theory.fene_radial_pdf(radius, **FENE)
    assert quad(weighted, 0., 3.)[0] == pytest.approx(3*.8, rel=2e-10)


def test_fene_hookean_large_extension_limit():
    result = theory.fene_mean_square_extension(spring_constant=2., q0=1e6,
                                                 thermal_energy=.8)
    assert result == pytest.approx(3*.8/2, rel=3e-12)
    assert theory.fene_potential(1., spring_constant=2., q0=1e200) == pytest.approx(1.)


def test_repulsive_cutoff_origin_and_radial_derivative():
    distance = np.array([0., .2, .8, 1., 2.])
    potential = theory.repulsive_potential(distance, epsilon=3., sigma=1.)
    force = theory.repulsive_force_magnitude(distance, epsilon=3., sigma=1.)
    np.testing.assert_allclose(potential, [2., 1.408, .112, 0., 0.])
    np.testing.assert_allclose(force, [3., 2.88, 1.08, 0., 0.])
    r = np.array([.1, .7, .99])*1.3
    delta = 1e-6
    derivative = (theory.repulsive_potential(r+delta, epsilon=4., sigma=1.3)
                  - theory.repulsive_potential(r-delta, epsilon=4., sigma=1.3))/(2*delta)
    np.testing.assert_allclose(-derivative,
        theory.repulsive_force_magnitude(r, epsilon=4., sigma=1.3), rtol=1e-9)
    assert theory.repulsive_force_magnitude(1.-1e-8, epsilon=3., sigma=1.) < 1e-6


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_conditional_pair_ideal_radial_volume(dimension):
    params = dict(epsilon=0., sigma=.7, thermal_energy=1., r_max=2., dimension=dimension)
    assert theory.repulsive_pair_radial_normalization(**params) == pytest.approx(2**dimension/dimension)
    r = np.array([0., .3, 1., 2.])
    expected = dimension*r**(dimension-1)/2**dimension
    np.testing.assert_allclose(theory.repulsive_pair_radial_pdf(r, **params), expected)
    assert theory.repulsive_pair_radial_pdf(2.1, **params) == 0.


def test_conditional_pair_normalization_and_boltzmann_ratio():
    pdf = lambda r: theory.repulsive_pair_radial_pdf(r, **PAIR)
    # Split at sigma so quadrature does not hide the cutoff.
    normalization = quad(pdf, 0., 1.3, epsabs=1e-10)[0] + quad(pdf, 1.3, 2., epsabs=1e-10)[0]
    assert normalization == pytest.approx(1., abs=2e-10)
    r1, r2 = .4, 1.6
    u1 = theory.repulsive_potential(r1, epsilon=4., sigma=1.3)
    u2 = theory.repulsive_potential(r2, epsilon=4., sigma=1.3)
    assert (pdf(r1)/r1**2)/(pdf(r2)/r2**2) == pytest.approx(math.exp(-(u1-u2)/.7))
    # Omitting the radial Jacobian gives the wrong ratio.
    assert not math.isclose(pdf(r1)/pdf(r2), math.exp(-(u1-u2)/.7), rel_tol=.01)


def test_pair_quadrature_handles_maximum_inside_force_cutoff_and_unit_change():
    params = {**PAIR, "r_max": .8}
    z = theory.repulsive_pair_radial_normalization(**params)
    direct = quad(lambda r: r*r*math.exp(-theory.repulsive_potential(r, epsilon=4., sigma=1.3)/.7),
                   0., .8, epsabs=1e-13)[0]
    assert z == pytest.approx(direct, rel=2e-11)
    # A length-unit change scales Z by length^d and a radial density by 1/length.
    scale = 7.
    changed = {**params, "sigma": params["sigma"]*scale, "r_max": params["r_max"]*scale}
    assert theory.repulsive_pair_radial_normalization(**changed) == pytest.approx(z*scale**3)
    assert theory.repulsive_pair_radial_pdf(.4*scale, **changed) == pytest.approx(
        theory.repulsive_pair_radial_pdf(.4, **params)/scale)


@pytest.mark.parametrize("override", [dict(diffusion=-1.), dict(mobility=0.),
    dict(mobility=-1.), dict(spring_constant=0.), dict(diffusion=True),
    dict(diffusion=float("nan")), dict(mobility=float("inf")), dict(spring_constant=1j)])
def test_hookean_invalid_physical_parameters(override):
    with pytest.raises(ValueError):
        theory.hookean_equilibrium_moment(**{**HOOKEAN, **override})


@pytest.mark.parametrize("moment", [[1., 2.], [[1., 2.], [0., 1.]],
    [[1., 2.], [2., 1.]], [[1., 0.], [0., float("nan")]], [], [[True]]])
def test_invalid_second_moment_rejected(moment):
    with pytest.raises(ValueError):
        theory.hookean_moment_rhs(moment, **HOOKEAN)


@pytest.mark.parametrize("call", [
    lambda: theory.brownian_msd(-1., diffusion=1.),
    lambda: theory.brownian_msd(1., diffusion=1., dimension=True),
    lambda: theory.hookean_steady_shear_moment(1., **HOOKEAN, dimension=1),
    lambda: theory.hookean_steady_shear_moment(float("nan"), **HOOKEAN),
    lambda: theory.hookean_em_moment([[1.]], True, .1, **HOOKEAN),
    lambda: theory.hookean_em_moment([[1.]], -1, .1, **HOOKEAN),
    lambda: theory.hookean_em_moment([[1.]], 1, 0., **HOOKEAN),
    lambda: theory.hookean_moment_rhs(np.eye(2), **HOOKEAN, velocity_gradient=np.eye(3)),
    lambda: theory.hookean_step_deformation(np.eye(2), np.eye(3)),
    lambda: theory.hookean_relaxation_moment(np.eye(2), [-1., 0.], **HOOKEAN),
    lambda: theory.hookean_saos_moduli(-1., **HOOKEAN),
    lambda: theory.hookean_polymer_stress(np.eye(2), **HOOKEAN, number_density=-1.),
    lambda: theory.fene_radial_pdf(-1., **FENE),
    lambda: theory.fene_mean_square_extension(**{**FENE, "thermal_energy": 0.}),
    lambda: theory.fene_parameter(**{**FENE, "q0": 0.}),
    lambda: theory.fene_parameter(**{**FENE, "spring_constant": -1.}),
    lambda: theory.fene_potential(3., spring_constant=2., q0=3.),
    lambda: theory.fene_force_magnitude(3.1, spring_constant=2., q0=3.),
    lambda: theory.repulsive_potential(-.1, epsilon=1., sigma=1.),
    lambda: theory.repulsive_force_magnitude(.1, epsilon=-1., sigma=1.),
    lambda: theory.repulsive_force_magnitude(.1, epsilon=1., sigma=0.),
    lambda: theory.repulsive_pair_radial_pdf(.1, **{**PAIR, "r_max": 0.}),
    lambda: theory.repulsive_pair_radial_pdf(.1, **{**PAIR, "thermal_energy": -1.}),
])
def test_domain_negative_controls(call):
    with pytest.raises(ValueError):
        call()


def test_zero_temperature_hookean_is_a_supported_deterministic_limit():
    params = {**HOOKEAN, "diffusion": 0.}
    np.testing.assert_array_equal(theory.hookean_equilibrium_moment(**params), np.zeros((3, 3)))
    result = theory.hookean_em_moment_step(np.eye(2), .1, **params)
    np.testing.assert_allclose(result, .64*np.eye(2))
    assert theory.hookean_saos_moduli(1., **params) == (0., 0.)


@pytest.mark.parametrize("call", [
    lambda: theory.hookean_equilibrium_moment(diffusion=1., mobility=1e308, spring_constant=1e308),
    lambda: theory.hookean_equilibrium_moment(diffusion=1., mobility=5e-324, spring_constant=.3),
    lambda: theory.hookean_steady_shear_moment(1e308, **HOOKEAN),
    lambda: theory.hookean_em_moment([[1.]], 10000, 2., **HOOKEAN),
    lambda: theory.brownian_msd(1e308, diffusion=1e308),
    lambda: theory.fene_parameter(spring_constant=1e308, q0=1e308, thermal_energy=1.),
    lambda: theory.fene_potential(1e200, spring_constant=2., q0=2e200),
])
def test_nonfinite_derived_quantities_rejected(call):
    with pytest.raises(ValueError):
        call()


def test_output_ownership_does_not_modify_inputs():
    moment = np.eye(2)
    deformation = np.array([[1., .3], [0., 1.]])
    moment_before, deformation_before = moment.copy(), deformation.copy()
    output = theory.hookean_step_deformation(moment, deformation)
    output[0, 0] = 123.
    np.testing.assert_array_equal(moment, moment_before)
    np.testing.assert_array_equal(deformation, deformation_before)
