"""Independent analytic oracles for the two new free-draining references.

These functions reconstruct standard Hookean/FENE dumbbell theory and the
declared soft repulsive potential; they do not call the native implementation.
``diffusion`` is a single bead's D, ``mobility`` is a single bead's mu, and
``spring_constant`` is H in U(q)=H*q**2/2.  Thermal consistency means kT=D/mu.
See docs/TWO_REFERENCE_THEORY.org for equations, units and source limitations.

NumPy is required.  General constant-flow propagation, FENE densities and pair
quadrature import SciPy only when called.  No stochastic simulation is performed.
"""

from __future__ import annotations

import math
from numbers import Integral, Real

import numpy as np

__all__ = [
    "thermal_energy", "brownian_msd", "hookean_relaxation_time",
    "hookean_equilibrium_moment", "hookean_moment_rhs",
    "hookean_em_moment_step", "hookean_em_moment",
    "hookean_em_equilibrium_moment", "hookean_moment_evolution",
    "hookean_steady_shear_moment", "hookean_step_deformation",
    "hookean_relaxation_moment", "hookean_polymer_stress",
    "hookean_saos_moduli", "fene_parameter", "fene_potential",
    "fene_force_magnitude", "fene_radial_pdf", "fene_mean_square_extension",
    "repulsive_potential", "repulsive_force_magnitude",
    "repulsive_pair_radial_normalization", "repulsive_pair_radial_pdf",
]


def _scalar(value, name: str, *, positive: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real scalar") from exc
    if not math.isfinite(result) or (result <= 0 if positive else result < 0):
        condition = "positive" if positive else "nonnegative"
        raise ValueError(f"{name} must be finite and {condition}")
    return result


def _signed_scalar(value, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real scalar") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real scalar")
    return result


def _integer(value, name: str, *, minimum: int = 1) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    if value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _array(value, name: str, *, nonnegative: bool = False) -> np.ndarray:
    try:
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf":
            raise ValueError(f"{name} must contain real numbers")
        result = raw.astype(float, copy=False)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must contain finite real numbers") from exc
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite real numbers")
    if nonnegative and np.any(result < 0):
        raise ValueError(f"{name} must be nonnegative")
    return result


def _finite(value, name: str):
    if not np.all(np.isfinite(value)):
        raise ValueError(f"derived {name} is nonfinite")
    return value


def _scalar_or_array(value):
    return float(value) if np.ndim(value) == 0 else value


def _moment(value) -> np.ndarray:
    result = _array(value, "moment")
    if result.ndim != 2 or not result.shape[0] or result.shape[0] != result.shape[1]:
        raise ValueError("moment must be a nonempty square matrix")
    scale = max(1.0, float(np.max(np.abs(result))))
    scaled = result / scale
    if not np.allclose(scaled, scaled.T, rtol=1e-12, atol=1e-14):
        raise ValueError("moment must be symmetric")
    if np.min(np.linalg.eigvalsh((scaled + scaled.T) / 2)) < -1e-12:
        raise ValueError("moment must be positive semidefinite")
    return result.copy()


def _gradient(value, dimension: int) -> np.ndarray:
    if value is None:
        return np.zeros((dimension, dimension))
    result = _array(value, "velocity_gradient")
    if result.shape != (dimension, dimension):
        raise ValueError("velocity_gradient must have the same shape as moment")
    return result


def _hookean(diffusion, mobility, spring_constant) -> tuple[float, float, float]:
    d = _scalar(diffusion, "diffusion")
    mu = _scalar(mobility, "mobility", positive=True)
    h = _scalar(spring_constant, "spring_constant", positive=True)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        mu_h = float(np.float64(mu) * h)
        rate = float(np.float64(4) * mu_h)
    if not math.isfinite(rate) or rate == 0 or mu_h == 0:
        raise ValueError("derived moment relaxation rate must be finite and positive")
    c = d / mu_h
    return d, rate, float(_finite(c, "equilibrium component"))


def thermal_energy(*, diffusion: float, mobility: float) -> float:
    """Return kT=D/mu, imposing the free-draining Einstein relation."""
    d = _scalar(diffusion, "diffusion")
    mu = _scalar(mobility, "mobility", positive=True)
    return float(_finite(d / mu, "thermal_energy"))


def brownian_msd(time, *, diffusion: float, dimension: int = 3,
                 relative: bool = False):
    """Unwrapped ideal MSD: 2*d*D*t, or 4*d*D*t for two-bead separation."""
    t = _array(time, "time", nonnegative=True)
    d = _scalar(diffusion, "diffusion")
    dim = _integer(dimension, "dimension")
    if not isinstance(relative, (bool, np.bool_)):
        raise ValueError("relative must be boolean")
    with np.errstate(over="ignore", invalid="ignore"):
        result = (4 if relative else 2) * dim * d * t
    return _scalar_or_array(_finite(result, "MSD"))


def hookean_relaxation_time(*, mobility: float, spring_constant: float) -> float:
    """Second-moment/stress relaxation time lambda=1/(4*mu*H)."""
    _, rate, _ = _hookean(0.0, mobility, spring_constant)
    return float(_finite(1 / rate, "relaxation time"))


def hookean_equilibrium_moment(*, diffusion: float, mobility: float,
                              spring_constant: float, dimension: int = 3):
    """Return E[q q.T]=D/(mu*H)*I for the continuum Hookean model."""
    _, _, c = _hookean(diffusion, mobility, spring_constant)
    return c * np.eye(_integer(dimension, "dimension"))


def hookean_moment_rhs(moment, *, diffusion: float, mobility: float,
                      spring_constant: float, velocity_gradient=None):
    """Continuum dM/dt=K M+M K.T-4*mu*H*M+4*D*I."""
    m = _moment(moment)
    d, rate, _ = _hookean(diffusion, mobility, spring_constant)
    k = _gradient(velocity_gradient, len(m))
    with np.errstate(over="ignore", invalid="ignore"):
        rhs = k @ m + m @ k.T - rate * m + (4 * d) * np.eye(len(m))
    return _finite(rhs, "moment derivative")


def _euler_transition(dimension, dt, diffusion, mobility, spring_constant,
                      velocity_gradient):
    dt = _scalar(dt, "dt", positive=True)
    d, rate, _ = _hookean(diffusion, mobility, spring_constant)
    k = _gradient(velocity_gradient, dimension)
    with np.errstate(over="ignore", invalid="ignore"):
        a = np.eye(dimension) + dt * (k - rate / 2 * np.eye(dimension))
        noise = (4 * d * dt) * np.eye(dimension)
    return _finite(a, "Euler transition"), _finite(noise, "Euler noise covariance")


def hookean_em_moment_step(moment, dt: float, *, diffusion: float, mobility: float,
                           spring_constant: float, velocity_gradient=None):
    """Exact one-step second moment of simultaneous-bead Euler--Maruyama.

    The formula remains an algebraic oracle for unstable dt; it does not claim
    that every accepted time step has an equilibrium distribution.
    """
    m = _moment(moment)
    a, noise = _euler_transition(len(m), dt, diffusion, mobility, spring_constant,
                                 velocity_gradient)
    with np.errstate(over="ignore", invalid="ignore"):
        result = a @ m @ a.T + noise
    return _finite(result, "Euler moment")


def hookean_em_moment(moment, steps: int, dt: float, *, diffusion: float,
                     mobility: float, spring_constant: float,
                     velocity_gradient=None):
    """Exact n-step Euler second moment, using a lifted affine matrix power.

    Work grows logarithmically with ``steps``.  This is the finite-step oracle,
    distinct from ``hookean_moment_evolution``'s continuum prediction.
    """
    m = _moment(moment)
    n = _integer(steps, "steps", minimum=0)
    a, noise = _euler_transition(len(m), dt, diffusion, mobility, spring_constant,
                                 velocity_gradient)
    size = m.size
    lifted = np.zeros((size + 1, size + 1))
    with np.errstate(over="ignore", invalid="ignore"):
        lifted[:size, :size] = np.kron(a, a)
        lifted[:size, size] = noise.ravel()
        lifted[size, size] = 1.0
        initial = np.r_[m.ravel(), 1.0]
        result = (np.linalg.matrix_power(lifted, n) @ initial)[:size]
    return _finite(result.reshape(m.shape), "Euler moment")


def hookean_em_equilibrium_moment(dt: float, *, diffusion: float, mobility: float,
                                 spring_constant: float, dimension: int = 3):
    """No-flow Euler stationary moment; requires 0 < mu*H*dt < 1."""
    dt = _scalar(dt, "dt", positive=True)
    _, rate, c = _hookean(diffusion, mobility, spring_constant)
    x = rate / 4 * dt
    if not 0 < x < 1:
        raise ValueError("Euler equilibrium requires 0 < mobility*spring_constant*dt < 1")
    return _finite((c / (1 - x)) * np.eye(_integer(dimension, "dimension")),
                   "Euler equilibrium moment")


def hookean_moment_evolution(moment, time, *, diffusion: float, mobility: float,
                            spring_constant: float, velocity_gradient=None):
    """Exact continuum second moment at nonnegative time under constant K.

    Scalar time returns (d,d); a time array returns time.shape+(d,d).
    A lifted affine generator handles both invertible and singular generators
    without subtracting two large stationary/initial contributions.
    """
    from scipy.linalg import expm

    m = _moment(moment)
    t = _array(time, "time", nonnegative=True)
    d, rate, _ = _hookean(diffusion, mobility, spring_constant)
    k = _gradient(velocity_gradient, len(m))
    b = k - rate / 2 * np.eye(len(m))
    size = m.size
    generator = np.zeros((size + 1, size + 1))
    generator[:size, :size] = np.kron(b, np.eye(len(m))) + np.kron(np.eye(len(m)), b)
    generator[:size, size] = (4 * d * np.eye(len(m))).ravel()
    _finite(generator, "moment generator")
    initial = np.r_[m.ravel(), 1.0]
    output = np.empty(t.shape + m.shape)
    with np.errstate(over="ignore", invalid="ignore"):
        for index in np.ndindex(t.shape):
            argument = _finite(generator * t[index], "moment propagator argument")
            result = (expm(argument) @ initial)[:size].reshape(m.shape)
            output[index] = _finite(result / 2 + result.T / 2, "continuum moment")
    return output


def hookean_steady_shear_moment(rate: float, *, diffusion: float, mobility: float,
                              spring_constant: float, dimension: int = 3):
    """Steady simple shear: K[x,y]=rate, with x=0 and y=1."""
    shear = _signed_scalar(rate, "rate")
    dim = _integer(dimension, "dimension", minimum=2)
    _, decay, c = _hookean(diffusion, mobility, spring_constant)
    result = c * np.eye(dim)
    if c == 0:
        return result
    with np.errstate(over="ignore", invalid="ignore"):
        wi = np.float64(shear) / decay
        result[0, 1] = result[1, 0] = c * wi
        result[0, 0] = c * (1 + 2 * wi**2)
    return _finite(result, "steady shear moment")


def hookean_step_deformation(moment, deformation):
    """Apply an instantaneous affine step: M(0+)=F M(0-) F.T."""
    m = _moment(moment)
    f = _array(deformation, "deformation")
    if f.shape != m.shape:
        raise ValueError("deformation must have the same shape as moment")
    with np.errstate(over="ignore", invalid="ignore"):
        result = f @ m @ f.T
    return _finite(result, "deformed moment")


def hookean_relaxation_moment(moment, time, *, diffusion: float, mobility: float,
                            spring_constant: float):
    """No-flow continuum relaxation, c I+(M0-c I)*exp(-t/lambda)."""
    m = _moment(moment)
    t = _array(time, "time", nonnegative=True)
    _, rate, c = _hookean(diffusion, mobility, spring_constant)
    equilibrium = c * np.eye(len(m))
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        result = equilibrium + np.exp(-rate * t)[..., None, None] * (m - equilibrium)
    return _finite(result, "relaxation moment")


def hookean_polymer_stress(moment, *, diffusion: float, mobility: float,
                           spring_constant: float, number_density: float = 1.0):
    """Polymer Kramers stress n*(H*M-kT*I), with tensile-positive sign.

    ``number_density`` counts dumbbells per volume.  Solvent stress is excluded.
    """
    m = _moment(moment)
    n = _scalar(number_density, "number_density")
    _hookean(diffusion, mobility, spring_constant)
    kt = thermal_energy(diffusion=diffusion, mobility=mobility)
    with np.errstate(over="ignore", invalid="ignore"):
        result = n * (spring_constant * m - kt * np.eye(len(m)))
    return _finite(result, "polymer stress")


def hookean_saos_moduli(angular_frequency, *, diffusion: float, mobility: float,
                        spring_constant: float, number_density: float = 1.0):
    """Return (G', G'') for the single-mode Hookean polymer contribution.

    gamma=gamma0*sin(omega*t) gives tau_xy/gamma0=G'*sin+G''*cos.
    The amplitude must be small for the stated SAOS interpretation.
    """
    omega = _array(angular_frequency, "angular_frequency", nonnegative=True)
    n = _scalar(number_density, "number_density")
    relaxation = hookean_relaxation_time(mobility=mobility,
                                         spring_constant=spring_constant)
    g = _finite(n * thermal_energy(diffusion=diffusion, mobility=mobility), "modulus")
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        x = omega * relaxation
        # hypot avoids overflow in 1+x*x; x=inf is the exact high-frequency limit.
        denominator = np.hypot(1.0, x)
        a = np.divide(x, denominator, out=np.ones_like(x), where=np.isfinite(x))
        inverse = 1 / denominator
        storage = g * a**2
        loss = g * a * inverse
    return (_scalar_or_array(_finite(storage, "storage modulus")),
            _scalar_or_array(_finite(loss, "loss modulus")))


def fene_parameter(*, spring_constant: float, q0: float,
                   thermal_energy: float) -> float:
    """True FENE extensibility b=H*Q0**2/kT (not a FENE-P closure)."""
    h = _scalar(spring_constant, "spring_constant", positive=True)
    q0 = _scalar(q0, "q0", positive=True)
    kt = _scalar(thermal_energy, "thermal_energy", positive=True)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        b = float(np.float64(h) * q0 * q0 / kt)
    if not math.isfinite(b) or b <= 0:
        raise ValueError("derived FENE parameter must be finite and positive")
    return b


def _fene_extension(extension, spring_constant, q0):
    q = _array(extension, "extension", nonnegative=True)
    h = _scalar(spring_constant, "spring_constant", positive=True)
    q0 = _scalar(q0, "q0", positive=True)
    if np.any(q >= q0):
        raise ValueError("FENE extension must be strictly less than q0")
    return q, h, q0


def fene_potential(extension, *, spring_constant: float, q0: float):
    """U(q)=-H*Q0**2/2*log(1-q**2/Q0**2), for 0 <= q < Q0."""
    q, h, q0 = _fene_extension(extension, spring_constant, q0)
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        x_squared = (q / q0)**2
        # Rewrite with q^2 to preserve the Hookean limit when q/Q0 is tiny;
        # the logarithmic ratio tends to one even if x_squared underflows.
        ratio = np.divide(-np.log1p(-x_squared), x_squared,
                          out=np.ones_like(q), where=x_squared > 0)
        result = 0.5 * h * q * q * ratio
    return _scalar_or_array(_finite(result, "FENE potential"))


def fene_force_magnitude(extension, *, spring_constant: float, q0: float):
    """Restoring-force magnitude H*q/(1-q**2/Q0**2), for 0 <= q < Q0."""
    q, h, q0 = _fene_extension(extension, spring_constant, q0)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        result = h * q / (1 - (q / q0)**2)
    return _scalar_or_array(_finite(result, "FENE force"))


def fene_radial_pdf(extension, *, spring_constant: float, q0: float,
                    thermal_energy: float, dimension: int = 3):
    """Normalized radial density, including the d-dimensional radial Jacobian.

    Density is zero at and above Q0.  Negative extension is rejected.  This is
    true FENE equilibrium, not the FENE-P preaveraged spring approximation.
    """
    from scipy.special import betaln

    q = _array(extension, "extension", nonnegative=True)
    dim = _integer(dimension, "dimension")
    b = fene_parameter(spring_constant=spring_constant, q0=q0,
                        thermal_energy=thermal_energy)
    q0 = float(q0)
    result = np.zeros_like(q)
    inside = (q > 0) & (q < q0)
    x = q[inside] / q0
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        log_pdf = (math.log(2) - math.log(q0) - betaln(dim / 2, b / 2 + 1)
                   + (dim - 1) * np.log(x) + b / 2 * np.log1p(-x * x))
        result[inside] = np.exp(log_pdf)
        if dim == 1:
            result[q == 0] = np.exp(math.log(2) - math.log(q0) - betaln(.5, b / 2 + 1))
    return _scalar_or_array(_finite(result, "FENE radial density"))


def fene_mean_square_extension(*, spring_constant: float, q0: float,
                               thermal_energy: float, dimension: int = 3) -> float:
    """Exact true-FENE equilibrium E[q**2]=d*Q0**2/(b+d+2)."""
    dim = _integer(dimension, "dimension")
    b = fene_parameter(spring_constant=spring_constant, q0=q0,
                        thermal_energy=thermal_energy)
    with np.errstate(over="ignore", invalid="ignore"):
        result = float(np.float64(dim) / (b + dim + 2) * q0 * q0)
    return float(_finite(result, "FENE mean square extension"))


def repulsive_potential(distance, *, epsilon: float, sigma: float):
    """Soft finite-range U=epsilon*(1-r/sigma)**2*(2+r/sigma)/3 for r<sigma."""
    r = _array(distance, "distance", nonnegative=True)
    eps = _scalar(epsilon, "epsilon")
    cutoff = _scalar(sigma, "sigma", positive=True)
    result = np.zeros_like(r)
    inside = r < cutoff
    x = r[inside] / cutoff
    with np.errstate(over="ignore", invalid="ignore"):
        result[inside] = eps * ((1 - x)**2 * (2 + x) / 3)
    return _scalar_or_array(_finite(result, "repulsive potential"))


def repulsive_force_magnitude(distance, *, epsilon: float, sigma: float):
    """Outward magnitude -dU/dr=epsilon*(1-(r/sigma)**2)/sigma within cutoff.

    At r=0 this reports the radial limit; a vector direction is undefined.
    The native pair solver is responsible for its explicit coincidence policy.
    """
    r = _array(distance, "distance", nonnegative=True)
    eps = _scalar(epsilon, "epsilon")
    cutoff = _scalar(sigma, "sigma", positive=True)
    result = np.zeros_like(r)
    inside = r < cutoff
    x = r[inside] / cutoff
    with np.errstate(over="ignore", invalid="ignore"):
        result[inside] = (eps / cutoff) * (1 - x * x)
    return _scalar_or_array(_finite(result, "repulsive force"))


def _pair_parameters(epsilon, sigma, thermal_energy, r_max, dimension):
    return (_scalar(epsilon, "epsilon"), _scalar(sigma, "sigma", positive=True),
            _scalar(thermal_energy, "thermal_energy", positive=True),
            _scalar(r_max, "r_max", positive=True), _integer(dimension, "dimension"))


def repulsive_pair_radial_normalization(*, epsilon: float, sigma: float,
                                        thermal_energy: float, r_max: float,
                                        dimension: int = 3) -> float:
    """Integral Z=integral_0^rmax r**(d-1)*exp(-U(r)/kT) dr.

    In a periodic cubic box, choosing r_max <= L/2 makes this a conditional
    spherical-shell oracle; it is not the complete cube's radial distribution.
    The angular surface measure cancels from the normalized density.
    """
    from scipy.integrate import quad

    eps, cutoff, kt, radius, dim = _pair_parameters(
        epsilon, sigma, thermal_energy, r_max, dimension)
    # Integrate a dimensionless radius to keep absolute quadrature tolerances
    # independent of the chosen length unit. Split at the force cutoff.
    split = min(1.0, cutoff / radius)

    def integrand(x):
        u = repulsive_potential(x * radius, epsilon=eps, sigma=cutoff)
        return x**(dim - 1) * math.exp(-u / kt)

    inside, _ = quad(integrand, 0.0, split, epsabs=1e-12, epsrel=1e-11)
    outside = (1 - split**dim) / dim
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        result = float(np.float64(radius)**dim * (inside + outside))
    if not math.isfinite(result) or result <= 0:
        raise ValueError("derived pair radial normalization must be finite and positive")
    return result


def repulsive_pair_radial_pdf(distance, *, epsilon: float, sigma: float,
                              thermal_energy: float, r_max: float,
                              dimension: int = 3):
    """Conditional equilibrium p(r)=r**(d-1)*exp(-U/kT)/Z on [0,r_max].

    This isolated-pair result does not predict the interacting many-particle
    radial distribution. Values beyond r_max have zero density.
    """
    r = _array(distance, "distance", nonnegative=True)
    eps, cutoff, kt, radius, dim = _pair_parameters(
        epsilon, sigma, thermal_energy, r_max, dimension)
    z = repulsive_pair_radial_normalization(epsilon=eps, sigma=cutoff,
        thermal_energy=kt, r_max=radius, dimension=dim)
    result = np.zeros_like(r)
    inside = r <= radius
    u = repulsive_potential(r[inside], epsilon=eps, sigma=cutoff)
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        result[inside] = r[inside]**(dim - 1) * np.exp(-u / kt) / z
    return _scalar_or_array(_finite(result, "pair radial density"))
