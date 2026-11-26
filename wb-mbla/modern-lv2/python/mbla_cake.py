"""One mBLA cake-resistance search on an already computed radial slice.

The outer fixed-point solver supplies the slice and applies under-relaxation.
This module preserves the inner search and offers a compiled Python backend.
"""

from dataclasses import dataclass
from functools import lru_cache
import math
from typing import Sequence


@dataclass(frozen=True)
class CakeSearchState:
    """Radial fields and scalars at the start of one resistance search."""

    wall_distance: Sequence[float]  # y/R, increasing from zero at the wall
    viscosity_integral: Sequence[float]  # normalized at the centerline
    inverse_diffusivity_integral: Sequence[float]  # integrated in y/R
    uz_ratio: float
    permeability: float  # m/(Pa s), clean membrane
    net_driving_pressure: float  # Pa, P - Pper - Pi at freeze
    reference_vw: float  # m/s
    epsilon_d: float  # D0/(R*reference_vw)
    area_factor: float  # (1 - cake thickness/R)^2 from the previous iteration
    inlet_flux_factor: float  # phi_bulk*(F2_0 - F2_modified_0)
    bulk_fraction: float
    freeze_fraction: float
    resistance_tolerance: float  # dimensionless Rcake/Rmem step tolerance


@dataclass(frozen=True)
class CakeSearchResult:
    """Carry resistance before relaxation, termination and the optional residual trace."""
    resistance: float  # Rcake/Rmem before the caller's under-relaxation
    steps: int
    termination: str
    trial_trace: tuple[tuple[float, float, float], ...]  # trial, step and normalized residual


def _search(y, i_eta, i_diff, uz, lp, net_pressure, reference_vw, epsilon, area,
            inlet_flux, bulk, freeze, tolerance, trace):
    # Geometry and F2 are independent of the trial resistance. Retain each
    # trapezoid sum's original increasing-y accumulation order.
    """Follow inherited trial order, separate radial sums and bounded refinement."""
    geometry = [(1. - y[j]) * (2. - y[j]) * i_eta[j] for j in range(len(y))]
    half_spacing = [0.5 * (y[j] - y[j - 1]) for j in range(1, len(y))]
    sum_f2 = 0.
    for j in range(1, len(y)):
        sum_f2 += half_spacing[j - 1] * (geometry[j - 1] + geometry[j])
    f2 = area * (sum_f2 * uz)

    trial = 0.
    step = 0.01
    # A negative initial residual brackets the first positive trial.
    residual = -1.
    history = []
    for iteration in range(1, 200001):
        previous_trial = trial
        previous_residual = residual
        trial += step
        vw = lp * (1. / (1. + trial)) * net_pressure / reference_vw
        alpha = vw / epsilon

        # F1 and modified F2 share their attenuation at each radial node.
        # Evaluate exp once and keep independent sums in source order.
        exponent = alpha * i_diff[0]
        prior_f1 = geometry[0] * math.exp(-exponent)
        prior_modified = prior_f1 * exponent
        sum_f1 = 0.
        sum_modified = 0.
        for j in range(1, len(y)):
            exponent = alpha * i_diff[j]
            local_f1 = geometry[j] * math.exp(-exponent)
            local_modified = local_f1 * exponent
            half_dy = half_spacing[j - 1]
            sum_f1 += half_dy * (prior_f1 + local_f1)
            sum_modified += half_dy * (prior_modified + local_modified)
            prior_f1 = local_f1
            prior_modified = local_modified

        f1 = area * (sum_f1 * uz)
        f2_modified = area * (sum_modified * uz)
        residual = ((freeze - bulk) * f1 + bulk * (f2 - f2_modified) - inlet_flux) / inlet_flux
        if not math.isfinite(residual):
            raise ArithmeticError("nonfinite cake residual")
        if trace:
            history.append((trial, step, residual))

        # The trend/range branch runs first. The second branch uses the
        # possibly restored residual and step from that decision.
        if trial > 100. or (residual != -1. and residual - previous_residual < 0.):
            if step <= tolerance:
                return previous_trial, iteration, True, history
            step /= 10.
            trial = 0.
            residual = previous_residual
        if residual * previous_residual < 0.:
            if step <= tolerance:
                return previous_trial, iteration, False, history
            step /= 10.
            trial = previous_trial
            residual = previous_residual
    raise ArithmeticError("cake resistance search did not terminate")


@lru_cache(maxsize=1)
def _compiled_backend():
    """Import optional dependencies only when the compiled backend is selected."""
    import numpy
    import numba

    return numpy, numba.njit(fastmath=False, cache=False)(_search)


def search_resistance(state: CakeSearchState, *, backend: str = "python",
                      trace: bool = False) -> CakeSearchResult:
    """Run a bounded search on explicit in-memory data, with or without Numba."""
    if backend not in ("python", "numba"):
        raise ValueError("backend must be 'python' or 'numba'")
    arrays = tuple(tuple(map(float, values)) for values in
                   (state.wall_distance, state.viscosity_integral, state.inverse_diffusivity_integral))
    radial = arrays[0]
    if not 3 <= len(radial) <= 4000 or any(len(array) != len(radial) for array in arrays):
        raise ValueError("three matched radial arrays with 3 to 4000 samples are required")
    if not all(math.isfinite(value) for array in arrays for value in array):
        raise ValueError("radial arrays must be finite")
    if radial[0] != 0. or radial[-1] > 1. or any(a >= b for a, b in zip(radial, radial[1:])):
        raise ValueError("wall distances must increase from zero to at most one")
    scalars = (
        state.uz_ratio, state.permeability, state.net_driving_pressure, state.reference_vw,
        state.epsilon_d, state.area_factor, state.inlet_flux_factor, state.bulk_fraction,
        state.freeze_fraction, state.resistance_tolerance,
    )
    if not all(map(math.isfinite, scalars)):
        raise ValueError("search parameters must be finite")
    if (min(state.permeability, state.reference_vw, state.epsilon_d, state.resistance_tolerance) <= 0.
            or state.area_factor < 0. or state.inlet_flux_factor == 0.
            or not 0. < state.bulk_fraction < state.freeze_fraction < 1.):
        raise ValueError("invalid scales, concentrations or zero inlet flux factor")

    if backend == "numba":
        numpy, search = _compiled_backend()
        inputs = tuple(numpy.asarray(array, dtype=numpy.float64) for array in arrays)
    else:
        search = _search
        inputs = arrays
    try:
        resistance, steps, warning, history = search(*inputs, *scalars, trace)
    except (OverflowError, ZeroDivisionError) as error:
        raise ArithmeticError("nonfinite cake residual") from error
    return CakeSearchResult(resistance, steps, "warning" if warning else "bracketed", tuple(history))
