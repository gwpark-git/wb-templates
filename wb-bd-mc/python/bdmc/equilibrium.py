"""Reduced equilibrium pressure observables for the independent EOS workflows."""

import math

import numpy as np

HARD_CONTACT_CAP = 1.05
HARD_SHELL_WIDTHS = (0.0025, 0.005, 0.01, 0.02, 0.04)


def carnahan_starling(phi):
    """Hard-sphere compressibility factor at nominal volume fraction phi."""
    phi = np.asarray(phi, dtype=float)
    if not np.isfinite(phi).all() or np.any((phi < 0) | (phi >= 1)):
        raise ValueError("volume fractions must be finite and in [0, 1)")
    return (1 + phi + phi**2 - phi**3) / (1 - phi)**3


def hard_sphere_survival(minimum_r2, *, particles, compression=0.25):
    """Indicators for a virtual compression with x = N log(V/V')."""
    if not isinstance(particles, int) or particles < 1:
        raise ValueError("positive particle count required")
    log_threshold = 2 * (compression / particles) / 3
    if not 0 < log_threshold <= 2 * math.log(HARD_CONTACT_CAP):
        raise ValueError("compression must be positive and within recorded distance coverage")
    radii = np.asarray(minimum_r2, dtype=float)
    if (radii.ndim != 1 or not radii.size or not np.isfinite(radii).all() or
            np.any((radii < 1) | (radii > HARD_CONTACT_CAP**2))):
        raise ValueError("expected nonoverlapping squared distances capped at 1.05 squared")
    return radii >= math.exp(log_threshold)


def virtual_compression_z(survival_fraction, compression=0.25):
    """Finite log-volume secant of the hard-sphere compressibility factor."""
    if not math.isfinite(compression) or not compression > 0:
        raise ValueError("positive finite compression required")
    probability = np.asarray(survival_fraction, dtype=float)
    if not np.isfinite(probability).all() or np.any((probability <= 0) | (probability > 1)):
        raise ValueError("finite nonzero survival probability required")
    return 1 - np.log(probability) / compression


def contact_shell_z(mean_unordered_pairs, *, particles, width):
    """Finite shell average with canonical N² normalization."""
    count = np.asarray(mean_unordered_pairs, dtype=float)
    if (not isinstance(particles, int) or particles < 1 or
            not 0 < width <= HARD_CONTACT_CAP - 1 or
            not np.isfinite(count).all() or np.any((count < 0) | (count > particles * (particles - 1) / 2))):
        raise ValueError("invalid hard-sphere shell population or width")
    return 1 + count / (particles * ((1 + width)**3 - 1))
