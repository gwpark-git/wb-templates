"""Explicitly normalized structural and temporal observers.

The native observers accept snapshots or fixed-N batches, never alter the
simulator and return independent readonly NumPy arrays. Directional results are
3D solid-angle densities; the planar map is a projection of 3D pair density.
"""

from __future__ import annotations

import numpy as np

from ._analysis import (
    autocorrelation,
    directional_pair_correlation,
    pair_correlation,
    planar_pair_correlation,
)


def coordination_number(correlation, cutoff=None, *, source="counts", allow_truncated=False):
    """Return mean distinct neighbors per particle up to radial bin edges.

    ``source='counts'`` uses the empirical ordered-pair counts. ``source='g'``
    integrates the piecewise constant RDF with exact shell volumes and its
    declared normalization density, including the finite-N correction.
    A scalar cutoff must be a stored edge; interpolation inside a histogram bin
    would invent information. A histogram starting above zero needs an explicit
    ``allow_truncated=True`` and then reports the shell integral from that edge.
    """
    if correlation.get("kind") != "radial":
        raise ValueError("coordination_number requires an isotropic radial result")
    edges = np.asarray(correlation["radial_edges"], dtype=float)
    if edges[0] != 0 and not allow_truncated:
        raise ValueError("radial histogram starts above zero; total coordination is unavailable")
    if source == "counts":
        values = np.cumsum(np.asarray(correlation["counts"], dtype=float)) / (
            correlation["n_frames"] * correlation["n_particles"]
        )
    elif source == "g":
        # g and its normalization density form a pair. Substituting N/V here
        # for an RDF normalized with (N-1)/V changes the finite-N neighbor count.
        values = np.cumsum(
            np.asarray(correlation["g"]) * correlation["bin_volumes"]
            * correlation["normalization_density"]
        )
    else:
        raise ValueError("source must be 'counts' or 'g'")
    if cutoff is None:
        values.setflags(write=False)
        return values
    if not np.isscalar(cutoff) or not np.isfinite(cutoff):
        raise ValueError("cutoff must be a finite scalar radial edge")
    matches = np.flatnonzero(edges == cutoff)
    if not len(matches):
        raise ValueError("cutoff must equal a stored radial edge; re-bin positions for another cutoff")
    index = int(matches[0])
    return 0.0 if index == 0 else float(values[index - 1])


__all__ = [
    "pair_correlation", "directional_pair_correlation", "planar_pair_correlation",
    "coordination_number", "autocorrelation",
]
