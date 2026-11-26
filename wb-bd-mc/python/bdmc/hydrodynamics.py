"""Free-space dense native HI references, separate from periodic simulation.

The Brownian engine uses GSL Gaussian noise and owns its state.  The Stokesian
operator is instantaneous and athermal.  Both return copied NumPy arrays; their
metadata records units, boundaries and the scientific limitations explicitly.
"""
from ._hydrodynamics import (
    HydrodynamicSpec,
    HydrodynamicSimulation,
    StokesianSpec,
    mobility,
    solve_stokesian,
)

__all__ = [
    "HydrodynamicSpec", "HydrodynamicSimulation", "StokesianSpec",
    "mobility", "solve_stokesian",
]
