"""Typed, independent reference models for particle BD and kinetic dumbbells.

These engines use explicit mobility/diffusion and a Gaussian noise default.
They supplement the original core; they do not change its clocks/checkpoints.
See the specification metadata for units and stress conventions.
"""
from ._reference import (DumbbellSimulation, DumbbellSpec,
                         RepulsiveSimulation, RepulsiveSpec, execution_backends)

__all__ = ["DumbbellSimulation", "DumbbellSpec",
           "RepulsiveSimulation", "RepulsiveSpec", "execution_backends"]
