"""Factored polynomial implementation for the synthetic documentation example."""


def properties(phi):
    """Return the same reduced quantities with a factored evaluation order."""
    diffusion = 1.0 + 1.2 * phi
    viscosity = 1.0 + phi * (2.5 + 2.0 * phi)
    osmotic_factor = 1.0 + phi * (4.0 + 10.0 * phi)
    return diffusion, viscosity, osmotic_factor
