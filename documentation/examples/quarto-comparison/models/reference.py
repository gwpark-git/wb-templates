"""Direct polynomial reference for the synthetic documentation example."""


def properties(phi):
    """Return reduced diffusivity, relative viscosity and osmotic factor."""
    diffusion = 1.0 + 1.2 * phi
    viscosity = 1.0 + 2.5 * phi + 2.0 * phi ** 2
    osmotic_factor = 1.0 + 4.0 * phi + 10.0 * phi ** 2
    return diffusion, viscosity, osmotic_factor
