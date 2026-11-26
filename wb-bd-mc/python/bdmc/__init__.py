"""BD–MC reference components; ML remains an optional, separate layer."""

__version__ = "0.1.0"


def __getattr__(name):
    # Allow legacy readers and ML tools to work without loading the C++ extension.
    if name in {"Simulation", "SimulationSpec", "neighbor_pairs"}:
        from . import _core
        return getattr(_core, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

