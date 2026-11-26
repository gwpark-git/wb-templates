"""Appearance-only researcher presets; scientific transformations stay outside."""
from __future__ import annotations

from contextlib import contextmanager
from importlib import metadata, resources
import json
import math

import matplotlib as mpl
from matplotlib import colors, font_manager
from cycler import cycler
import numpy as np

__version__ = "0.2.0"
PRESETS = ("personal",)
COLOR_MODES = ("color", "monochrome")


def resolve_style(name):
    """The personal style has presentation modes, not author/year editions."""
    if not isinstance(name, str) or name not in PRESETS:
        raise ValueError(f"unknown style {name!r}; choose {', '.join(PRESETS)}")
    return name

# Friendly names and allowlists keep style overrides inside the appearance API.
_ALIASES = {"figsize": ("figure.figsize",), "dpi": ("figure.dpi", "savefig.dpi"),
            "linewidth": ("lines.linewidth",), "markersize": ("lines.markersize",),
            "labelsize": ("axes.labelsize",),
            "tick_labelsize": ("xtick.labelsize", "ytick.labelsize"),
            "bbox_inches": ("savefig.bbox",)}
_NUMERIC_RC = {"lines.linewidth", "lines.markersize", "axes.labelsize", "xtick.labelsize",
               "ytick.labelsize", "legend.fontsize", "font.size", "axes.linewidth",
               "xtick.major.size", "ytick.major.size", "xtick.major.width", "ytick.major.width",
               "xtick.minor.size", "ytick.minor.size", "xtick.minor.width", "ytick.minor.width"}
_COLOR_RC = {"lines.color", "axes.edgecolor", "axes.labelcolor", "text.color", "xtick.color", "ytick.color"}
_BASELINE_KEYS = {"xtick.direction", "ytick.direction", "xtick.major.size", "ytick.major.size",
                  "xtick.major.width", "ytick.major.width", "xtick.minor.size", "ytick.minor.size",
                  "xtick.minor.width", "ytick.minor.width", "axes.linewidth",
                  "legend.fontsize", "legend.frameon", "font.size", "axes.prop_cycle",
                  "lines.markerfacecolor", "lines.markeredgewidth"}


def _preset(name):
    """Read the strict JSON resource as a fresh nested mapping for one supported style."""
    return json.loads(resources.files(__package__).joinpath("presets.json").read_text())[resolve_style(name)]


def plot_defaults(name="personal") -> dict:
    """Return generic presentation kwargs only; no x/y limits or physics."""
    return _preset(name)["defaults"]


def dash_pattern(pattern, preset="personal"):
    """Return a zero-offset Matplotlib dash tuple with point lengths from presets.json."""
    config = _preset(preset)
    if "dash_patterns" not in config or pattern not in config["dash_patterns"]:
        raise ValueError(f"no source custom dash {pattern!r} for preset {preset!r}")
    return (0, tuple(config["dash_patterns"][pattern]))


def batlow_colors(positions):
    """Read authentic installed cmcrameri batlow data without import-time registration.

    Importing cmcrameri itself registers all its colormaps globally. Reading its
    declared package data avoids that unrelated effect and retains real colors.
    """
    positions = np.asarray(positions, dtype=float)
    if not np.isfinite(positions).all() or np.any((positions < 0) | (positions > 1)):
        raise ValueError("batlow positions must be finite and in [0,1]")
    distribution = metadata.distribution("cmcrameri")
    data = np.loadtxt(distribution.locate_file("cmcrameri/cmaps/batlow.txt"))
    return colors.ListedColormap(data, name="batlow")(positions)


def _positive(value, name):
    """Reject bool, non-finite and nonpositive values for point and size settings."""
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")


def _overrides(overrides, defaults):
    """Translate validated appearance keys into rcParams and effective generic defaults."""
    if overrides is None:
        return {}, defaults
    if not isinstance(overrides, dict):
        raise ValueError("style overrides must be a dictionary")
    params = {}
    for key, value in overrides.items():
        if key in _ALIASES:
            if key == "figsize":
                if not isinstance(value, (list, tuple)) or len(value) != 2:
                    raise ValueError("figsize must contain width and height")
                for item in value:
                    _positive(item, "figsize")
            elif key == "bbox_inches":
                if value not in (None, "tight"):
                    raise ValueError("bbox_inches must be None or 'tight'")
            else:
                _positive(value, key)
            for rc_key in _ALIASES[key]:
                params[rc_key] = value
            defaults[key] = value
        elif key in _NUMERIC_RC:
            _positive(value, key)
            params[key] = value
        elif key in _COLOR_RC:
            if not colors.is_color_like(value):
                raise ValueError(f"invalid color for {key}")
            params[key] = value
        elif key in ("xtick.direction", "ytick.direction") and value in ("in", "out", "inout"):
            params[key] = value
        elif key == "legend.frameon" and isinstance(value, bool):
            params[key] = value
        else:
            raise ValueError(f"unsupported appearance override {key!r}")
    # Reflect direct rc overrides in the small generic-defaults view as well.
    for friendly, rc_keys in _ALIASES.items():
        if friendly not in overrides and all(key in params for key in rc_keys):
            values = [params[key] for key in rc_keys]
            if all(value == values[0] for value in values):
                defaults[friendly] = values[0]
    return params, defaults


@contextmanager
def style_context(name="personal", overrides=None, *, color_mode="color"):
    """Apply a preset and validated appearance overrides, restoring caller state.

    Yields name/version/font/text-engine/notes and effective generic defaults.
    Presets do not select scientific variables, scales, limits or transforms.
    """
    requested_name = name
    name = resolve_style(name)
    if not isinstance(color_mode, str) or color_mode not in COLOR_MODES:
        raise ValueError(f"color_mode must be one of {COLOR_MODES}")
    preset = _preset(name)
    changes, defaults = _overrides(overrides, preset["defaults"])
    # Reset inherited style-dependent settings inside a context before applying the preset.
    base = {key: mpl.rcParamsDefault[key] for key in _BASELINE_KEYS}
    with mpl.rc_context(base):
        requested = preset["font"]
        try:
            font_manager.findfont(requested, fallback_to_default=False)
            selected = requested
        except ValueError:
            selected = "DejaVu Sans"
        config = dict(preset["rc"])
        config.update({"font.family": "sans-serif", "font.sans-serif": [selected],
                       "text.usetex": False, "lines.color": "black"})
        if color_mode == "color":
            config["axes.prop_cycle"] = cycler(color=batlow_colors([.2, .7]))
        else:
            config["axes.prop_cycle"] = (cycler(color=["black"] * 4)
                                         + cycler(linestyle=["-", "--", "-.", ":"]))
        with mpl.rc_context(config):
            rendering = {"font": {"requested": requested, "selected": selected,
                                   "fallback": requested != selected}, "text_engine": "Mathtext"}
            for friendly, rc_keys in _ALIASES.items():
                for key in rc_keys:
                    mpl.rcParams[key] = defaults[friendly]
            # Explicit rc overrides take precedence over the baseline and friendly defaults.
            mpl.rcParams.update(changes)
            yield {"name": requested_name, "canonical_name": name, "color_mode": color_mode,
                   "version": __version__, "defaults": defaults,
                   **rendering, "notes": preset["notes"]}
