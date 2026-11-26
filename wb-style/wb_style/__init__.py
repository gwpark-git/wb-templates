"""Settings derived from Gunwoo Park's researcher-authored plotting scripts.

The original paper plotting style and scripts are credited in NOTICE.md.
"""
from __future__ import annotations

from contextlib import contextmanager
from importlib.resources import files
import json
from shutil import which

import matplotlib as mpl
from matplotlib import font_manager

__version__ = "0.2.0"


def configuration() -> dict:
    """Return a fresh editable copy of the measured source style settings."""
    return json.loads(files(__package__).joinpath("style.json").read_text(encoding="utf-8"))


def dash_pattern(name: str):
    """Return zero offset and point lengths from style.json; an unknown name raises KeyError."""
    return (0, tuple(configuration()["dash_patterns"][name]))


def font_selection() -> dict:
    """Select the requested installed font or the configured fallback and report both."""
    config = configuration()
    original = config["original_font"]
    try:
        font_manager.findfont(original, fallback_to_default=False)
        selected = original
    except ValueError:
        selected = config["portable_fallback_font"]
    return {"requested": original, "selected": selected, "fallback": selected != original}


@contextmanager
def style_context(*, tex: bool = False):
    """Apply the source style temporarily; restore caller rcParams on exit.

    Portable default uses Mathtext; opt-in TeX needs an installed working TeX
    distribution including mathrsfs. Neither mode promises paper pixel identity.
    """
    if tex and which("latex") is None:
        raise RuntimeError("TeX mode requires latex and mathrsfs; use the portable default")
    config = configuration()
    font = font_selection()
    params = dict(config["rc_params"])
    params.update({"font.sans-serif": [font["selected"]], "text.usetex": tex})
    if tex:
        params["text.latex.preamble"] = r"\usepackage{mathrsfs}"
    # Use Matplotlib context restoration for normal, exceptional and nested use.
    with mpl.rc_context(params):
        yield {"font": font, "text_engine": "TeX" if tex else "Mathtext"}
