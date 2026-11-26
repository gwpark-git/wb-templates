# wb plot-styles

One reusable personal Matplotlib style, with **color** and **monochrome** modes.
Both modes share geometry, typography, line widths and inward ticks. Source
author/year labels belong to provenance, not to style selection.

```python
import matplotlib.pyplot as plt
from wb_plot_styles import style_context

with style_context(color_mode="monochrome") as rendering:
    fig, ax = plt.subplots(figsize=rendering["defaults"]["figsize"])
    ax.plot([0, 1], [0, 1], marker="o", label="A")
    ax.plot([0, 1], [1, 0], marker="s", label="B")
    ax.set(xlabel="x", ylabel="y")
    ax.legend()
    fig.savefig("figure.png", dpi=rendering["defaults"]["dpi"], bbox_inches="tight")
```

`style_context()` defaults to `color_mode="color"` and the common `personal`
style. Color uses authentic batlow colors. Monochrome uses black and cycles
solid, dashed, dash-dot and dotted lines; callers select markers/fills explicitly.
A mode supplies defaults: explicit artist colors/styles still win.
`plot_defaults()` returns size/export defaults. For example,
`style_context(overrides={"markersize": 10}, color_mode="color")` changes marker
size. Caller settings are restored on exit, including after an exception.

The retained common baseline is 6×6 inches, 2.5 pt lines, 15 pt labels/ticks,
inward ticks and batlow .2/.7. The 300 dpi export setting was added explicitly;
it was not recovered from original scripts. This is a reusable default selected
from the researcher's helpers, not a claim that every historical figure matched.

Version 0.2.0 removes author/year selectors and aliases. Earlier scripts that
explicitly selected an edition need migration to the common style plus the
desired color mode and explicit size overrides. Python package/import names stay
`wb-plot-styles` / `wb_plot_styles`. No plotting companion or agent is required.

Install this child with Python 3.10+, NumPy 1.24+, Matplotlib 3.6+ and cmcrameri
1.10+ in the chosen environment. Then, from an independent copy:

```sh
python3 -m pip install .
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
PYTHONDONTWRITEBYTECODE=1 python3 examples/render_examples.py --output /tmp/wb-style-demos
```

The output directory must be new. Synthetic examples show both color modes,
vertically merged axes, multiple symbols per legend entry and a matrix legend.
They do not reproduce scientific calculations. For CSV parsing use the separate
`wb-scientific-io` adapter. Data meaning, units, transformations, scales and axis
ranges remain explicit standard Matplotlib/user choices.

[Developer/API](docs/DEVELOPER.org) / [Provenance](provenance.json) /
[Attribution](NOTICE.md) / [MIT license](LICENSE)

Project-owned style and adaptations use the collection MIT grant. The Matplotlib
and cmcrameri notices carry their upstream terms.

Gunwoo Park authored the original scripts and style and directed the reusable adaptation.

[Style/paper and batlow references](REFERENCES.md). The actual checked cmcrameri1.10 MIT notice is retained in [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES/cmcrameri-1.10-LICENSE.txt).
