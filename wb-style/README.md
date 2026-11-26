# wb-style

Reusable Matplotlib appearance derived from **Gunwoo Park's researcher-authored
Python plotting scripts** used for the paper cited below. See [NOTICE](NOTICE.md)
for the source and adaptation history.

This package owns inward ticks, spines, typography, line/marker defaults and
explicit dash patterns. It does not load a research archive, interpret scientific
variables or reproduce paper data. The independent `plot-styles` package provides
the common personal color/monochrome API; neither package requires the other.

## Install and run a complete neutral example

Python 3.10+, NumPy1.24+ and Matplotlib3.6+ are required. From this child:

```sh
python3 -m pip install .
python3 -B examples/plot_example.py --output /tmp/new-wb-style-demo
python3 -B -m unittest discover -s tests -v
```

The example creates deterministic **synthetic** `series.csv` with columns
`x,series_a,series_b`, then PNG/PDF and `render.json` in a new output directory.
Expected: two labelled scalar series, inward ticks, hollow comparison markers,
`data_origin: synthetic`, and no scientific simulation. Existing destinations
are refused. Use another new directory for another run. To plot a deliberately
selected neutral CSV with that same header:

```sh
python3 -B examples/plot_example.py --input /selected/series.csv --output /tmp/new-wb-style-csv
```

The selected file is read only. Rows/values are plotted unchanged; no units,
normalization, physical parameters or archive filename conventions are inferred.

```python
import matplotlib.pyplot as plt
from wb_style import configuration, dash_pattern, style_context

with style_context():
    defaults = configuration()["plot_defaults"]
    fig, ax = plt.subplots(figsize=defaults["figsize_inches"])
    ax.plot([0, 1], [0, 1], color="black", lw=defaults["line_width_pt"],
            ls=dash_pattern("densely dashed"))
    ax.set(xlabel="x (synthetic)", ylabel="signal (synthetic)")
    fig.savefig("figure.pdf", bbox_inches="tight")
```

Caller rcParams are restored after normal/error exit. Arial is a source preference;
DejaVu Sans is the portable fallback, with Mathtext by default. `--tex` opts into
an existing TeX/mathrsfs installation. No font is bundled and paper-pixel identity
is not claimed.

## Reference and scope

G. W. Park, J. K. G. Dhont and G. Nägele, *The Journal of Chemical Physics*
**164**, 244119 (2026), [doi:10.1063/5.0336005](https://doi.org/10.1063/5.0336005).
Citation identifies the style's research source; the example data are new synthetic
values. See [CITATION.cff](CITATION.cff), [references and color licenses](REFERENCES.md),
[style/testing contract](docs/STYLE_AND_VALIDATION.org),
[provenance](provenance.json)
and the collection [MIT grant](LICENSE).

The distribution/import names are `wb-style` / `wb_style`. Replace earlier import
names deliberately; this source revision provides no alias or second runtime.
