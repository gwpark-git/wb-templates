"""Neutral synthetic CSV example, unrelated to publication datasets.

The original researcher-authored paper style is attributed in NOTICE.md.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from wb_style import configuration, dash_pattern, style_context

HEADER = ("x", "series_a", "series_b")


def synthetic_series():
    """Create an 81-row array of x and two deterministic, unitless display series."""
    x = np.linspace(0, 4, 81)
    return np.column_stack((x, np.exp(-x), .7 * np.exp(-x / 2) + .1 * np.cos(2 * np.pi * x)))


def load_series(path):
    """Read the exact three-column UTF-8 CSV contract and reject missing or non-finite cells."""
    with Path(path).open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        if tuple(next(reader, ())) != HEADER:
            raise ValueError("Expected exactly x,series_a,series_b CSV columns")
        try:
            data = np.asarray([[float(v) for v in row] for row in reader], dtype=float)
        except (ValueError, TypeError):
            raise ValueError("CSV rows must contain three finite numbers") from None
    if data.ndim != 2 or data.shape[1] != 3 or len(data) < 2 or not np.isfinite(data).all():
        raise ValueError("CSV needs at least two rows of three finite numbers")
    return data


def marker_indices(n_rows, stride=8):
    """Select every stride-th row plus the final row without resampling coordinates."""
    if n_rows < 1 or stride < 1:
        raise ValueError("Row count and stride must be positive")
    return np.unique(np.append(np.arange(0, n_rows, stride), n_rows - 1))


def render(output, input_path=None, *, tex=False):
    """Validate input before creating output; write CSV, PNG, PDF and a hashed run record."""
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError("Select a new output directory")
    if output.is_relative_to(Path(__file__).resolve().parents[1]):
        raise ValueError("Keep generated outputs outside the maintained package")
    data = load_series(input_path) if input_path is not None else synthetic_series()
    defaults = configuration()["plot_defaults"]
    with style_context(tex=tex) as rendering:
        fig, ax = plt.subplots(figsize=defaults["figsize_inches"])
        try:
            # Marker x coordinates come from the selected rows, including nonuniform grids.
            selected = marker_indices(len(data))
            for col, label, dash, marker in [(1, "Series A", "densely dashed", "o"), (2, "Series B", "dashdotted", "s")]:
                ax.plot(data[:, 0], data[:, col], color=defaults["line_color"], lw=defaults["line_width_pt"], ls=dash_pattern(dash), label=label)
                ax.plot(data[selected, 0], data[selected, col], ls="none", marker=marker, color=defaults["line_color"], mfc=defaults["marker_face_color"], ms=defaults["marker_size_pt"], mew=defaults["marker_edge_width_pt"])
            ax.set_xlabel("x (synthetic)" if input_path is None else "x", fontsize=defaults["label_size_pt"])
            ax.set_ylabel("signal (synthetic)" if input_path is None else "signal", fontsize=defaults["label_size_pt"])
            ax.tick_params(labelsize=defaults["tick_label_size_pt"])
            ax.legend()
            output.mkdir(parents=True)
            # Use 17 significant digits to preserve finite binary64 values on CSV readback.
            np.savetxt(output / "series.csv", data, delimiter=",", header=",".join(HEADER), comments="", fmt="%.17g")
            fig.savefig(output / "series.png", dpi=150, bbox_inches="tight")
            fig.savefig(output / "series.pdf", bbox_inches="tight")
        finally:
            plt.close(fig)
    summary = {"data_origin": "synthetic" if input_path is None else "explicit_selected_csv", "rows": len(data), "columns": list(HEADER), "series_count": 2, "row_order_and_values_transformed": False, "rendering": rendering, "scientific_simulation_rerun": False, "pixel_equivalence": "not claimed", "csv_sha256": hashlib.sha256((output / "series.csv").read_bytes()).hexdigest()}
    (output / "render.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    """Parse the optional explicit data file and launch the selected rendering mode."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="Explicit neutral three-column CSV; default creates synthetic data")
    parser.add_argument("--output", type=Path, required=True, help="New output directory outside maintained source")
    parser.add_argument("--tex", action="store_true", help="Opt into existing TeX/mathrsfs")
    args = parser.parse_args()
    print(json.dumps(render(args.output, args.input, tex=args.tex), indent=2))


if __name__ == "__main__":
    main()
