"""Render an explicit synthetic CSV using the common personal style and color mode.

No scientific model, unit inference, original data or originating path is used.
"""
from pathlib import Path
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wb_scientific_io.csv_table import read_csv
from wb_plot_styles import style_context


def main():
    """Write a neutral CSV, load explicit x/y columns and render paired output using one style mode."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new output directory")
    parser.add_argument("--color-mode", choices=("color", "monochrome"), default="color")
    args = parser.parse_args()
    # Resolve style/dependencies before output creation; never overwrite an old run.
    with style_context(color_mode=args.color_mode) as rendering:
        args.output.mkdir(parents=True, exist_ok=False)
        data_path = args.output / "synthetic.csv"
        data_path.write_text("x,y,label\n0,0,synthetic\n0.5,0.25,synthetic\n1,1,synthetic\n", encoding="utf-8")
        data = read_csv(data_path, columns=["x", "y"])
        fig, ax = plt.subplots(figsize=rendering["defaults"]["figsize"])
        ax.plot(data["columns"]["x"], data["columns"]["y"], marker="o")
        ax.set(xlabel="x (synthetic)", ylabel="y (synthetic)")
        for suffix in ("png", "pdf"):
            fig.savefig(args.output / f"figure.{suffix}", dpi=rendering["defaults"]["dpi"], bbox_inches="tight")
        plt.close(fig)
        (args.output / "rendering.json").write_text(json.dumps({"synthetic": True, "style": rendering}, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
