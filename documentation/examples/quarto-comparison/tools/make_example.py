"""Generate the synthetic two-implementation data and comparison assets."""

from pathlib import Path
import argparse
import csv
import difflib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wb_plot_styles import style_context

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "models"))
import candidate
import reference

VOLUME_FRACTIONS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.4)
QUANTITIES = ("D/D0", "eta/eta_s", "Z")
DISPLAY_QUANTITIES = (r"$D/D_0$", r"$\eta/\eta_s$", r"$Z$")
TOLERANCE = 1e-14


def compare():
    """Pair the two implementations at each concentration and check their residuals."""
    cases = [(phi, reference.properties(phi), candidate.properties(phi))
             for phi in VOLUME_FRACTIONS]
    residuals = [max(abs(current[i] - expected[i]) for _, expected, current in cases)
                 for i in range(len(QUANTITIES))]
    if any(delta > TOLERANCE for delta in residuals):
        raise ValueError("Synthetic arithmetic regression exceeds the example tolerance")
    return cases, residuals


def write_tables(destination, cases, residuals):
    """Store all paired values and signed residuals and summarize maxima by quantity."""
    with (destination / "data/comparison.csv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["phi", *(f"{name}_{item}" for name in QUANTITIES
                                  for item in ("reference", "candidate", "delta"))])
        for phi, expected, current in cases:
            values = [value for a, b in zip(expected, current) for value in (a, b, b - a)]
            writer.writerow([f"{value:.16e}" for value in (phi, *values)])
    rows = [
        "## Physical-quantity comparison {#physical-quantity-comparison}", "",
        "The personal monochrome preset sets the figure style. Black open circles show the reference and red crosses show the candidate.", "",
        "![Three synthetic reduced quantities versus volume fraction.](figures/comparison.png)", "",
        "[Numerical CSV](data/comparison.csv) / [Figure PDF](figures/comparison.pdf)", "",
        "| Reduced quantity | Maximum absolute difference | Tolerance |", "|---|---:|---:|",
        *(f"| {name} | {error:.3e} | {TOLERANCE:.1e} |" for name, error in zip(DISPLAY_QUANTITIES, residuals)),
        "",
    ]
    (destination / "_results.md").write_text("\n".join(rows), encoding="utf-8")


def write_diff(destination):
    """Generate the unified diff from the named files that produced the data."""
    left, right = (ROOT / "models" / f"{name}.py" for name in ("reference", "candidate"))
    diff = difflib.unified_diff(left.read_text().splitlines(), right.read_text().splitlines(),
                                fromfile="models/reference.py", tofile="models/candidate.py",
                                lineterm="")
    rows = [
        "## Source-code comparison {#source-code-comparison}", "",
        "The changed operations are the viscosity and osmotic-factor arithmetic.",
        "Compare the input, outputs and residual tolerance before changing a solver.", "",
        "[Reference source](models/reference.py) / [Candidate source](models/candidate.py)", "",
        "```diff", *diff, "```", "",
    ]
    (destination / "_code.md").write_text("\n".join(rows), encoding="utf-8")


def write_figures(destination, cases):
    """Render both formats with the same preset and Mathtext notation."""
    x = [case[0] for case in cases]
    with style_context("personal", {"figsize": (13, 4.4)}, color_mode="monochrome") as style:
        fig, axs = plt.subplots(1, len(QUANTITIES), layout="constrained")
        for index, (ax, quantity) in enumerate(zip(axs, DISPLAY_QUANTITIES)):
            ax.plot(x, [expected[index] for _, expected, _ in cases],
                    "o-", color="black", markerfacecolor="white", label="Reference")
            ax.plot(x, [current[index] for _, _, current in cases],
                    "x", color="#b72020", label="Candidate")
            ax.set(xlabel=r"Volume fraction $\phi$", ylabel=quantity)
        axs[0].legend()
        for suffix in ("png", "pdf"):
            fig.savefig(destination / f"figures/comparison.{suffix}")
        plt.close(fig)
        print(f"Figure style: {style['canonical_name']} / {style['color_mode']} / {style['text_engine']}")


def main():
    """Write a complete set of example assets to a new destination."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    cases, residuals = compare()
    args.destination.mkdir()
    for name in ("data", "figures"):
        (args.destination / name).mkdir()
    write_tables(args.destination, cases, residuals)
    write_diff(args.destination)
    write_figures(args.destination, cases)


if __name__ == "__main__":
    main()
