"""Synthetic source-derived layouts and ordinary preset comparisons; no research data."""
from pathlib import Path
import argparse
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from wb_plot_styles import COLOR_MODES, batlow_colors, plot_defaults, style_context
from wb_plot_styles.helpers import HandlerXoffset, draw_matrix_legend


def save(fig, directory, name, dpi=150):
    """Write a raster PNG and a vector PDF for a named figure, then release it."""
    fig.savefig(directory / f"{name}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(directory / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def vertical(directory):
    """Render aligned synthetic sine and cosine panels and return the applied style metadata."""
    with style_context() as context:
        palette = batlow_colors([.2, .7])
        x = np.linspace(0, 2 * np.pi, 100)
        fig, axes = plt.subplots(2, 1, sharex=True, figsize=(6, 6), gridspec_kw={"hspace": 0})
        for ax, values, label, color in zip(axes, (np.sin(x), np.cos(x)), ("$f_1(x)$", "$f_2(x)$"), palette):
            ax.plot(x / np.pi, values, color=color, lw=2.5)
            ax.set_ylabel(label)
            ax.set(xlim=(0, 2), ylim=(-1.1, 1.1))
        axes[-1].set_xlabel(r"$x/\pi$ (synthetic)")
        axes[-1].set_xticks([0, .5, 1, 1.5, 2])
        save(fig, directory, "vertical_merged")
    return context


def multi_symbol(directory):
    """Demonstrate two grouped marker legends using signed handler offsets."""
    with style_context(color_mode="monochrome"):
        fig, ax = plt.subplots(figsize=(8, 6))
        x = np.linspace(0, 1, 40)
        q, b = [], []
        for i, marker in enumerate(("o", "s", "^", "d")):
            q.append(ax.plot(x[::5], (x ** 2 + i / 12)[::5], ls="none", marker=marker, mfc="white", mec="black")[0])
            b.append(ax.plot(x[::5], (1 - x ** 2 - i / 12)[::5], ls="none", marker=marker, mfc="black", mec="black")[0])
        p1, = ax.plot(x, x ** 2, "k-")
        p2, = ax.plot(x, 1 - x ** 2, "k--")
        first = ax.legend(list(zip(q, b)), [f"$a=a_{i}$" for i in range(4)],
                           handler_map={**{item: HandlerXoffset(x_offset=6) for item in q},
                                        **{item: HandlerXoffset(x_offset=-6) for item in b}},
                           loc="upper right", fontsize=15, frameon=False,
                           handlelength=1, markerscale=1, handletextpad=1)
        ax.add_artist(first)
        offsets = (18, 6, -6, -18)
        ax.legend([tuple(q), p1, tuple(b), p2],
                   ["$y_1$: symbols", "$y_1$: curve", "$y_2$: symbols", "$y_2$: curve"],
                   handler_map={item: HandlerXoffset(x_offset=offset) for group in (q, b) for item, offset in zip(group, offsets)},
                   loc="upper left", fontsize=15, frameon=False, handlelength=3, handletextpad=.5)
        # Reserve synthetic whitespace above the largest symbols for both legends.
        ax.set(xlabel="$x$ (synthetic)", ylabel="$y$ (synthetic)", ylim=(-.4, 2))
        save(fig, directory, "multi_symbol_legend")


def matrix(directory):
    """Render synthetic decays and a factor-by-series inset legend; x uses a log axis."""
    with style_context():
        fig, ax = plt.subplots(figsize=(6, 4))
        x = np.logspace(-3, 2, 60)
        a_colors = batlow_colors(np.linspace(.1, .4, 3))
        b_colors = batlow_colors(1 - np.linspace(.1, .4, 3))
        for i, tau in enumerate((.1, 1, 10)):
            ax.plot(x, np.exp(-(x / tau) ** .7), ls="none", marker="o", ms=8, mec="k", mew=1, color=a_colors[i])
            ax.plot(x, np.exp(-(x / tau) ** 1.5), ls="none", marker="s", ms=7, mec="k", mew=1, color=b_colors[i])
        ax.set(xscale="log", xlabel="$x$ (synthetic)", ylabel="$f(x)$ (synthetic)")
        legend_ax = inset_axes(ax, width=1.8, height=1.2, loc="upper left",
                               bbox_to_anchor=(.08, .38), bbox_transform=ax.transAxes, borderpad=0)
        legend_ax.set_axis_off()
        draw_matrix_legend(legend_ax, ["0.1", "1.0", "10.0"], ["A", "B"], ["o", "s"],
                            [[a_colors[i], b_colors[i]] for i in range(3)],
                            fontsize=12, xscale=.55, yscale=3.3, gaps=(.1, .26, .25), row_header_text=r"$\tau$")
        save(fig, directory, "matrix_legend")


def comparisons(directory):
    """Render the same synthetic curves in each color mode and record each context."""
    metadata = []
    x = np.linspace(0, 1, 40)
    for mode in COLOR_MODES:
        defaults = plot_defaults()
        with style_context(color_mode=mode) as context:
            fig, ax = plt.subplots(figsize=defaults["figsize"])
            ax.plot(x, x, label="curve A", marker="o", markevery=5)
            ax.plot(x, x ** 2, label="curve B", ls="--", marker="s", markevery=5)
            ax.set(xlabel="$x$ (synthetic)", ylabel="$y$ (synthetic)")
            ax.legend(loc="upper left")
            save(fig, directory, f"personal_{mode}", dpi=defaults["dpi"])
            metadata.append(context)
    return metadata


def main():
    """Require a fresh output directory, write the figures and emit rendering.json."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new output directory")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    metadata = {"synthetic": True, "vertical": vertical(args.output)}
    multi_symbol(args.output)
    matrix(args.output)
    metadata["presets"] = comparisons(args.output)
    (args.output / "rendering.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
