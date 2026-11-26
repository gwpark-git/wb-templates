"""Compare separate equilibrium simulation results, using explicit data paths."""

import argparse
from contextlib import nullcontext
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from bdmc.equilibrium import (
    HARD_SHELL_WIDTHS, carnahan_starling, contact_shell_z, hard_sphere_survival, virtual_compression_z,
)

COMPRESSIONS = (0.025, 0.05, 0.1, 0.25, 0.5, 1.0)


def sem(values):
    return float(np.std(values, ddof=1) / math.sqrt(len(values)))


def block_means(series, length):
    return np.concatenate([
        trace[:len(trace) // length * length].reshape(-1, length).mean(axis=1)
        for trace in series
    ])


def finite_read(path, *, shape, delimiter=None, skiprows=0):
    raw = np.loadtxt(path, delimiter=delimiter, ndmin=2, skiprows=skiprows)
    if raw.shape != shape or not np.isfinite(raw).all():
        raise ValueError(f"unexpected shape or nonfinite observations in {path}")
    return raw


def group_cases(plan, fields):
    groups = {}
    for case in json.loads(plan.read_text())["cases"]:
        if case["family"] == "main":
            key = tuple(case[field] for field in fields)
            groups.setdefault(key, []).append(case)
    if not groups:
        raise ValueError(f"no main cases in {plan}")
    for cases in groups.values():
        if len(cases) < 2 or len({(c["seed"], c.get("stream")) for c in cases}) != len(cases):
            raise ValueError("main cases require independent replica seeds")
    return sorted(groups.items())


def soft_rows(plan, runs):
    for (phi, a, side, dt), cases in group_cases(plan, ("phi", "A", "nside", "dt")):
        density = 6 * phi / math.pi
        traces = []
        for case in cases:
            folder = runs / case["id"]
            expected = {key: case[key] for key in ("seed", "phi", "A", "nside", "dt", "equilibration", "production", "sample_interval")}
            state = json.loads((folder / "state.json").read_text())
            if not state["complete"] or any(state[key] != value for key, value in expected.items()):
                raise ValueError(f"unverified or mismatched LAMMPS run {folder}")
            count = round(case["production"] / case["sample_interval"])
            stride = round(case["sample_interval"] / dt)
            first = round(case["equilibration"] / dt)
            raw = np.loadtxt(folder / "observables.dat", ndmin=2)
            # fix ave/time may also emit the equilibration endpoint.
            if raw.shape == (count + 1, 10) and raw[0, 0] == first:
                raw = raw[1:]
            if (raw.shape != (count, 10) or not np.isfinite(raw).all() or
                    not np.array_equal(raw[:, 0], first + stride * np.arange(1, count + 1))):
                raise ValueError(f"sample clock mismatch in {folder}")
            np.testing.assert_allclose(raw[:, 1], raw[:, 2:5].mean(axis=1), atol=2e-11, rtol=0)
            # With kBT=1, Z=1+pair pressure/rho. The ideal term is exact.
            traces.append(1 + raw[:, 1] / density)
        sampling = cases[0]["sample_interval"]
        z = float(np.mean([trace.mean() for trace in traces]))
        se = max(sem([trace.mean() for trace in traces]),
                 *(sem(block_means(traces, round(width / sampling))) for width in (0.5, 1, 2)))
        yield {"model": f"Soft BD A={a:g}", "edition": "new", "phi": phi, "N": side**3,
               "replicas": len(cases), "dt": dt, "Z": z, "Z_SE": se}


def hard_rows(plan, runs):
    for (phi, side, move), cases in group_cases(plan, ("phi", "nside", "step")):
        n = side**3
        tables = []
        for case in cases:
            stride = case["sample_stride"]
            count = case["production"] // stride
            raw = finite_read(runs / f"{case['id']}.csv", shape=(count, 9), delimiter=",", skiprows=1)
            if not np.array_equal(raw[:, 0], case["equilibration"] + stride * np.arange(1, count + 1)):
                raise ValueError("hard-sphere sample clock mismatch")
            if np.any(raw[:, 6] < 0) or np.any(np.diff(raw[:, 2:7], axis=1) < 0) or np.any((raw[:, 8] < 0) | (raw[:, 8] > 1)):
                raise ValueError("invalid hard-sphere shell counts or move acceptance")
            tables.append(raw)
        compression = {}
        uncertainties = {}
        for x in COMPRESSIONS:
            indicators = [hard_sphere_survival(table[:, 1], particles=n, compression=x) for table in tables]
            probabilities = np.asarray([trace.mean() for trace in indicators])
            estimates = virtual_compression_z(probabilities, x)
            compression[str(x)] = float(estimates.mean())
            spacing = cases[0]["sample_stride"]
            uncertainties[str(x)] = max(sem(estimates), *(
                sem(block_means(indicators, width // spacing)) / (x * probabilities.mean())
                for width in (500, 1000, 2000) if width >= spacing))
        contacts = {str(width): float(contact_shell_z(
            np.mean([table[:, k + 2].mean() for table in tables]), particles=n, width=width))
            for k, width in enumerate(HARD_SHELL_WIDTHS)}
        yield {"model": "Hard-sphere MC", "edition": "new", "phi": phi, "N": n,
               "replicas": len(cases), "move": move, "Z": compression["0.25"], "Z_SE": uncertainties["0.25"],
               "compression_Z": compression, "compression_SE": uncertainties, "contact_Z": contacts,
               "acceptance": float(np.mean([table[:, 8].mean() for table in tables]))}


def archive_rows(soft, hard):
    for path, is_soft in ((soft, True), (hard, False)):
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                if row["family"] == "main":
                    yield {"model": f"Soft BD A={float(row['A']):g}" if is_soft else "Hard-sphere MC",
                           "edition": "archive", "phi": float(row["phi"]), "N": int(row["N"]),
                           "replicas": int(row["replicas"]), "Z": float(row["Z"]),
                           "Z_SE": float(row["Z_SE"] if is_soft else row["SE"])}


def figure(rows, output, use_wb_style):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    style = nullcontext()
    if use_wb_style:
        from wb_style import style_context
        style = style_context()
    markers = {"Soft BD A=10": "o", "Soft BD A=100": "s", "Hard-sphere MC": "^"}
    colors = {"archive": "0.2", "new": "#c62828"}
    with style:
        fig, (top, bottom) = plt.subplots(2, 1, figsize=(7, 7), height_ratios=(1.5, 1), sharex=True)
        line = np.linspace(0, 0.325, 250)
        top.plot(line, carnahan_starling(line), color="black", linewidth=1.3)
        bottom.axhline(0, color="black", linewidth=1.0)
        for model, marker in markers.items():
            for edition, color in colors.items():
                series = sorted((row for row in rows if row["model"] == model and row["edition"] == edition),
                                key=lambda row: row["phi"])
                offset = -0.0016 if edition == "archive" else 0.0016
                x = np.array([row["phi"] for row in series])
                z = np.array([row["Z"] for row in series])
                se = np.array([row["Z_SE"] for row in series])
                opts = dict(marker=marker, linestyle="none", color=color, mfc="white", capsize=3, markersize=6)
                top.errorbar(x + offset, z, yerr=se, **opts)
                cs = carnahan_starling(x)
                bottom.errorbar(x + offset, 100 * (z / cs - 1), yerr=100 * se / cs, **opts)
        first = top.legend(handles=[Line2D([], [], color="black", label="Carnahan–Starling"),
                           *(Line2D([], [], marker=marker, ls="none", color="black",
                                    mfc="white", label=model) for model, marker in markers.items())],
                           loc="upper left", frameon=False)
        top.add_artist(first)
        top.legend(handles=[Line2D([], [], marker="o", mfc="white", ls="none", color=color, label=label)
                            for color, label in [(colors["archive"], "Archived"), (colors["new"], "New")]],
                   loc="lower right", frameon=False)
        top.set_ylabel(r"$Z=\Pi/(\rho k_\mathrm{B}T)$")
        bottom.set(xlabel=r"Nominal volume fraction $\phi$", ylabel="Deviation from CS (%)", xlim=(0, 0.325))
        for axes in (top, bottom):
            axes.tick_params(which="both", direction="in", top=True, right=True)
        fig.text(0.5, 0.01, "N = 216. Bars show one descriptive SE. Horizontal offsets separate editions.",
                 ha="center", fontsize=9)
        fig.tight_layout(rect=(0, 0.035, 1, 1))
        fig.savefig(output / "equation_of_state.pdf")
        fig.savefig(output / "equation_of_state.png", dpi=300)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("soft-plan", "soft-runs", "hard-plan", "hard-runs", "soft-archive", "hard-archive", "output"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--wb-style", action="store_true")
    args = parser.parse_args()
    rows = [*soft_rows(args.soft_plan, args.soft_runs), *hard_rows(args.hard_plan, args.hard_runs),
            *archive_rows(args.soft_archive, args.hard_archive)]
    keys = {(row["model"], row["edition"], row["phi"], row["N"]) for row in rows}
    if len(rows) != len(keys) or any((model, "archive" if edition == "new" else "new", phi, n) not in keys
                                    for model, edition, phi, n in keys):
        raise ValueError("new and archived comparison points do not match one to one")
    for row in rows:
        rho = 6 * row["phi"] / math.pi
        row.update(CS_Z=float(carnahan_starling(row["phi"])), pressure_reduced=rho * row["Z"],
                   pressure_SE=rho * row["Z_SE"])
        row["deviation_CS_percent"] = 100 * (row["Z"] / row["CS_Z"] - 1)
    args.output.mkdir()
    fields = ("model", "edition", "phi", "N", "replicas", "Z", "Z_SE", "CS_Z",
              "deviation_CS_percent", "pressure_reduced", "pressure_SE")
    with (args.output / "equation_of_state.csv").open("w", newline="") as stream:
        csv_out = csv.DictWriter(stream, fields, extrasaction="ignore")
        csv_out.writeheader()
        csv_out.writerows(rows)
    manifest = {}
    for name in ("soft_plan", "hard_plan", "soft_archive", "hard_archive"):
        path = getattr(args, name)
        manifest[name] = {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    (args.output / "equation_of_state.json").write_text(
        json.dumps({"editions": manifest, "points": rows}, allow_nan=False, indent=2) + "\n")
    figure(rows, args.output, args.wb_style)


if __name__ == "__main__":
    main()
