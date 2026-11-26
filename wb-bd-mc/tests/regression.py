#!/usr/bin/env python3
"""Bounded portable-build regressions, not a legacy Intel reference comparison."""

import argparse
from collections import Counter
import datetime
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


CODE = Path(__file__).resolve().parents[1]
COMMON = {
    "Np": "64", "Nt": "101", "dt/tauR": "0.001", "dt/tauB": "0.001",
    "box_dimension": "6", "LBx": "6", "LBy": "6", "LBz": "6",
    "N_THREADS_BD": "1", "N_THREADS_SS": "1", "Rt": "0.001",
    "N_skip_ener": "10", "N_skip_file": "10", "record_RDIST": "FALSE",
    "friction_junction": "DEFAULT", "probability_area": "FALSE",
    "filename_base": "case", "MC_LOG": "FALSE",
}
CASES = {
    "brownian": {"Method": "BROWNIAN", "Np": "128", "Nt": "102",
                 "N_skip_file": "1", "N_skip_ener": "1"},
    "dumbbell": {"Method": "DUMBBELL", "cell_list": "FALSE"},
    "repulsive": {"Method": "REPULSIVE_BROWNIAN"},
    "association": {"Method": "NAPLE_ASSOCIATION", "MC_LOG": "TRUE"},
    "equilibration": {"Method": "NAPLE_ASSOCIATION", "Step": "EQUILIBRATION"},
    "association_bd2": {"Method": "NAPLE_ASSOCIATION", "N_THREADS_BD": "2"},
    "association_no_cells": {"Method": "NAPLE_ASSOCIATION", "cell_list": "FALSE"},
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def table(path):
    rows = [[float(x) for x in line.split()] for line in path.read_text().splitlines()]
    require(rows and all(row for row in rows), f"empty output: {path}")
    require(all(math.isfinite(x) for row in rows for x in row), f"non-finite output: {path}")
    return rows


def integers(row):
    require(all(x == int(x) for x in row), "non-integral connectivity")
    return [int(x) for x in row]


def topology(run, parameters, frames, energy):
    n = int(parameters["Np"])
    nc = int(parameters["N_chains_per_particle"])
    chains = table(run / "case.chain")
    hashes = table(run / "case.hash")
    weights = table(run / "case.weight")
    require(len(chains) == frames, "chain frame count")
    require(len(hashes) == len(weights) == n * frames, "connectivity frame count")
    max_bridges = 0
    for frame, row in enumerate(chains):
        ends = integers(row)
        require(len(ends) == 2 * n * nc, "chain ends were lost")
        require(all(0 <= p < n for p in ends), "invalid chain endpoint")
        # Reconstruct adjacency independently from the labelled chain endpoints.
        expected = [Counter() for _ in range(n)]
        bridges = 0
        for a, b in zip(ends[:n * nc], ends[n * nc:]):
            expected[a][b] += 1
            expected[b][a] += 1
            bridges += a != b
        max_bridges = max(max_bridges, bridges)
        require(energy[frame][4] == bridges, "energy/chain bridge count mismatch")
        for particle in range(n):
            neighbors = integers(hashes[frame * n + particle])
            counts = integers(weights[frame * n + particle])
            require(len(neighbors) == len(counts), "hash/weight length mismatch")
            require(len(set(neighbors)) == len(neighbors), "duplicate adjacency entry")
            require(neighbors[0] == particle, "missing self entry")
            require(all(0 <= p < n for p in neighbors), "invalid adjacency index")
            require(all(w >= 0 for w in counts), "negative connection weight")
            actual = Counter({p: w for p, w in zip(neighbors, counts) if w})
            require(actual == expected[particle], "adjacency/labelled-chain mismatch")
    require(max_bridges > 0, "association test never formed a bridge")
    return {"frames": frames, "max_bridges": max_bridges,
            "chain_end_count": 2 * n * nc}


def validate(run, parameters):
    trajectory = table(run / "case.traj")
    energy = table(run / "case.ener")
    n = int(parameters["Np"])
    dim = int(parameters["N_dimension"])
    skip = int(parameters["N_skip_file"])
    steps = range(0, int(parameters["Nt"]) - 1, skip)
    require(len(trajectory) == len(energy) == len(steps), "record count mismatch")
    method = parameters["Method"]
    association = method == "NAPLE_ASSOCIATION" and parameters["Step"] != "EQUILIBRATION"
    columns = 31 if association else 18 if method == "BROWNIAN" else 24
    require(all(len(row) == columns for row in energy), "energy schema changed")
    dt = float(parameters["dt/tauB" if method in ("BROWNIAN", "DUMBBELL") else "dt/tauR"])
    for row, step in zip(trajectory, steps):
        require(len(row) == 1 + 2 * n * dim, "trajectory schema changed")
        require(math.isclose(row[0], step * dt, abs_tol=1e-8), "trajectory timestamp")
        for p in range(n):
            for d, axis in enumerate(("LBx", "LBy", "LBz")):
                require(-1e-6 <= row[1 + 2 * dim * p + d] <= float(parameters[axis]) + 1e-6,
                        "position outside periodic box")
    require(trajectory[0][1:] != trajectory[-1][1:], "positions did not evolve")
    result = {"frames": len(trajectory), "energy_columns": columns}
    if method == "BROWNIAN":
        # The legacy uniform noise is variance-normalised: <|delta r|^2> = 2*d*dt.
        squared = []
        for previous, current in zip(trajectory, trajectory[1:]):
            for p in range(n):
                delta2 = 0.
                for d, axis in enumerate(("LBx", "LBy", "LBz")):
                    length = float(parameters[axis])
                    dx = current[1 + 2 * dim * p + d] - previous[1 + 2 * dim * p + d]
                    dx -= length * round(dx / length)
                    delta2 += dx * dx
                squared.append(delta2)
        ratio = sum(squared) / len(squared) / (2 * dim * dt)
        require(abs(ratio - 1) < 0.06, f"Brownian displacement variance: {ratio}")
        result["brownian_variance_ratio"] = ratio
    if association:
        result.update(topology(run, parameters, len(trajectory), energy))
    return result


def execute(executable, name, overrides, output):
    run = output / name
    run.mkdir()
    parameters = {}
    lines = []
    replacements = {**COMMON, **overrides, "output_path": str(run)}
    for line in (CODE / "short_test.inp").read_text().splitlines():
        if "=" in line and not line.startswith("#"):
            key, value = line.split("=", 1)
            value = replacements.get(key, value)
            parameters[key] = value
            line = key + "=" + value
        lines.append(line)
    for key, value in replacements.items():
        if key not in parameters:
            lines.append(key + "=" + value)
            parameters[key] = value
    inp = run / "input.inp"
    inp.write_text("\n".join(lines) + "\n")
    (run / "parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    env = {**os.environ, "OMP_DYNAMIC": "FALSE", "OMP_THREAD_LIMIT": "8",
           "GSL_RNG_TYPE": "mt19937"}
    with (run / "stdout.log").open("w") as log:
        process = subprocess.run([str(executable), str(inp)], cwd=run, env=env,
                                 stdout=log, stderr=subprocess.STDOUT, timeout=60)
    require(process.returncode == 0, f"{name}: exit {process.returncode}; see {run}/stdout.log")
    text = (run / "stdout.log").read_text()
    require("Bad condition" not in text and "ERR:" not in text, f"{name}: input/runtime diagnostic")
    return run, validate(run, parameters)


def compare(left, right):
    for suffix in ("traj", "ener", "chain", "hash", "weight"):
        a, b = left / ("case." + suffix), right / ("case." + suffix)
        require(a.exists() == b.exists(), f"{suffix}: missing counterpart")
        if not a.exists() or not a.stat().st_size:
            continue
        aa, bb = table(a), table(b)
        require(len(aa) == len(bb), f"{suffix}: record count")
        for row_a, row_b in zip(aa, bb):
            require(len(row_a) == len(row_b), f"{suffix}: row length")
            for column, (x, y) in enumerate(zip(row_a, row_b)):
                if suffix == "ener" and column == 5:
                    continue  # Only this field contains wall-clock time.
                require(math.isclose(x, y, rel_tol=2e-6, abs_tol=1e-9),
                        f"{suffix}: numeric difference at column {column}: {x} vs {y}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--run-dir", type=Path, default=CODE / "runs/regression")
    parser.add_argument("--compare-executable", type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve()
    args.run_dir.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix="run-", dir=args.run_dir.resolve()))
    report = {"created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "executable": str(executable), "cases": {}, "comparisons": [],
              "scope": "Portable-build checks; no saved legacy Intel numerical reference available."}
    try:
        help_result = subprocess.run([str(executable)], capture_output=True, text=True, timeout=10)
        require(help_result.returncode == 0 and "USAGE" in help_result.stdout, "help interface")
        directories = {}
        for name, overrides in CASES.items():
            directories[name], report["cases"][name] = execute(executable, name, overrides, output)
            print("PASS:", name, report["cases"][name], flush=True)
        # Fixed-thread repeatability; compare all scientific fields, excluding timings.
        for name in ("brownian", "association"):
            repeat, _ = execute(executable, name + "_repeat", CASES[name], output)
            compare(directories[name], repeat)
            report["comparisons"].append(name + ": repeated seeded run")
        if args.compare_executable:
            reference = args.compare_executable.resolve()
            report["compare_executable"] = str(reference)
            for name, overrides in CASES.items():
                if name == "association_bd2":
                    continue  # Parallel reduction order is not part of this comparison.
                other, _ = execute(reference, name + "_compare", overrides, output)
                compare(directories[name], other)
                report["comparisons"].append(name + ": other compiler")
        report["status"] = "passed"
        print("PASS: repeatability and requested compiler comparisons", flush=True)
    except Exception as exc:
        report["status"] = "failed"
        report["error"] = str(exc)
        raise
    finally:
        (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        print("Report:", output / "report.json", flush=True)


if __name__ == "__main__":
    main()
