#!/usr/bin/env python3
"""Frozen current-port/GSL regressions; not independent physical/Intel truth.

``create`` is an explicit, non-overwriting operation. ``verify`` authenticates
the saved files first, then runs the executable in a new directory. Neither
operation changes existing reference data. Only energy column 5 is timing.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading

CODE = Path(__file__).resolve().parents[1]
PROJECT = CODE

SCHEMA = "bdmc.current-port-reference.v1"
SCOPE = ("Frozen current-port/GSL behavior; not an independent physical oracle "
         "or the deferred historical SSD/Intel short_test.inp comparison.")
POLICY = {"relative_tolerance": 2e-6, "absolute_tolerance": 1e-9,
          "ignored_numeric_fields": {"ener": [5]},
          "integer_tables": ["chain", "hash", "weight"],
          "rationale": "Legacy tables print limited precision; only energy column 5 is wall-clock time."}
ENV = {"OMP_DYNAMIC": "FALSE", "OMP_THREAD_LIMIT": "8", "OMP_NUM_THREADS": "1",
       "OMP_SCHEDULE": "static", "GSL_RNG_TYPE": "mt19937", "GSL_RNG_SEED": "0",
       "LC_ALL": "C"}
COMMON = {"Np": "64", "Nt": "101", "dt/tauR": "0.001", "dt/tauB": "0.001",
          "N_dimension": "3", "box_dimension": "6", "LBx": "6", "LBy": "6", "LBz": "6",
          "N_THREADS_BD": "1", "N_THREADS_SS": "1", "Rt": "0.001",
          "N_skip_ener": "10", "N_skip_file": "10", "record_RDIST": "FALSE",
          "friction_junction": "DEFAULT", "probability_area": "FALSE",
          "filename_base": "case", "MC_LOG": "FALSE", "output_path": ".",
          "basic_random_seed": "100", "basic_random_seed_SS": "150",
          "N_chains_per_particle": "10", "tracking_individual_chain": "TRUE",
          "MECHANICAL_PERTURBATION": "FALSE", "shear_axis": "0", "shear_grad_axis": "1"}
CASES = {
    "brownian": {"Method": "BROWNIAN", "Np": "128", "Nt": "102",
                 "N_skip_ener": "1", "N_skip_file": "1"},
    "dumbbell": {"Method": "DUMBBELL", "cell_list": "FALSE"},
    "repulsive": {"Method": "REPULSIVE_BROWNIAN"},
    "equilibration": {"Method": "NAPLE_ASSOCIATION", "Step": "EQUILIBRATION"},
    "network": {"Method": "NAPLE_ASSOCIATION", "cell_list": "FALSE"},
    "network_simple_positive": {"Method": "NAPLE_ASSOCIATION", "cell_list": "FALSE",
                                "MECHANICAL_PERTURBATION": "SIMPLE_SHEAR", "Wi_tau_C": "0.5"},
    "network_simple_negative": {"Method": "NAPLE_ASSOCIATION", "cell_list": "FALSE",
                                "MECHANICAL_PERTURBATION": "SIMPLE_SHEAR", "Wi_tau_C": "-0.5"},
    "network_step_positive": {"Method": "NAPLE_ASSOCIATION", "cell_list": "FALSE",
                              "MECHANICAL_PERTURBATION": "STEP_SHEAR", "gamma_0": "0.5"},
    "network_step_negative": {"Method": "NAPLE_ASSOCIATION", "cell_list": "FALSE",
                              "MECHANICAL_PERTURBATION": "STEP_SHEAR", "gamma_0": "-0.5"},
}
SUFFIXES = ("traj", "ener", "chain", "hash", "weight")
CASE_BYTES = 4_000_000  # Conservative per-process bound, including stdout.
CORPUS_BYTES = 24_000_000
DIAGNOSTIC = re.compile(r"Bad condition|\bERR\b|\bERROR\b|\bFATAL\b|runtime error:|"
                        r"AddressSanitizer|UndefinedBehaviorSanitizer|Assertion.*failed|"
                        r"terminate called|Segmentation fault", re.IGNORECASE)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def parse_input(text):
    result = {}
    for line in text.splitlines():
        if line.strip() and not line.startswith("#") and "=" in line:
            key, value = line.split("=", 1)
            require(key not in result, f"duplicate input parameter: {key}")
            result[key] = value
    return result


def input_text(overrides):
    parameters = parse_input((CODE / "short_test.inp").read_text())
    parameters.update(COMMON)
    parameters.update(overrides)
    return "# Frozen current-port reference input; output path is relative to run cwd.\n" + "".join(
        f"{key}={value}\n" for key, value in sorted(parameters.items()))


def table(path, *, empty=False):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), f"missing table: {path}")
    if not path.stat().st_size:
        require(empty, f"empty table: {path}")
        return []
    rows = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        require(bool(line.strip()), f"blank table row: {path}:{number}")
        try:
            row = [float(value) for value in line.split()]
        except ValueError as exc:
            raise AssertionError(f"nonnumeric table: {path}:{number}") from exc
        require(all(math.isfinite(value) for value in row), f"nonfinite table: {path}:{number}")
        rows.append(row)
    return rows


def integer_row(row):
    require(all(value == int(value) for value in row), "nonintegral topology")
    return [int(value) for value in row]


def is_network(parameters):
    return parameters["Method"] == "NAPLE_ASSOCIATION" and parameters["Step"] != "EQUILIBRATION"


def validate_run(run, parameters, behavior="legacy-v1"):
    """Validate schemas and invariants independently of saved numerical values."""
    require(behavior in ("legacy-v1", "repaired-v2"), "unsupported behavioral revision")
    run = Path(run)
    log = (run / "stdout.log").read_text()
    match = DIAGNOSTIC.search(log)
    require(match is None, f"input/runtime diagnostic: {match.group(0) if match else ''}")
    require("START SIMULATION" in log, "missing simulation-start diagnostic")
    require("Total simulation time" in log, "missing simulation-completion diagnostic")
    network = is_network(parameters)
    n, dim = int(parameters["Np"]), int(parameters["N_dimension"])
    steps = list(range(0, int(parameters["Nt"]) - 1, int(parameters["N_skip_file"])))
    require(parameters["N_skip_file"] == parameters["N_skip_ener"], "suite requires matched recording intervals")
    require(n > 0 and dim == 3 and len(steps) > 1, "unsupported reference shape")
    arrays = {suffix: table(run / f"case.{suffix}", empty=suffix in ("chain", "hash", "weight"))
              for suffix in SUFFIXES}
    trajectory, energy = arrays["traj"], arrays["ener"]
    columns = 31 if network else 18 if parameters["Method"] == "BROWNIAN" else 24
    require(len(trajectory) == len(energy) == len(steps), "table frame count mismatch")
    require(all(len(row) == 1 + 2 * n * dim for row in trajectory), "trajectory schema mismatch")
    require(all(len(row) == columns for row in energy), "energy schema mismatch")
    dt_key = "dt/tauB" if parameters["Method"] in ("BROWNIAN", "DUMBBELL") else "dt/tauR"
    dt = float(parameters[dt_key])
    lengths = [max(float(parameters["box_dimension"]), float(parameters[axis]))
               for axis in ("LBx", "LBy", "LBz")]
    for frame, step in enumerate(steps):
        require(math.isclose(trajectory[frame][0], step * dt, rel_tol=0, abs_tol=1e-8),
                "trajectory timestamp mismatch")
        # Preserve v1's zero Brownian clock when authenticating that corpus;
        # repaired v2 requires the actual pre-step force/frame timestamp.
        energy_time = (0. if parameters["Method"] == "BROWNIAN" and behavior == "legacy-v1"
                       else step * dt)
        require(math.isclose(energy[frame][0], energy_time, rel_tol=0, abs_tol=1e-8),
                "energy time-field mismatch")
        require(energy[frame][5] >= 0, "negative wall-clock time")
        for particle in range(n):
            for axis, length in enumerate(lengths):
                value = trajectory[frame][1 + 2 * dim * particle + axis]
                require(-1e-6 <= value <= length + 1e-6, "position outside periodic box")
    require(trajectory[0][1:] != trajectory[-1][1:], "positions did not evolve")
    summary = {"frames": len(steps), "trajectory_columns": 1 + 2 * n * dim,
               "energy_columns": columns}
    if parameters["Method"] == "BROWNIAN":
        summary["energy_time_field"] = ("zero in current implementation" if behavior == "legacy-v1"
                                        else "pre-step force/frame timestamp")
    if network:
        nc = int(parameters["N_chains_per_particle"])
        chains, hashes, weights = arrays["chain"], arrays["hash"], arrays["weight"]
        require(len(chains) == len(steps), "chain frame count mismatch")
        require(len(hashes) == len(weights) == n * len(steps), "adjacency frame count mismatch")
        bridge_counts = []
        for frame, row in enumerate(chains):
            ends = integer_row(row)
            require(len(ends) == 2 * n * nc, "chain endpoint count mismatch")
            require(all(0 <= end < n for end in ends), "chain endpoint outside particle range")
            expected = [Counter() for _ in range(n)]
            bridges = 0
            for a, b in zip(ends[:n * nc], ends[n * nc:]):
                expected[a][b] += 1
                expected[b][a] += 1
                bridges += a != b
            bridge_counts.append(bridges)
            require(energy[frame][4] == bridges, "energy/labelled-chain bridge mismatch")
            different_pairs = sum(sum(other != p and count > 0 for other, count in adjacency.items())
                                  for p, adjacency in enumerate(expected)) // 2
            require(energy[frame][30] == different_pairs, "energy/labelled-chain distinct-edge mismatch")
            for p in range(n):
                neighbors = integer_row(hashes[frame * n + p])
                counts = integer_row(weights[frame * n + p])
                require(len(neighbors) == len(counts) and bool(neighbors), "adjacency row schema mismatch")
                require(neighbors[0] == p and len(set(neighbors)) == len(neighbors), "adjacency self/duplicate mismatch")
                require(all(0 <= other < n for other in neighbors), "adjacency particle outside range")
                require(all(count >= 0 for count in counts), "negative adjacency weight")
                actual = Counter({other: count for other, count in zip(neighbors, counts) if count})
                require(actual == expected[p], "adjacency/labelled-chain mismatch")
        require(max(bridge_counts) > 0, "network never formed a bridge")
        require(any(row != chains[0] for row in chains[1:]), "recorded chain topology never evolved")
        summary.update(max_bridges=max(bridge_counts), chain_end_count=2 * n * nc,
                       distinct_topology_frames=len({tuple(row) for row in chains}))
    elif parameters["Method"] == "DUMBBELL":
        require(n % 2 == 0 and not arrays["chain"] and not arrays["weight"], "dumbbell topology schema mismatch")
        require(len(arrays["hash"]) == n * len(steps), "dumbbell adjacency frame count mismatch")
        for index, row in enumerate(arrays["hash"]):
            particle = index % n
            require(integer_row(row) == [particle, (particle + n // 2) % n], "dumbbell partner mismatch")
        summary["fixed_dumbbells"] = n // 2
    else:
        require(not arrays["chain"] and not arrays["hash"] and not arrays["weight"],
                "unexpected topology output for non-network case")
    protocol = parameters["MECHANICAL_PERTURBATION"]
    if protocol != "FALSE":
        require(network and parameters["cell_list"] == "FALSE", "unsupported shear reference configuration")
        require(parameters["shear_axis"] == "0" and parameters["shear_grad_axis"] == "1", "unsupported shear axes")
        key = "Wi_tau_C" if protocol == "SIMPLE_SHEAR" else "gamma_0"
        require(float(parameters[key]) != 0, "zero shear in shear fixture")
        marker = "APPLY SIMPLE_SHEAR" if protocol == "SIMPLE_SHEAR" else "APPLY STEP DEFORMATION"
        require(marker in log, "requested shear path was not activated")
        require("measure_minimum_distance_simple_shear_fixed_axis" in log, "shear geometry path was not activated")
        summary.update(shear_protocol=protocol, shear_parameter=float(parameters[key]))
    if parameters["Method"] == "BROWNIAN":
        squared = []
        for previous, current in zip(trajectory, trajectory[1:]):
            for p in range(n):
                total = 0.
                for axis, length in enumerate(lengths):
                    dx = current[1 + 2 * dim * p + axis] - previous[1 + 2 * dim * p + axis]
                    dx -= length * round(dx / length)
                    total += dx * dx
                squared.append(total)
        ratio = sum(squared) / len(squared) / (2 * dim * dt)
        require(abs(ratio - 1) < .06, f"Brownian variance ratio outside bounded expectation: {ratio}")
        summary["brownian_variance_ratio"] = ratio
    return summary


def validate_protocol_effects(directories, parameters):
    """Paired seeded controls prove the frozen drive affects actual coordinates."""
    control = table(Path(directories["network"]) / "case.traj")
    effect = {}
    for name in CASES:
        if not name.startswith("network_"):
            continue
        driven = table(Path(directories[name]) / "case.traj")
        p = parameters[name]
        if p["MECHANICAL_PERTURBATION"] == "STEP_SHEAR":
            gamma = float(p["gamma_0"])
            length = float(p["LBx"])
            for particle in range(int(p["Np"])):
                offset = 1 + 6 * particle
                expected = (control[0][offset] + gamma * control[0][offset + 1]) % length
                require(math.isclose(driven[0][offset], expected, rel_tol=2e-6, abs_tol=1.5e-5),
                        "step-shear initial affine deformation mismatch")
        require(driven[-1][1:] != control[-1][1:], f"{name}: shear has no measured trajectory effect")
        effect[name] = "nonzero drive path and paired-coordinate effect checked"
    return effect


def compare_runs(reference, candidate, *, exact=False):
    """Compare every numeric field; integer topology stays exact in all modes."""
    compared = 0
    for suffix in SUFFIXES:
        aa = table(Path(reference) / f"case.{suffix}", empty=suffix in POLICY["integer_tables"])
        bb = table(Path(candidate) / f"case.{suffix}", empty=suffix in POLICY["integer_tables"])
        require(len(aa) == len(bb), f"{suffix}: frame count differs")
        for frame, (a, b) in enumerate(zip(aa, bb)):
            require(len(a) == len(b), f"{suffix}: row schema differs")
            for column, (x, y) in enumerate(zip(a, b)):
                if suffix == "ener" and column == 5:
                    continue
                integer = suffix in POLICY["integer_tables"]
                agrees = x == y if exact or integer else math.isclose(
                    x, y, rel_tol=POLICY["relative_tolerance"], abs_tol=POLICY["absolute_tolerance"])
                require(agrees, f"{suffix}: scientific difference at row {frame}, column {column}: {x} vs {y}")
                compared += 1
    return {"scientific_values_compared": compared, "comparison": "exact" if exact else "printed-precision tolerance"}


def native_source_manifest():
    paths = [CODE / "Makefile", CODE / "short_test.inp"]
    for directory in (CODE / "lib", CODE / "src"):
        paths.extend(path for path in directory.iterdir() if path.suffix in (".h", ".cpp"))
    return {str(path.relative_to(CODE)): sha256(path) for path in sorted(paths)}


def command_output(command):
    try:
        process = subprocess.run(command, cwd=CODE, capture_output=True, text=True, timeout=10)
        return {"command": command, "returncode": process.returncode,
                "stdout": process.stdout.strip(), "stderr": process.stderr.strip()}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"command": command, "error": str(exc)}


def provenance(executable, compiler):
    return {"created_at_utc": datetime.now(timezone.utc).isoformat(), "scope": SCOPE,
            "executable": {"path": str(executable), "sha256": sha256(executable)},
            "source": {"git_head": command_output(["git", "rev-parse", "HEAD"]),
                       "git_status": command_output(["git", "status", "--porcelain"]),
                       "native_file_sha256": native_source_manifest(),
                       "scope": "All immediate lib/src .cpp/.h plus Makefile and upstream input; dirty working tree retained."},
            "compiler": command_output([compiler, "--version"]),
            "build_flags": {"default_Makefile": "-O2 -g -std=c++17 -Wall -Wextra -fopenmp",
                            "qualification": "Compiler supplied by operator; flags shown are source defaults, not extracted from binary."},
            "gsl": command_output(["pkg-config", "--modversion", "gsl"]),
            "linked_libraries": command_output(["ldd", str(executable)]),
            "environment_overrides": ENV,
            "process_containment": ("outer-runner process group" if os.environ.get("BDMC_REFERENCE_MANAGED_GROUP") == "1"
                                    else "standalone simulator process group"),
            "runtime_environment": {key: value for key, value in os.environ.items()
                                    if key.startswith(("OMP_", "GSL_")) or key in ("LANG", "LC_ALL")},
            "platform": {"system": platform.system(), "release": platform.release(),
                         "machine": platform.machine(), "python": platform.python_version()},
            "harness_sha256": sha256(Path(__file__)),
            "seed_policy": "Input basic_random_seed=100; basic_random_seed_SS=150; GSL mt19937; per-thread seed=base+thread; libc random() seeds initial positions."}


def data_budget(path):
    """Load storage only for a writer, preserving the caller's import path."""
    original_path = sys.path[:]
    try:
        sys.path.append(str(CODE / "python"))
        from bdmc.storage import DataBudget
    finally:
        sys.path[:] = original_path
    return DataBudget(path)


def run_process(command, *, cwd, env, log, timeout=60):
    """Bound a simulator and reap its ordinary process-group descendants.

    Standalone calls own a simulator group. The unified runner explicitly sets
    BDMC_REFERENCE_MANAGED_GROUP=1: this harness must then be its group leader,
    and simulator descendants inherit that group for the runner to clean.
    """
    def limit_file_size():
        import resource
        resource.setrlimit(resource.RLIMIT_FSIZE, (CASE_BYTES, CASE_BYTES))

    managed = env.get("BDMC_REFERENCE_MANAGED_GROUP") == "1"
    require(not managed or os.getpgrp() == os.getpid(), "managed reference harness must be process-group leader")
    process = subprocess.Popen(command, cwd=cwd, env=env, stdout=log,
                               stderr=subprocess.STDOUT, preexec_fn=limit_file_size,
                               start_new_session=not managed)
    previous_handler = None
    install_handler = threading.current_thread() is threading.main_thread()
    if install_handler:
        def terminate(signum, frame):
            raise SystemExit(128 + signum)
        previous_handler = signal.signal(signal.SIGTERM, terminate)
    try:
        process.wait(timeout=timeout)
        return process
    finally:
        # Managed descendants stay in the outer runner's group. Killing that
        # group here would kill this harness; the supervisor owns its cleanup.
        try:
            if managed:
                if process.poll() is None:
                    process.kill()
            else:
                # Also stop descendants left after a successful leader exit,
                # including processes that deliberately ignore SIGTERM.
                os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        finally:
            process.wait()
            if install_handler:
                signal.signal(signal.SIGTERM, previous_handler)


def run_case(executable, output, name, text, budget, behavior="repaired-v2"):
    run = Path(output) / name
    with budget.writer(run, planned_bytes=CASE_BYTES):
        run.mkdir()
        (run / "input.inp").write_text(text)
        env = {**os.environ, **ENV}
        # A regressed executable must not fill the disk while its timeout runs.
        # Seven expected files each have this independent hard cap; the actual
        # bounded inputs remain well below the 4 MB combined-case estimate.
        with (run / "stdout.log").open("x") as log:
            process = run_process([str(executable), str(run / "input.inp")], cwd=run,
                                  env=env, log=log, timeout=60)
        require(process.returncode == 0, f"{name}: exit {process.returncode}; see {run}/stdout.log")
        require(sum(path.stat().st_size for path in run.iterdir() if path.is_file()) < CASE_BYTES,
                f"{name}: exceeded bounded case output estimate")
        summary = validate_run(run, parse_input(text), behavior)
    return run, summary


def load_dataset(dataset):
    dataset = Path(dataset).resolve(strict=True)
    checksums = json.loads((dataset / "CHECKSUMS.json").read_text())
    require(checksums.get("algorithm") == "sha256", "unsupported checksum algorithm")
    actual_files = {str(path.relative_to(dataset)) for path in dataset.rglob("*") if path.is_file()}
    require(not any(path.is_symlink() for path in dataset.rglob("*")), "symlink in frozen dataset")
    require(actual_files == set(checksums["files"]) | {"CHECKSUMS.json"}, "frozen dataset file inventory mismatch")
    for name, digest in checksums["files"].items():
        require(not Path(name).is_absolute() and ".." not in Path(name).parts, "unsafe checksum path")
        require(sha256(dataset / name) == digest, f"frozen fixture checksum mismatch: {name}")
    manifest = json.loads((dataset / "manifest.json").read_text())
    require(manifest.get("schema") == SCHEMA and manifest.get("scope") == SCOPE, "unsupported reference manifest schema/scope")
    behavior = manifest.get("behavioral_revision", "legacy-v1")
    require(behavior in ("legacy-v1", "repaired-v2"), "unsupported behavioral revision")
    require(manifest.get("comparison_policy") == POLICY, "reference comparison policy changed")
    require(set(manifest["cases"]) == set(CASES), "reference case inventory changed")
    directories, parameters = {}, {}
    for name, entry in manifest["cases"].items():
        run = dataset / name
        params = parse_input((run / "input.inp").read_text())
        require(params == entry["parameters"], f"{name}: input/manifest parameter mismatch")
        for key, value in {**COMMON, **CASES[name]}.items():
            require(params.get(key) == value, f"{name}: frozen protocol changed: {key}")
        summary = validate_run(run, params, behavior)
        require(summary == entry["validation"], f"{name}: frozen validation summary differs")
        directories[name], parameters[name] = run, params
    validate_protocol_effects(directories, parameters)
    return manifest


def new_run_directory(parent, budget):
    with budget.writer(parent, planned_bytes=CASE_BYTES * len(CASES)):
        Path(parent).mkdir(parents=True, exist_ok=True)
        return Path(tempfile.mkdtemp(prefix="run-", dir=Path(parent).resolve()))


def create(executable, dataset, run_dir, budget, compiler="g++", behavior="repaired-v2"):
    dataset = Path(dataset).resolve()
    if dataset.exists():
        raise FileExistsError(f"reference dataset already exists; choose a new explicit version: {dataset}")
    require(behavior in ("legacy-v1", "repaired-v2"), "unsupported behavioral revision")
    budget.require_output(dataset, planned_bytes=CORPUS_BYTES)
    output = new_run_directory(run_dir, budget)
    record = provenance(executable, compiler)
    manifest = {"schema": SCHEMA, "scope": SCOPE, "comparison_policy": POLICY,
                "behavioral_revision": behavior,
                "provenance": record, "cases": {}, "raw_run_directory": str(output)}
    directories, parameters = {}, {}
    # Capture complete effective inputs before running so upstream input edits
    # cannot silently produce different cases partway through generation.
    inputs = {name: input_text(overrides) for name, overrides in CASES.items()}
    for name, text in inputs.items():
        run, summary = run_case(executable, output, name, text, budget, behavior)
        directories[name], parameters[name] = run, parse_input(text)
        manifest["cases"][name] = {"parameters": parameters[name], "validation": summary}
        print(f"PASS create: {name} {summary}", flush=True)
    manifest["protocol_checks"] = validate_protocol_effects(directories, parameters)
    require(native_source_manifest() == record["source"]["native_file_sha256"], "native source changed during reference generation")
    require(sha256(executable) == record["executable"]["sha256"], "executable changed during reference generation")
    with budget.writer(dataset, planned_bytes=CORPUS_BYTES):
        dataset.mkdir(parents=True)
        for name, run in directories.items():
            target = dataset / name
            target.mkdir()
            for filename in ("input.inp", "stdout.log", *(f"case.{suffix}" for suffix in SUFFIXES)):
                shutil.copyfile(run / filename, target / filename)
        write_json(dataset / "manifest.json", manifest)
        (dataset / "README.md").write_text(
            "# Frozen BD–MC current-port reference " + dataset.name + "\n\n"
            "## Scope\n" + SCOPE + "\n\n"
            "Behavioural revision: `" + behavior + "`. The format schema remains v1; "
            "behaviour and file-format versions are separate.\n\n"
            "Nine bounded serial cases cover Brownian, dumbbell, repulsive, equilibration, "
            "active transient network and its positive/negative SIMPLE_SHEAR/STEP_SHEAR paths. "
            "Shear uses cubic boxes, x/y axes and no cell list. Initial state and short stochastic "
            "evolution are regression evidence, not long-time physical validation.\n\n"
            "## Integrity and comparison\nCHECKSUMS.json records all other files including manifest.json. "
            "This detects accidental fixture edits; it is not a cryptographic signature. "
            "Verification never regenerates these files. Scientific values use relative tolerance "
            "2e-6 and absolute tolerance 1e-9 because legacy text has limited precision; integer "
            "topology is exact. The only ignored numeric field is zero-based energy column 5 "
            "(wall-clock time). --exact additionally requires exact serial scientific values. "
            "Stdout is retained for diagnostics/provenance, not compared as a numeric observable.\n\n"
            "## Reproduction\nRun tests/reference_suite.py verify with explicit --executable, --dataset "
            "and --run-dir. --multithread-invariants additionally checks two-BD-thread network "
            "runs using invariants; their trajectories are not compared with serial data. "
            "Use an explicitly new version path with create to establish another reviewed baseline.\n")
        write_json(dataset / "CHECKSUMS.json", {"algorithm": "sha256", "files": {
            str(path.relative_to(dataset)): sha256(path) for path in sorted(dataset.rglob("*")) if path.is_file()}})
        require(sum(path.stat().st_size for path in dataset.rglob("*") if path.is_file()) < CORPUS_BYTES,
                "corpus exceeded bounded output estimate")
    load_dataset(dataset)
    with budget.writer(output / "report.json", planned_bytes=250_000):
        write_json(output / "report.json", {"status": "passed", "operation": "create",
                                             "dataset": str(dataset), "scope": SCOPE})
    print(f"Report: {output / 'report.json'}", flush=True)
    return output


def verify(executable, dataset, run_dir, budget, *, exact=False, multithread=False, compiler="g++"):
    manifest = load_dataset(dataset)  # Check fixtures before starting a process.
    behavior = manifest.get("behavioral_revision", "legacy-v1")
    output = new_run_directory(run_dir, budget)
    report = {"operation": "verify", "scope": SCOPE, "dataset": str(Path(dataset).resolve()),
              "behavioral_revision": behavior,
              "dataset_checksums_sha256": sha256(Path(dataset) / "CHECKSUMS.json"),
              "provenance": provenance(executable, compiler), "cases": {}, "multithread_invariants": {}}
    report["capture_identity"] = {
        "native_source_matches": report["provenance"]["source"]["native_file_sha256"] ==
                                 manifest["provenance"]["source"]["native_file_sha256"],
        "executable_matches": report["provenance"]["executable"]["sha256"] ==
                              manifest["provenance"]["executable"]["sha256"],
        "harness_matches": report["provenance"]["harness_sha256"] ==
                           manifest["provenance"]["harness_sha256"],
        "policy": "Identity is reported separately; legitimate source changes still undergo numerical regression."}
    directories, parameters = {}, {}
    try:
        for name in CASES:
            text = (Path(dataset) / name / "input.inp").read_text()
            run, summary = run_case(executable, output, name, text, budget, behavior)
            comparison = compare_runs(Path(dataset) / name, run, exact=exact)
            report["cases"][name] = {"validation": summary, **comparison}
            directories[name], parameters[name] = run, parse_input(text)
            print(f"PASS verify: {name} ({comparison['scientific_values_compared']} values)", flush=True)
        report["protocol_checks"] = validate_protocol_effects(directories, parameters)
        if multithread:
            for name in CASES:
                if not name.startswith("network"):
                    continue
                params = {**parameters[name], "N_THREADS_BD": "2"}
                text = "".join(f"{key}={value}\n" for key, value in sorted(params.items()))
                _, summary = run_case(executable, output, name + "_bd2", text, budget, behavior)
                report["multithread_invariants"][name] = summary
                print(f"PASS invariants only: {name} (BD=2, SS=1)", flush=True)
        report["status"] = "passed"
    except Exception as exc:
        report.update(status="failed", error=str(exc))
        raise
    finally:
        with budget.writer(output / "report.json", planned_bytes=250_000):
            write_json(output / "report.json", report)
        print(f"Report: {output / 'report.json'}", flush=True)
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for operation in ("create", "verify"):
        command = subparsers.add_parser(operation)
        command.add_argument("--executable", required=True, type=Path)
        command.add_argument("--dataset", required=True, type=Path)
        command.add_argument("--run-dir", required=True, type=Path)
        command.add_argument("--budget", type=Path, default=PROJECT / "DATA_BUDGET.json")
        command.add_argument("--compiler", default="g++")
        if operation == "verify":
            command.add_argument("--exact", action="store_true")
            command.add_argument("--multithread-invariants", action="store_true")
        else:
            command.add_argument("--behavior", choices=("legacy-v1", "repaired-v2"), default="repaired-v2",
                                 help="explicit output semantics; never changes an existing corpus")
    args = parser.parse_args(argv)
    executable = args.executable.resolve(strict=True)
    budget = data_budget(args.budget)
    # Keep this bounded batch together. Nested case/report writers remain
    # checked and reentrant in this process. Reacquiring the shared lock per
    # tiny case can otherwise queue behind a different session's long run.
    with budget.writer(args.run_dir, planned_bytes=CASE_BYTES * (len(CASES) + 5)):
        if args.command == "create":
            output = create(executable, args.dataset, args.run_dir, budget, args.compiler, args.behavior)
        else:
            output = verify(executable, args.dataset, args.run_dir, budget,
                            exact=args.exact, multithread=args.multithread_invariants, compiler=args.compiler)
    print(f"PASS {args.command}: {output}", flush=True)


if __name__ == "__main__":
    main()
