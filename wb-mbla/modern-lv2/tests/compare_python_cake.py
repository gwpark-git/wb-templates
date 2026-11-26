#!/usr/bin/env python3
"""Compare and time the Python cake search on exact inputs captured from C++."""

import argparse
from dataclasses import replace
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from mbla_cake import CakeSearchState, search_resistance


# Instrument only an isolated copy. Capture a bounded sequence of searches
# starting early and continuing through the selected fixture's run.
SETUP = r'''
  static long mbla_call = 0;
  ++mbla_call;
  const bool mbla_capture = mbla_call <= 6 || (mbla_call % 139 == 0 && mbla_call < 1800);
  static const bool mbla_trace = std::getenv("MBLA_COMPARE_TRIALS") != nullptr;
  std::vector<std::array<double, 3>> mbla_trials;
  bool mbla_warning = false;
  const auto mbla_start = mbla_capture ? std::chrono::steady_clock::now()
                                       : std::chrono::steady_clock::time_point();
'''

END = r'''
  if (mbla_capture) {
    const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::steady_clock::now() - mbla_start).count();
    std::clog << "sample " << mbla_call << ' ' << elapsed << ' ' << search_steps << ' '
              << mbla_warning << std::hexfloat << ' ' << previous_trial << ' '
              << uz_div_u0 << ' ' << PARS.Lp << ' '
              << (pressure_ratio * PARS.DLP - PARS.Pper - osmotic_ratio * PARS.DLP) << ' '
              << PARS.vw0 << ' ' << PARS.epsilon_d << ' ' << area_factor << ' '
              << inlet_flux << ' ' << PARS.phi_bulk << ' ' << PARS.phi_freeze << ' '
              << PARS.TOL_delta_f << '\n';
    for (long j = 0; j < PARS.Ny; ++j) {
      std::clog << "radial " << yt_arr(j) << ' ' << profiles.Ieta_yt_arr(index_z, j)
                << ' ' << profiles.ID_yt_arr(index_z, j) << '\n';
    }
    for (const auto& row : mbla_trials) {
      std::clog << "trial " << row[0] << ' ' << row[1] << ' ' << row[2] << '\n';
    }
    std::clog << "end\n" << std::defaultfloat;
  }
'''


def insert(text, marker, addition):
    if text.count(marker) != 1:
        raise ValueError(f"Expected one C++ capture anchor: {marker!r}")
    return text.replace(marker, addition + marker)


def build_capture(work):
    source = work / "capture-source"
    shutil.copytree(ROOT / "src", source / "src")
    shutil.copyfile(ROOT / "CMakeLists.txt", source / "CMakeLists.txt")
    file = source / "src/ported_lib.cpp"
    code = file.read_text()
    start = code.index("double cal_Rcake_div_Rmem(")
    stop = code.index("double process_at_z_modi(", start)
    kernel = code[start:stop]
    kernel = insert(kernel, "  while (true) {", SETUP)
    kernel = insert(kernel, "    if (trial > search_upper_limit",
                    "    if (mbla_capture && mbla_trace) mbla_trials.push_back({trial, step, residual});\n")
    kernel = insert(kernel, '        printf("Warning:',
                    "        mbla_warning = true;\n")
    kernel = insert(kernel, "  return previous_trial;", END)
    file.write_text("#include <array>\n#include <chrono>\n#include <cstdlib>\n#include <vector>\n"
                    + code[:start] + kernel + code[stop:])
    build = work / "capture-build"
    with (work / "capture-build.log").open("wb") as log:
        subprocess.run(["cmake", "-S", str(source), "-B", str(build), "-DCMAKE_BUILD_TYPE=Release",
                        "-DBUILD_TESTING=OFF"], stdout=log, stderr=subprocess.STDOUT, check=True, timeout=60)
        subprocess.run(["cmake", "--build", str(build), "--parallel", "2"],
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=90)
    return build / "mbla"


def run(executable, directory, config, capture_trials, expected_status):
    directory.mkdir()
    (directory / "case.inp").write_text(config)
    environment = os.environ.copy()
    environment.pop("MBLA_COMPARE_TRIALS", None)
    if capture_trials:
        environment["MBLA_COMPARE_TRIALS"] = "1"
    with (directory / "solver.log").open("wb") as log, (directory / "samples.txt").open("wb") as samples:
        process = subprocess.run([str(executable), "case.inp"], cwd=directory, env=environment,
                                 stdout=log, stderr=samples, timeout=60)
    if process.returncode != expected_status:
        raise RuntimeError(f"{directory.name}: expected {expected_status}, observed {process.returncode}")
    return {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (directory / "solver.log", directory / "result.dat")}


def read_samples(path):
    with path.open() as file:
        for line in file:
            label, sample_id, elapsed_ns, steps, warning, resistance, *scalars = line.split()
            assert label == "sample"
            radial = []
            trials = []
            for record in file:
                kind, *values = record.split()
                if kind == "end":
                    break
                row = tuple(map(float.fromhex, values))
                assert len(row) == 3
                if kind == "radial":
                    radial.append(row)
                elif kind == "trial":
                    trials.append(row)
                else:
                    raise ValueError(f"Unknown capture record: {kind}")
            arrays = tuple(zip(*radial, strict=True))
            state = CakeSearchState(*arrays, *map(float.fromhex, scalars))
            yield (int(sample_id), int(elapsed_ns), int(steps), bool(int(warning)),
                   float.fromhex(resistance), state, tuple(trials))


def compare_candidate(sample, backend):
    _, _, steps, warning, resistance, state, trials = sample
    result = search_resistance(state, backend=backend, trace=True)
    assert (result.resistance, result.steps, result.termination) == (
        resistance, steps, "warning" if warning else "bracketed"), (sample[0], result)
    assert len(result.trial_trace) == len(trials) == steps
    max_residual_difference = 0.
    for expected, found in zip(trials, result.trial_trace, strict=True):
        assert expected[:2] == found[:2], (sample[0], expected, found)
        max_residual_difference = max(max_residual_difference, abs(expected[2] - found[2]))
        # exp() may differ by one ULP across toolchains. The measured path
        # and resistance are held exact even if a residual is rounded.
        assert math.isclose(expected[2], found[2], rel_tol=1e-13, abs_tol=1e-14), (sample[0], expected, found)
    return max_residual_difference


def elapsed_seconds(samples, backend):
    trials = []
    for _ in range(3):
        start = time.perf_counter_ns()
        for sample in samples:
            search_resistance(sample[5], backend=backend)
        trials.append((time.perf_counter_ns() - start) / 1e9)
    return statistics.median(trials)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True, help="Ordinary strict-FP C++ executable")
    parser.add_argument("--work", type=Path, required=True, help="New directory for build and measurements")
    parser.add_argument("--numba", action="store_true", help="Include the optional CPU compiler")
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir()
    executable = args.exe.resolve(strict=True)
    instrumented = build_capture(work)

    cases = (
        ("mild", "cake-mild", 120, None, 0),
        ("demanding", "cake", 120, None, 0),
        ("wide-mild", "cake-mild", 400, 6, 3),
    )
    results = []
    warmup_seconds = None
    for name, fixture, nr, n_iter, status in cases:
        config = (ROOT / "examples" / f"{fixture}.inp").read_text().replace("Nr=120", f"Nr={nr}")
        # Capture the same inlet edition as the archived inner-search benchmark.
        config += "inlet_velocity_order=1\n"
        if n_iter:
            lines = [f"N_iter={n_iter}" if line.startswith("N_iter=") else line
                     for line in config.splitlines()]
            config = "\n".join(lines) + "\n"
        expected = run(executable, work / (name + "-original"), config, False, status)
        traced_dir = work / (name + "-traced")
        assert run(instrumented, traced_dir, config, True, status) == expected
        traced = list(read_samples(traced_dir / "samples.txt"))
        assert traced
        cpp_seconds = []
        for repeat in range(3):
            measured_dir = work / f"{name}-measured-{repeat}"
            assert run(instrumented, measured_dir, config, False, status) == expected
            measured = list(read_samples(measured_dir / "samples.txt"))
            assert [(s[0], s[2:6]) for s in traced] == [(s[0], s[2:6]) for s in measured]
            cpp_seconds.append(sum(s[1] for s in measured) / 1e9)
        summary = {
            "case": name, "radial_nodes": nr, "searches": len(traced),
            "trial_count": sum(s[2] for s in traced), "warning_exits": sum(s[3] for s in traced),
            "cpp_seconds": statistics.median(cpp_seconds), "solver_sha256": expected,
        }
        backends = ["python"]
        if args.numba:
            backends.append("numba")
            if warmup_seconds is None:
                start = time.perf_counter()
                search_resistance(traced[0][5], backend="numba")
                warmup_seconds = time.perf_counter() - start
        for backend in backends:
            summary[backend] = {
                "max_residual_abs_difference": max(compare_candidate(s, backend) for s in traced),
                "median_seconds": elapsed_seconds(traced, backend),
            }
        results.append(summary)
        print(json.dumps(summary), flush=True)

    # Validate the Python boundary independently of C++ replay.
    warning_case = CakeSearchState((0., .1, .3), (0., .1, 1.), (0., .5, 1.),
                                   1., 1e-9, -1., 1e-7, .1, 1., 5., .001, .494, .01)
    for backend in backends:
        warning = search_resistance(warning_case, backend=backend, trace=True)
        assert (warning.resistance, warning.steps, warning.termination) == (.01, 2, "warning")
        assert warning.trial_trace[1][2] < warning.trial_trace[0][2]
    invalid = replace(warning_case, wall_distance=(0., 0., 1.))
    try:
        search_resistance(invalid)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid radial grid was accepted")
    report = {"scope": "Inner searches on captured fixed-point states, not whole-solver performance.",
              "cpp_executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
              "numba_first_call_seconds": warmup_seconds,
              "synthetic_warning_branch": "PASS", "invalid_radial_grid": "PASS", "cases": results}
    (work / "RESULT.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
