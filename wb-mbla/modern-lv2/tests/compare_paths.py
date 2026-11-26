#!/usr/bin/env python3
"""Compare Level 1 and Level 2 intermediate states in disposable builds."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
TRACE_HEADER = r'''
#pragma once
#include <cstdint>
#include <cstring>
#include <iostream>
#include "matrix.h"

struct StateDigest {
  std::uint64_t value = 14695981039346656037ULL;

  void number(double x) {
    static_assert(sizeof(double) == sizeof(std::uint64_t));
    std::uint64_t word;
    std::memcpy(&word, &x, sizeof(word));
    value = (value ^ word) * 1099511628211ULL;
  }

  void matrix(const MATRIX& x) {
    number(x.rows);
    number(x.cols);
    for (long i = 0; i < x.size; ++i) {
      number(x.data[i]);
    }
  }
};

template <class... Matrices>
void trace_state(const char* label, long iteration, const Matrices&... states) {
  StateDigest digest;
  (digest.matrix(states), ...);
  std::clog << label << ' ' << iteration << ' ' << digest.value << '\n';
}

template <class... Values>
void trace_branch(const char* label, long index, Values... values) {
  std::clog << label << ' ' << index << std::hexfloat;
  ((std::clog << ' ' << values), ...);
  std::clog << std::defaultfloat << '\n';
}
'''


def splice(source, anchor, insertion, *, before=False):
    if source.count(anchor) != 1:
        raise ValueError(f"Trace requires one source anchor: {anchor!r}")
    return source.replace(anchor, insertion + anchor if before else anchor + insertion)


def instrument_driver(source, baseline):
    if baseline:
        source = splice(source, "update_Pi_div_DLP_arr(Pi_div_DLP_arr, phiw_set_1,",
                        '\n      trace_state("saved", n, phiw_set_1, Rcake_div_Rmem_arr, '
                        'AUX_h1.shear_rate_nodim_pre, AUX_h1.shear_stress_nodim_pre);\n      ', before=True)
        source = splice(source, "if (IDENT_brake && n>30)",
                        'trace_branch("stop", n, ind_onset_cake, err, err_delta, IDENT_brake);\n      ',
                        before=True)
    else:
        source = splice(source, "update_Pi_div_DLP_arr(osmotic_pressure, previous_wall,",
                        '\n    trace_state("saved", iteration, previous_wall, previous_resistance, '
                        'flow.shear_rate_nodim_pre, flow.shear_stress_nodim_pre);\n    ', before=True)
        source = splice(source, "if (stopping.should_stop(iteration))",
                        'trace_branch("stop", iteration, cake_onset, wall_error, '
                        'stopping.outlet_resistance_change, stopping.latched);\n    ', before=True)
    return '#include "level2_trace.h"\n' + source


def instrument_kernel(source, baseline):
    # The sweep, cake search and final analysis contain some identical calls.
    # Work within the function boundaries so each insertion has a unique anchor.
    start = source.index("double cal_Rcake_div_Rmem(")
    stop = source.index("double process_at_z_modi(", start)
    search = source[start:stop]
    search = splice(search, "long search_steps = 0;",
                    "\n  StateDigest trial_trace;")
    if baseline:
        search = splice(search, 'if (!std::isfinite(f_i)) throw std::runtime_error("nonfinite cake residual");',
                        "\n      trial_trace.number(tmp_nodim_R_2);\n"
                        "      trial_trace.number(DEL_nodim_R_tmp);\n"
                        "      trial_trace.number(f_i);\n")
        search = splice(search, 'if (DEL_nodim_R_tmp <=PARS.TOL_delta_f)',
                        'trial_trace.number(1.);\n          ', before=True)
        search = splice(search, 'if (DEL_nodim_R_tmp <= PARS.TOL_delta_f)',
                        'trial_trace.number(2.);\n          ', before=True)
        search = splice(search, "return Rcake_div_Rmem;",
                        'trace_branch("search", index_z, search_steps, trial_trace.value, Rcake_div_Rmem);\n  ',
                        before=True)
    else:
        search = splice(search, "if (trial > search_upper_limit",
                        "trial_trace.number(trial);\n    trial_trace.number(step);\n"
                        "    trial_trace.number(residual);\n    ", before=True)
        search = splice(search, "if (trial > search_upper_limit || (residual != -1 && residual - previous_residual < 0)) {",
                        "\n      trial_trace.number(1.);")
        search = splice(search, "if (residual * previous_residual < 0) {",
                        "\n      trial_trace.number(2.);")
        search = splice(search, "return previous_trial;",
                        'trace_branch("search", index_z, search_steps, trial_trace.value, previous_trial);\n  ',
                        before=True)
    source = source[:start] + search + source[stop:]

    start = source.index("long phiw_update(")
    stop = source.index("int update_shear_rate_nodim_ZY(", start)
    sweep = source[start:stop]
    if baseline:
        sweep = splice(sweep, "get_vw_div_vw0_arr_with_fouling(AUX_h1.vw_div_vw0_arr, AUX_h1.P_div_DLP_arr, Pi_div_DLP_arr, AUX_h1.Ln_arr, Rcake_div_Rmem_arr, PARS);",
                       '\n  trace_state("hydraulic", n_iter, Pi_div_DLP_arr, AUX_h1.Ln_arr, '
                       'AUX_h1.P_div_DLP_arr, AUX_h1.uz_div_u0_arr, AUX_h1.vw_div_vw0_arr);')
        sweep = splice(sweep, "safety_rule = FALSE;\n      // discontinuity",
                       "const double candidate_before_guard = phiw_set_2(i);\n      ", before=True)
        sweep = splice(sweep, "if (tmp_phiw > phi_f_div_phi_b",
                       'trace_branch("wall", i, candidate_before_guard, safety_rule, tmp_phiw, '
                       'tmp_phiw >= phi_f_div_phi_b);\n      ', before=True)
        sweep = splice(sweep, "return ind_onset_cake;",
                       'trace_state("updated", n_iter, phiw_set_2, Rcake_div_Rmem_new_arr, '
                       'AUX_h1.F_info_arr, AUX_h1.shear_rate_nodim, AUX_h1.shear_stress_nodim, '
                       'AUX_h2.phi_yt_div_phib_arr);\n  ', before=True)
    else:
        sweep = splice(sweep, "get_vw_div_vw0_arr_with_fouling(flow.vw_div_vw0_arr, flow.P_div_DLP_arr, osmotic_pressure,\n"
                       "                                  flow.Ln_arr, previous_resistance, parameters);",
                       '\n  trace_state("hydraulic", iteration, osmotic_pressure, flow.Ln_arr, '
                       'flow.P_div_DLP_arr, flow.uz_div_u0_arr, flow.vw_div_vw0_arr);')
        sweep = splice(sweep, "if (update.relaxed > freeze",
                       'trace_branch("wall", i, candidate, update.used_upstream, update.relaxed, '
                       'update.relaxed >= freeze);\n    ', before=True)
        sweep = splice(sweep, "return cake_onset;",
                       'trace_state("updated", iteration, wall, resistance, flow.F_info_arr, '
                       'flow.shear_rate_nodim, flow.shear_stress_nodim, profiles.phi_yt_div_phib_arr);\n  ',
                       before=True)
    return '#include "level2_trace.h"\n' + source[:start] + sweep + source[stop:]


def build_trace(origin, directory, baseline):
    source = directory / "source"
    source.mkdir(parents=True)
    shutil.copytree(origin / "src", source / "src")
    shutil.copyfile(origin / "CMakeLists.txt", source / "CMakeLists.txt")
    src = source / "src"
    for name, generate in (("calculation_CP_profile.cpp", instrument_driver),
                           ("ported_lib.cpp", instrument_kernel)):
        target = src / name
        target.write_text(generate(target.read_text(), baseline))
    (src / "level2_trace.h").write_text(TRACE_HEADER)
    build = directory / "build"
    with (directory / "build.log").open("wb") as log:
        for argv in (
            ["cmake", "-S", str(source), "-B", str(build), "-DCMAKE_BUILD_TYPE=Release", "-DBUILD_TESTING=OFF"],
            ["cmake", "--build", str(build), "--parallel", "2"],
        ):
            subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=90)
    return build / "mbla"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def run(directory, exe, case, instrumented):
    directory.mkdir()
    config = (ROOT / "examples" / f"{case}.inp").read_text()
    assert config.count("DETAILED_REPORT=FALSE") == 1
    (directory / "case.inp").write_text(config.replace("DETAILED_REPORT=FALSE", "DETAILED_REPORT=TRUE"))
    with (directory / "stdout.log").open("wb") as stdout, (directory / "trace.log").open("wb") as trace:
        subprocess.run([str(exe), "case.inp"], cwd=directory, stdout=stdout, stderr=trace,
                       check=True, timeout=40)
    assert ("MBLA_STATUS=converged" in (directory / "stdout.log").read_text())
    if not instrumented:
        assert not (directory / "trace.log").stat().st_size
    return {name: digest((directory / name).read_bytes())
            for name in ("stdout.log", "result.dat", "result.dat.Rout", "result.dat.Fval")}


def compare_trace(a, b):
    counts = {"saved": 0, "hydraulic": 0, "wall": 0, "search": 0, "updated": 0, "stop": 0}
    guard_hits = 0
    search_steps = 0
    digest_ = hashlib.sha256()
    with a.open("rb") as first, b.open("rb") as second:
        for line, (before, after) in enumerate(itertools.zip_longest(first, second), 1):
            if before != after:
                raise AssertionError(f"Trace divergence at line {line}: {before!r} vs {after!r}")
            fields = before.split()
            event = fields[0].decode()
            counts[event] += 1
            if event == "wall":
                guard_hits += int(fields[3])
            if event == "search":
                search_steps += int(fields[2])
            digest_.update(before)
    assert counts["saved"] == counts["hydraulic"] == counts["updated"] == counts["stop"] > 0
    assert counts["wall"] == 24 * counts["saved"]
    return {**counts, "guard_hits": guard_hits, "search_steps": search_steps}, digest_.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True, help="Unmodified Level 1 source folder")
    parser.add_argument("--baseline-exe", type=Path, required=True)
    parser.add_argument("--candidate-exe", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True, help="New directory for builds and trace evidence")
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir()
    baseline = build_trace(args.baseline.resolve(strict=True), work / "lv1", baseline=True)
    candidate = build_trace(ROOT, work / "lv2", baseline=False)
    results = []
    for case in ("noncake", "cake-mild", "cake"):
        locations = [work / f"{case}-{name}" for name in ("baseline", "candidate", "traced-lv1", "traced-lv2")]
        binaries = [args.baseline_exe.resolve(strict=True), args.candidate_exe.resolve(strict=True),
                    baseline, candidate]
        outputs = [run(directory, exe, case, i >= 2)
                   for i, (directory, exe) in enumerate(zip(locations, binaries, strict=True))]
        assert all(output == outputs[0] for output in outputs[1:]), (case, outputs)
        counts, trace_sha256 = compare_trace(*(path / "trace.log" for path in locations[2:]))
        assert (counts["search"] > 0) == (case != "noncake"), (case, counts)
        row = {"case": case, "counts": counts, "trace_sha256": trace_sha256, "output_sha256": outputs[0]}
        results.append(row)
        print(case, "PASS", counts, flush=True)
    (work / "RESULT.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
