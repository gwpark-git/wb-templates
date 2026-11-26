#!/usr/bin/env python3
"""Compare CPU and optional Apple Metal on identical fixed-volume transport jobs."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess


SEEDS = (57131, 57267, 57901, 58129)
STUDENT_95_DF3 = 3.1824463


def paired_interval(values):
    mean = statistics.mean(values)
    half_width = STUDENT_95_DF3 * statistics.stdev(values) / math.sqrt(len(values))
    return [mean - half_width, mean + half_width]


def compare(rows, candidate, reference):
    old = sorted((r for r in rows if r["backend"] == reference), key=lambda r: r["seed"])
    new = sorted((r for r in rows if r["backend"] == candidate), key=lambda r: r["seed"])
    speed = [a["wall_seconds"] / b["wall_seconds"] for a, b in zip(old, new)]
    deltas = {key: [b[key] - a[key] for a, b in zip(old, new)] for key in
              ("Z", "D_over_D0", "eta_over_eta_s")}
    limits = {
        "D_over_D0": 0.02 * abs(statistics.mean(r["D_over_D0"] for r in old)),
        "eta_over_eta_s": max(2e-4, 0.05 * abs(statistics.mean(r["eta_over_eta_s"] - 1 for r in old))),
    }
    interval = {key: paired_interval(value) for key, value in deltas.items()}
    speed_ci = [math.exp(v) for v in paired_interval([math.log(v) for v in speed])]
    return {
        "reference": reference, "candidate": candidate, "paired_speedups": speed,
        "speedup_median": statistics.median(speed), "speedup_geometric_mean_95_CI": speed_ci,
        "signed_differences": deltas, "paired_mean_difference_95_CI": interval,
        "screening_limits_absolute": limits,
        "material_speed_screen": statistics.median(speed) >= 1.5 and speed_ci[0] > 1,
        "numerical_screen": all(max(map(abs, interval[key])) <= limit for key, limit in limits.items()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu-only", action="store_true")
    parser.add_argument("--extended", action="store_true", help="also compare N=216")
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    capabilities = subprocess.check_output([executable, "--capabilities"], text=True, timeout=15)
    metal = next(line.split("\t", 4) for line in capabilities.splitlines() if line.startswith("metal\t"))
    output = args.output.resolve()
    output.mkdir()
    report = {
        "status": "running", "created_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(), "processor": platform.processor(),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "capabilities": capabilities.splitlines(), "requested": "CPU-only" if args.cpu_only else "CPU and Metal",
        "finite_sampling": "4 discarded and 16 sampled diameter Brownian times, dt=0.002, epsilon/kBT=12",
        "measurement": "native complete property call, including construction, backend setup, integration and estimation",
        "seeds": SEEDS, "cases": [], "runs": [],
    }

    def save():
        (output / "benchmark.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    if not args.cpu_only and metal[3] != "1":
        report.update(status="unavailable", reason=metal[4])
        save()
        return 77

    candidates = ["serial_brute", "serial_cell"]
    if not args.cpu_only:
        candidates.append("metal")
    try:
        particle_counts = (27, 64, 216) if args.extended else (27, 64)
        for n in particle_counts:
            for phi in (0.03, 0.494):
                rows = []
                for index, seed in enumerate(SEEDS):
                    # Paired timing across matching seeds with alternating backend order.
                    for backend in candidates[::(-1 if index % 2 else 1)]:
                        command = [executable, backend, str(n), str(phi), str(seed)]
                        job = subprocess.run(command, text=True, capture_output=True, timeout=300, check=False)
                        if job.returncode:
                            raise RuntimeError(f"{command}: exited {job.returncode}: {job.stderr.strip()}")
                        elapsed, z, d, eta = map(float, next(csv.reader(io.StringIO(job.stdout))))
                        if not all(math.isfinite(value) and value > 0 for value in (elapsed, z, d, eta)):
                            raise ValueError(f"invalid benchmark output from {backend}")
                        row = dict(backend=backend, n=n, phi=phi, seed=seed, wall_seconds=elapsed,
                                   Z=z, D_over_D0=d, eta_over_eta_s=eta)
                        rows.append(row)
                        report["runs"].append(row)
                comparisons = [compare(rows, "serial_cell", "serial_brute")]
                if not args.cpu_only:
                    comparisons.extend(compare(rows, "metal", cpu) for cpu in candidates[:2])
                report["cases"].append({"n": n, "phi": phi, "comparisons": comparisons})
                save()
        report.update(status="completed",
                      metal_speed_and_accuracy_measured=not args.cpu_only)
    except (OSError, RuntimeError, ValueError, subprocess.TimeoutExpired) as error:
        report.update(status="failed", reason=str(error))
    save()
    return 0 if report["status"] == "completed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
