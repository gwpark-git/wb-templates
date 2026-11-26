"""Bounded native MPI ensembles of complete, independently seeded reference runs.

Only the Python caller writes files, under the existing fixed DataBudget. MPI
ranks return compact final observations; they never write trajectories. Each
job's seed, initialization, numerical law and identity are independent of rank
count. This does not implement single-system spatial domain decomposition.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import selectors
import shutil
import signal
import subprocess
import time
from typing import Sequence

from .storage import DataBudget


SCHEMA = "bdmc.mpi-ensemble.v1"
MANIFEST_SCHEMA = "bdmc.mpi-jobs.v1"
COLUMNS = ("run_id", "condition_id", "lineage_id", "model", "seed", "parameter", "count",
           "steps", "dt", "diffusion", "mobility", "box_length", "sigma", "q0", "shear_rate", "warmup")
BACKENDS = ("serial", "openmp", "cuda", "metal", "accelerate")
BACKEND_PRECISION = {"serial": "float64", "openmp": "float64", "cuda": "float64",
                     "metal": "float32", "accelerate": "float64"}
MAX_OUTPUT_BYTES = 16 * 1024 * 1024
MAX_STDOUT_BYTES = 4 * 1024 * 1024 + 4096
MAX_STDERR_BYTES = 1024 * 1024


def _integer(value, name, low, high):
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise ValueError(f"{name} must be an integer in [{low}, {high}]")
    return value


def _number(value, name, low, high, *, positive=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or not low <= value <= high or (positive and value <= 0)):
        raise ValueError(f"{name} must be finite and within [{low}, {high}]")
    return float(value)


@dataclass(frozen=True)
class MPIJob:
    """One explicit condition/seed realization; no automatic seed generation.

    parameter denotes epsilon for repulsive particles and spring_constant for
    dumbbells. count denotes particles or dumbbells respectively. Repulsive
    initialization is the reference regular grid; dumbbells start at zero.
    Dumbbells use simple shear (including zero rate); all jobs use Gaussian
    noise and the same explicit Euler implementation as bdmc.reference.
    """

    run_id: str
    condition_id: str
    model: str
    seed: int
    parameter: float
    count: int = 8
    steps: int = 5
    dt: float = .001
    diffusion: float = .5
    mobility: float = 1.
    box_length: float = 4.
    sigma: float = 1.
    q0: float = 4.
    shear_rate: float = 0.
    warmup: int = 0
    lineage_id: str | None = None

    def __post_init__(self):
        if self.lineage_id is None:
            object.__setattr__(self, "lineage_id", self.run_id)
        for name in ("run_id", "condition_id", "lineage_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_.:-]{1,128}", value):
                raise ValueError(f"{name} requires 1..128 ASCII letters, digits, _, ., :, or -")
        if self.model not in ("repulsive", "hookean", "fene"):
            raise ValueError("model must be repulsive, hookean, or fene")
        _integer(self.seed, "seed", 0, 2**64-1)
        _integer(self.count, "count", 1, 512 if self.model == "repulsive" else 4096)
        _integer(self.steps, "steps", 1, 100000)
        _integer(self.warmup, "warmup", 0, 100000)
        ranges = {"parameter": (0., 1e6), "dt": (0., 1.), "diffusion": (0., 1e6),
                  "mobility": (0., 1e6), "box_length": (0., 1e6), "sigma": (0., 1e6),
                  "q0": (0., 1e6), "shear_rate": (-1e6, 1e6)}
        for name, (low, high) in ranges.items():
            positive = name not in {"diffusion", "shear_rate"} and not (name == "parameter" and self.model == "repulsive")
            object.__setattr__(self, name, _number(getattr(self, name), name, low, high, positive=positive))
        if self.model == "repulsive" and (self.sigma > self.box_length/2 or self.shear_rate != 0):
            raise ValueError("repulsive jobs require sigma<=box_length/2 and no shear")
        if self.model != "repulsive" and self.mobility*self.parameter*self.dt >= 1:
            raise ValueError("dumbbells require mobility*spring_constant*dt<1")


def serialize_jobs(jobs: Sequence[MPIJob]) -> str:
    # MPI redistributes complete explicit realizations. Seed and condition are
    # job properties, so changing rank count must not generate a new ensemble or
    # masquerade as spatial decomposition of one interacting particle system.
    jobs = tuple(jobs)
    if not 1 <= len(jobs) <= 256 or any(not isinstance(job, MPIJob) for job in jobs):
        raise ValueError("MPI batches require 1..256 explicit MPIJob objects")
    if len({job.run_id for job in jobs}) != len(jobs):
        raise ValueError("MPI run_id values must be unique")
    if len({job.seed for job in jobs}) != len(jobs):
        raise ValueError("independent-run MPI simulation seeds must be unique across the manifest")
    work = sum((job.steps+job.warmup)*(job.count**2 if job.model == "repulsive" else 2*job.count)
               for job in jobs)
    if work > 300_000_000:
        raise ValueError("MPI batch exceeds its 300 million update-work bound")
    lines = [MANIFEST_SCHEMA, "\t".join(COLUMNS)]
    for job in jobs:
        record = asdict(job)
        lines.append("\t".join(str(record[name]) for name in COLUMNS))
    manifest = "\n".join(lines)+"\n"
    if len(manifest.encode()) > 1024*1024:
        raise ValueError("MPI job manifest exceeds 1MiB")
    return manifest


def _matches_job(result, job):
    if not isinstance(result, dict) or result.get("status") != "success":
        return False
    for name in ("run_id", "condition_id", "lineage_id", "model", "seed"):
        if result.get(name) != getattr(job, name):
            return False
    expected = {name: getattr(job, name) for name in COLUMNS if name not in
                {"run_id", "condition_id", "lineage_id", "model", "seed"}}
    return (result.get("spec") == expected and result.get("step") == job.steps+job.warmup
            and result.get("time") == (job.steps+job.warmup)*job.dt
            and result.get("observation_lag") == job.steps*job.dt)


class MPIExecutionError(RuntimeError):
    def __init__(self, message, *, report, directory):
        super().__init__(message)
        self.report = report
        self.directory = Path(directory)


def _terminate_group(process):
    # Children inherit the new session. Killing only mpiexec can leave local
    # workers running after a timeout; terminate its complete local group.
    for signum in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(process.pid, signum)
        except ProcessLookupError:
            return
        if signum == signal.SIGTERM:
            try:
                process.wait(timeout=.5)
            except subprocess.TimeoutExpired:
                pass
            # Kill surviving group members even when the launcher exited.
    process.wait(timeout=5)


def _bounded_process(command, *, directory, environment, timeout):
    stdout, stderr = bytearray(), bytearray()
    problem = None
    started = time.monotonic()
    process = subprocess.Popen(command, cwd=directory, env=environment,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               start_new_session=True)
    try:
        with selectors.DefaultSelector() as selector:
            selector.register(process.stdout, selectors.EVENT_READ, (stdout, MAX_STDOUT_BYTES))
            selector.register(process.stderr, selectors.EVENT_READ, (stderr, MAX_STDERR_BYTES))
            while selector.get_map():
                remaining = timeout-(time.monotonic()-started)
                if remaining <= 0:
                    problem = f"MPI batch exceeded {timeout:g}-second timeout; process group terminated"
                    break
                for key, _ in selector.select(min(.2, remaining)):
                    chunk = os.read(key.fileobj.fileno(), 65536)
                    if not chunk:
                        selector.unregister(key.fileobj)
                        continue
                    buffer, bound = key.data
                    if len(buffer)+len(chunk) > bound:
                        buffer.extend(chunk[:bound-len(buffer)])
                        problem = "MPI launcher output exceeded its declared bound; process group terminated"
                        break
                    buffer.extend(chunk)
                if problem:
                    break
        if problem:
            _terminate_group(process)
        else:
            try:
                process.wait(timeout=max(.01, timeout-(time.monotonic()-started)))
            except subprocess.TimeoutExpired:
                problem = f"MPI batch exceeded {timeout:g}-second timeout; process group terminated"
                _terminate_group(process)
    except BaseException:
        _terminate_group(process)
        raise
    finally:
        process.stdout.close()
        process.stderr.close()
    return process.returncode, bytes(stdout), bytes(stderr), problem


def run_mpi_ensemble(jobs: Sequence[MPIJob], *, executable, output, budget: DataBudget,
                     ranks=2, backend="serial", threads=1, launcher="mpiexec",
                     launcher_args: Sequence[str] = (), timeout=60.) -> dict:
    """Launch actual native MPI ranks and persist a reviewable bounded record.

    output must be a fresh directory beneath a registered fixed-baseline data
    root. The caller can supply e.g. Open MPI --host/--hostfile/--oversubscribe
    or a scheduler-specific launcher without shell evaluation. Rank0 alone
    reads the input; remote nodes need the same native executable/libraries.
    Any failed job, timeout or launcher failure raises MPIExecutionError with
    the saved report. Nothing is retried or silently excluded.
    """
    jobs = tuple(jobs)
    manifest = serialize_jobs(jobs)
    _integer(ranks, "ranks", 1, 64)
    _integer(threads, "threads", 1, 256)
    if ranks*threads > 256:
        raise ValueError("MPI ranks*threads must not exceed 256")
    if backend not in BACKENDS:
        raise ValueError(f"backend must be one of {BACKENDS}")
    _number(timeout, "timeout", 0., 3600., positive=True)
    if not isinstance(budget, DataBudget):
        raise TypeError("an existing fixed-baseline DataBudget is required")
    executable = Path(executable).resolve(strict=True)
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError("MPI executable must be an executable regular file")
    resolved_launcher = shutil.which(os.fspath(launcher))
    if resolved_launcher is None:
        raise FileNotFoundError(f"MPI launcher is unavailable: {launcher}")
    if isinstance(launcher_args, (str, bytes)) or any(not isinstance(arg, str) or "\0" in arg for arg in launcher_args):
        raise ValueError("launcher_args must be a sequence of literal argument strings")
    output = Path(output).resolve()
    command = [resolved_launcher, *launcher_args, "-n", str(ranks), str(executable),
               "--jobs-file", str(output/"jobs.tsv"), "--backend", backend, "--threads", str(threads)]
    report = {"schema": SCHEMA, "status": "running", "jobs": [asdict(job) for job in jobs],
              "created_at_utc": datetime.now(timezone.utc).isoformat(),
              "executable": str(executable), "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
              "manifest_sha256": hashlib.sha256(manifest.encode()).hexdigest(), "command": command,
              "limits": {"max_jobs": 256, "update_work": 300_000_000, "output_bytes": MAX_OUTPUT_BYTES,
                         "timeout_seconds": float(timeout)},
              "failure_policy": "raise on unsuccessful batch; no retries, dropped runs or replacement seeds"}
    with budget.writer(output, planned_bytes=MAX_OUTPUT_BYTES):
        output.mkdir(parents=True, exist_ok=False)
        (output/"jobs.tsv").write_text(manifest)
        (output/"report.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    environment = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE",
                   "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if hasattr(os, "geteuid") and os.geteuid() == 0:
        environment.update(OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    started = time.monotonic()
    stdout, stderr = b"", b""
    error = None
    try:
        returncode, stdout, stderr, error = _bounded_process(
            command, directory=output, environment=environment, timeout=float(timeout))
        report["returncode"] = returncode
        if error is None:
            try:
                native = json.loads(stdout)
            except (ValueError, UnicodeError) as parse_error:
                error = f"MPI launcher did not return a valid JSON report: {parse_error}"
            else:
                report["native"] = native
                if not isinstance(native, dict) or native.get("schema") != SCHEMA:
                    error = "MPI executable returned an unsupported report schema"
                elif returncode != 0 or native.get("status") != "success":
                    error = "MPI execution failed; inspect native job errors and stderr.log"
                elif [row.get("run_id") for row in native.get("results", [])] != [job.run_id for job in jobs]:
                    error = "MPI results do not preserve the explicit manifest order and identities"
                elif any(row.get("status") != "success" for row in native["results"]):
                    error = "MPI executable returned an unsuccessful job"
                else:
                    for job, row in zip(jobs, native["results"]):
                        if not _matches_job(row, job):
                            error = "MPI result seed, specification, time or source identity differs from manifest"
                            break
    except Exception as exception:
        error = f"MPI launch failed: {type(exception).__name__}: {exception}"
    finally:
        report["status"] = "failed" if error else "success"
        if error:
            report["error"] = error
        report["elapsed_seconds"] = time.monotonic()-started
        report["storage_after_compute"] = budget.check()
        encoded = (json.dumps(report, indent=2, allow_nan=False)+"\n").encode()
        planned = len(encoded)+len(stdout)+len(stderr)+3*4096
        if planned > MAX_OUTPUT_BYTES:
            raise RuntimeError("MPI saved output exceeds its predeclared 16MiB bound")
        with budget.writer(output, planned_bytes=planned):
            (output/"stdout.json").write_bytes(stdout)
            (output/"stderr.log").write_bytes(stderr)
            (output/"report.json").write_bytes(encoded)
    if error:
        raise MPIExecutionError(error, report=report, directory=output)
    return report


def mpi_dataset(report: dict, *, model: str, target_names: Sequence[str],
                feature_names: Sequence[str] = ("parameter", "observation_lag")):
    """Adapt final scalar observations into the existing grouped ML Dataset.

    Select one physical model explicitly. Different models retain different
    observable contracts. Seeds/identifiers are provenance, never predictive
    features. Condition and lineage grouping remain available to grouped_split.
    The adapter accepts successful verified Python launcher reports only.
    """
    import numpy as np
    from .ml import Dataset, RegressionTask

    if not isinstance(report, dict) or report.get("schema") != SCHEMA or report.get("status") != "success":
        raise ValueError("MPI dataset requires a successful launcher report")
    native = report.get("native", {})
    if native.get("schema") != SCHEMA or native.get("status") != "success":
        raise ValueError("MPI dataset has no successful native report")
    execution = native.get("execution", {})
    backend, backend_threads = execution.get("backend"), execution.get("threads_per_rank")
    if backend not in BACKEND_PRECISION:
        raise ValueError("MPI dataset is missing its explicit execution backend")
    _integer(backend_threads, "threads_per_rank", 1, 256)
    if backend != "openmp" and backend_threads != 1:
        raise ValueError("MPI result threads>1 require the OpenMP backend")
    if model not in ("repulsive", "hookean", "fene"):
        raise ValueError("select repulsive, hookean, or fene explicitly")
    features, targets = tuple(feature_names), tuple(target_names)
    physical = {"parameter", "count", "steps", "dt", "diffusion", "mobility", "box_length", "sigma", "q0", "shear_rate", "warmup", "time", "observation_lag"}
    units = {"potential_energy_per_particle": "energy", "finite_lag_msd": "length_squared",
             "connector_second_moment": "length_squared", "polymer_shear_stress_per_dumbbell": "energy",
             "first_normal_stress_per_dumbbell": "energy", "mean_energy": "energy"}
    supported = ({"potential_energy_per_particle", "finite_lag_msd"} if model == "repulsive" else
                 {"connector_second_moment", "polymer_shear_stress_per_dumbbell", "first_normal_stress_per_dumbbell", "mean_energy"})
    if not features or len(set(features)) != len(features) or any(name not in physical for name in features):
        raise ValueError("features must be unique explicit physical fields; seeds/identities are provenance")
    if not targets or len(set(targets)) != len(targets) or any(name not in supported for name in targets):
        raise ValueError("targets must be distinct scalar observable names for the selected model")
    if (not report.get("manifest_sha256") or not report.get("executable_sha256")
            or not isinstance(report.get("jobs"), list)):
        raise ValueError("MPI report is missing manifest/executable source provenance")
    expected = [MPIJob(**job) for job in report["jobs"]]
    results = native.get("results", [])
    if len(results) != len(expected):
        raise ValueError("MPI result count differs from the explicit manifest")
    if hashlib.sha256(serialize_jobs(expected).encode()).hexdigest() != report["manifest_sha256"]:
        raise ValueError("MPI manifest identity differs from its saved SHA256")
    rows, labels, groups, conditions, lineages, provenance = [], [], [], [], [], []
    for job, result in zip(expected, results):
        if not _matches_job(result, job):
            raise ValueError("MPI result has an unsuccessful or mismatched physical/source identity")
        if job.model != model:
            continue
        values = {**result["spec"], "time": result["time"], "observation_lag": result["observation_lag"]}
        feature_row = [values[name] for name in features]
        target_row = [result["summary"][name] for name in targets]
        if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
               for value in (*feature_row, *target_row)):
            raise ValueError("MPI dataset requires finite scalar features/observables")
        rows.append(feature_row); labels.append(target_row)
        groups.append(job.run_id); conditions.append(job.condition_id); lineages.append(job.lineage_id)
        provenance.append({**asdict(job), "initial_state_fnv1a64": result["initial_state_fnv1a64"],
                           "origin_state_fnv1a64": result["origin_state_fnv1a64"],
                           "final_state_fnv1a64": result["final_state_fnv1a64"]})
    if not rows:
        raise ValueError("MPI report contains no rows for the selected model")
    metadata = {"schema": "bdmc.mpi-dataset.v1", "model": model, "runs": provenance,
                "manifest_sha256": report["manifest_sha256"], "executable_sha256": report["executable_sha256"],
                "execution": {"backend": backend, "threads_per_rank": backend_threads,
                              "arithmetic_precision": BACKEND_PRECISION[backend],
                              "rank_count_note": "operational rank count remains in launcher report; scientific rows preserve rank-invariant identities"},
                "initialization": "native regular grid for repulsive; zero connectors for dumbbells",
                "noise": "Gaussian, independent explicitly seeded per job; rank count does not alter seed/order",
                "sampling": "one final observation per run after the declared warmup; not an equilibrium guarantee",
                "observable_definitions": {
                    "potential_energy_per_particle": "final pair potential energy divided by particle count",
                    "finite_lag_msd": "mean squared unwrapped displacement from post-warmup state",
                    "connector_second_moment": "trace of final mean Q outer Q",
                    "polymer_shear_stress_per_dumbbell": "positive Cauchy mean Qx Fy (off diagonal)",
                    "first_normal_stress_per_dumbbell": "positive Cauchy stress xx minus yy",
                    "mean_energy": "final spring energy averaged over dumbbells"},
                "limitations": "finite realization response labels; no single-system decomposition or scientific validation implied"}
    return Dataset(np.asarray(rows), np.asarray(labels), np.asarray(groups),
                   feature_names=features, target_names=targets,
                   group_keys={"condition": np.asarray(conditions), "continuation": np.asarray(lineages)},
                   task=RegressionTask(name=f"{model}_mpi_reference", target_kind="realization",
                                       target_units=tuple(units[name] for name in targets)), metadata=metadata)
