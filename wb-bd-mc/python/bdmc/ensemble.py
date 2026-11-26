"""Bounded independent-run CPU ensembles; transport JSON specs, never native state.

This is a local spawn-process execution seam, not a distributed/GPU backend.
Workers produce compact named observables and never write files themselves.
"""
from __future__ import annotations

from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
import hashlib
import json
import math
import multiprocessing
import os
from pathlib import Path
import time
import traceback
from types import MappingProxyType
from typing import Any, Mapping, Sequence

OBSERVABLES = {
    "msd": {"definition": "particle mean of squared unwrapped displacement from this job's initial state",
            "units": "length_squared"},
    "connector_second_moment": {"definition": "fixed-pair mean squared minimum-image connector length at observation time",
                                "units": "length_squared"},
    "potential_energy": {"definition": "current core total potential energy; not per-particle energy", "units": "core_energy"},
    "mean_unwrapped_position": {"definition": "particle mean unwrapped position per coordinate", "units": "length"},
}
_THREAD_CONTROL = None


def _integer(value: Any, name: str, *, zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < (0 if zero else 1):
        raise ValueError(f"{name} must be a {'nonnegative' if zero else 'positive'} integer")
    return value


def _identity(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > 256:
        raise ValueError(f"{name} must be a nonempty string with at most 256 characters")
    return value


@dataclass(frozen=True)
class JobSpec:
    # run_id names an execution, condition_id its physical setting, and lineage_id
    # its shared continuation history. Callers supply those scientific relations;
    # distinct strings or seeds alone are not proof of independent realizations.
    run_id: str
    spec: Any
    steps: int
    stride: int
    observables: tuple[str, ...]
    condition_id: str
    lineage_id: str
    initial_positions: tuple[tuple[float, ...], ...] | None = None
    particle_ids: tuple[int, ...] | None = None
    include_initial: bool = True

    def __post_init__(self) -> None:
        from . import SimulationSpec

        for name in ("run_id", "condition_id", "lineage_id"):
            _identity(getattr(self, name), name)
        _integer(self.steps, "steps", zero=True)
        _integer(self.stride, "stride")
        if not isinstance(self.include_initial, bool) or (self.steps == 0 and not self.include_initial):
            raise ValueError("include_initial must be boolean and the job must record at least one frame")
        source = self.spec.to_dict() if hasattr(self.spec, "to_dict") else dict(self.spec)
        for name in ("n_particles", "dimension", "seed"):
            if name in source:
                _integer(source[name], name, zero=name == "seed")
        # Core validation normalises all defaults. Store a detached immutable snapshot.
        normalised = SimulationSpec(**source).to_dict()
        normalised["box_lengths"] = tuple(normalised["box_lengths"])
        object.__setattr__(self, "spec", MappingProxyType(normalised))
        names = tuple(self.observables)
        if not names or len(set(names)) != len(names) or any(name not in OBSERVABLES for name in names):
            raise ValueError("observables must be distinct supported names")
        if "connector_second_moment" in names and normalised["model"] != "dumbbell":
            raise ValueError("connector_second_moment requires fixed dumbbell pairs")
        object.__setattr__(self, "observables", names)
        if self.initial_positions is not None:
            positions = tuple(tuple(float(value) for value in row) for row in self.initial_positions)
            if (len(positions) != normalised["n_particles"]
                    or any(len(row) != normalised["dimension"] for row in positions)
                    or any(not math.isfinite(value) for row in positions for value in row)):
                raise ValueError("initial_positions must be finite and match particle/dimension shape")
            object.__setattr__(self, "initial_positions", positions)
        if self.particle_ids is not None:
            ids = tuple(_integer(value, "particle ID", zero=True) for value in self.particle_ids)
            if (len(ids) != normalised["n_particles"] or len(set(ids)) != len(ids)
                    or any(value > 2**64 - 1 for value in ids)):
                raise ValueError("particle IDs must be unique uint64 values, one per particle")
            object.__setattr__(self, "particle_ids", ids)

    def to_dict(self) -> dict[str, Any]:
        spec = dict(self.spec)
        spec["box_lengths"] = list(spec["box_lengths"])
        return {"run_id": self.run_id, "condition_id": self.condition_id, "lineage_id": self.lineage_id,
                "spec": spec, "steps": self.steps, "stride": self.stride,
                "observables": list(self.observables), "include_initial": self.include_initial,
                "initial_positions": self.initial_positions, "particle_ids": self.particle_ids}


EnsembleJob = JobSpec


@dataclass(frozen=True)
class EnsembleLimits:
    max_jobs: int = 32
    max_workers: int = 4
    max_timeout_seconds: int = 300
    max_particles: int = 10000
    max_steps: int = 100000
    max_frames_per_job: int = 1001
    max_particle_steps: int = 10000000
    max_estimated_output_bytes: int = 8000000
    max_estimated_worker_array_bytes: int = 32000000

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            _integer(value, name)


def _estimates(jobs: Sequence[JobSpec], limits: EnsembleLimits) -> dict[str, Any]:
    if not jobs or len(jobs) > limits.max_jobs:
        raise ValueError("ensemble must contain 1..max_jobs jobs")
    estimates = []
    for job in jobs:
        n, d = job.spec["n_particles"], job.spec["dimension"]
        frames = (job.steps + job.stride - 1) // job.stride + int(job.include_initial)
        # Includes conservative numeric JSON expansion, metadata and one-frame arrays.
        width = sum(d if name == "mean_unwrapped_position" else 1 for name in job.observables)
        output = 65536 + len(json.dumps(job.to_dict(), allow_nan=False).encode()) * 8 + frames * (2 + width) * 80
        array_bytes = n * d * 8 * 16 + frames * (2 + width) * 8
        if (n > limits.max_particles or job.steps > limits.max_steps
                or frames > limits.max_frames_per_job or n * job.steps > limits.max_particle_steps
                or array_bytes > limits.max_estimated_worker_array_bytes):
            raise ValueError(f"job {job.run_id} exceeds advertised particle/step/frame/array bounds")
        estimates.append({"run_id": job.run_id, "frames": frames, "estimated_output_bytes": output,
                          "estimated_worker_array_bytes": array_bytes, "particle_steps": n * job.steps})
    total = 131072 + sum(item["estimated_output_bytes"] for item in estimates)
    if total > limits.max_estimated_output_bytes:
        raise ValueError("ensemble exceeds max_estimated_output_bytes")
    return {"estimated_output_bytes": total, "jobs": estimates,
            "memory_estimate_scope": "numeric native/NumPy arrays and compact summaries; excludes interpreter/shared-library RSS"}


def _worker_initialise() -> None:
    # This module intentionally has no top-level NumPy/core import: these limits
    # are set before the spawned worker imports native numerical libraries.
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"


def _base_result(job: Mapping[str, Any], expected_hash: str) -> dict[str, Any]:
    return {"run_id": job["run_id"], "condition_id": job["condition_id"], "lineage_id": job["lineage_id"],
            "seed": job["spec"]["seed"], "spec": job["spec"], "job": dict(job),
            "status": "pending", "worker_pid": None, "run_seconds": None,
            "core_binary_sha256": None, "expected_core_binary_sha256": expected_hash,
            "core_hash_verified": False, "summary": {"steps": [], "times": [],
                "values": {name: [] for name in job["observables"]}},
            "observable_definitions": {name: OBSERVABLES[name] for name in job["observables"]},
            "error": None}


def _execute_payload(payload: str) -> dict[str, Any]:
    # Spawn workers construct their own simulator from a value specification.
    # Sharing/pickling a live native object would couple locks, RNG ownership and
    # process lifetime; compact observations are the explicit return boundary.
    envelope = json.loads(payload)
    job = envelope["job"]
    result = _base_result(job, envelope["expected_hash"])
    result["worker_pid"] = os.getpid()
    started = time.perf_counter()
    try:
        if time.monotonic() >= envelope["deadline"]:
            raise TimeoutError("ensemble deadline passed before this job started")
        import numpy as np
        from . import Simulation, SimulationSpec, _core
        digest = hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest()
        result["core_binary_sha256"] = digest
        if digest != envelope["expected_hash"]:
            raise RuntimeError("core binary changed between ensemble preflight and worker import")
        result["core_hash_verified"] = True
        np.seterr(over="raise", invalid="raise", divide="raise")

        global _THREAD_CONTROL
        try:
            from threadpoolctl import threadpool_info, threadpool_limits
            _THREAD_CONTROL = threadpool_limits(limits=1)
            result["native_thread_pools"] = threadpool_info()
            result["native_thread_limit_verification"] = "threadpoolctl runtime control and inspection"
            if any(pool.get("num_threads", 1) > 1 for pool in result["native_thread_pools"]):
                raise RuntimeError("a detected native numerical pool exceeds the one-thread worker limit")
        except ImportError:
            result["native_thread_pools"] = None
            result["native_thread_limit_verification"] = "environment requested; runtime libraries not inspected"
        spec = SimulationSpec(**job["spec"])
        simulation = Simulation(spec, initial_positions=job["initial_positions"], particle_ids=job["particle_ids"])
        initial = simulation.observe()
        origin = np.array(initial["unwrapped_positions"], copy=True)
        ids = initial["particle_ids"].tolist()
        id_to_row = {particle_id: i for i, particle_id in enumerate(ids)}
        bonds = initial["fixed_bonds"].tolist()
        pair_rows = [(id_to_row[a], id_to_row[b]) for a, b in bonds]
        box = np.asarray(spec.box_lengths)
        result["metadata"] = dict(simulation.metadata)
        result["native_thread_limit"] = 1

        def record(frame):
            values = {}
            for name in job["observables"]:
                if name == "msd":
                    value = float(np.square(frame["unwrapped_positions"] - origin).sum(axis=1).mean())
                elif name == "connector_second_moment":
                    q = np.array([frame["positions"][b] - frame["positions"][a] for a, b in pair_rows])
                    q = np.remainder(q + box / 2, box) - box / 2
                    value = float(np.square(q).sum(axis=1).mean())
                elif name == "potential_energy":
                    value = float(frame["potential_energy"])
                else:
                    value = frame["unwrapped_positions"].mean(axis=0).tolist()
                if not np.isfinite(value).all():
                    raise ValueError(f"nonfinite observable {name}")
                values[name] = value
            # Commit a frame only when every requested observable is finite.
            result["summary"]["steps"].append(int(frame["step"]))
            result["summary"]["times"].append(float(frame["time"]))
            for name, value in values.items():
                result["summary"]["values"][name].append(value)

        if job["include_initial"]:
            record(initial)
        completed = 0
        while completed < job["steps"]:
            if time.monotonic() >= envelope["deadline"]:
                raise TimeoutError("ensemble deadline passed during this job")
            chunk = min(job["stride"], job["steps"] - completed)
            simulation.advance(chunk)
            completed += chunk
            record(simulation.observe())
        result["status"] = "success"
    except Exception as error:
        result["status"] = "timeout" if isinstance(error, TimeoutError) else "failed"
        result["error"] = {"type": type(error).__name__, "message": str(error)[:4096],
                           "traceback": traceback.format_exc(limit=6)[-8192:]}
    result["run_seconds"] = time.perf_counter() - started
    result["last_observed_step"] = result["summary"]["steps"][-1] if result["summary"]["steps"] else None
    return result


def _terminate_pool(pool: ProcessPoolExecutor) -> None:
    if hasattr(pool, "terminate_workers"):
        pool.terminate_workers()  # Public API on Python 3.14+.
        return
    # Python 3.13 has no public termination method. Bound the wait and explicitly
    # kill unfinished owned processes, then cancel queued jobs without resubmission.
    processes = list(getattr(pool, "_processes", {}).values())
    for process in processes:
        if process.is_alive():
            process.terminate()
    for process in processes:
        process.join(timeout=0.5)
        if process.is_alive():
            process.kill()
            process.join(timeout=0.5)
    pool.shutdown(wait=False, cancel_futures=True)


@dataclass
class EnsembleReport:
    metadata: dict[str, Any]
    results: list[dict[str, Any]]

    def to_dict(self) -> dict[str, Any]:
        return {"schema": "bdmc.ensemble.v1", "metadata": self.metadata, "results": self.results}

    def save(self, directory: str | Path, *, budget) -> Path:
        if budget is None:
            raise ValueError("saving an ensemble requires an explicit fixed DataBudget")
        directory = Path(directory)
        estimate = self.metadata["estimates"]["estimated_output_bytes"]
        budget.require_output(directory, planned_bytes=estimate)
        with budget.writer(directory, planned_bytes=estimate):
            return self._save_locked(directory, budget, estimate)

    def _save_locked(self, directory: Path, budget, estimate: int) -> Path:
        before = budget.require_output(directory, planned_bytes=estimate)
        if directory.exists():
            raise FileExistsError(directory)
        directory.mkdir(parents=True, exist_ok=False)
        from .storage import allocated_bytes

        # Recheck the whole conservative remaining estimate immediately before
        # the sole central report write; workers never create output files.
        budget.require_output(directory, planned_bytes=estimate)
        self.metadata["storage"] = {"before": before, "estimated_output_bytes": estimate,
                                    "actual_report_file_bytes": 0, "actual_output_allocated_bytes": 0}
        path = directory / "ensemble.json"
        for _ in range(5):
            # Stabilise the self-reported byte count (its decimal digits are
            # themselves part of the JSON). Allocation may grow by one block
            # when the after-check is first added, so remeasure between writes.
            for _ in range(5):
                text = json.dumps(self.to_dict(), indent=2, allow_nan=False) + "\n"
                size = len(text.encode())
                if self.metadata["storage"]["actual_report_file_bytes"] == size:
                    break
                self.metadata["storage"]["actual_report_file_bytes"] = size
            if size > estimate:
                raise ValueError("actual compact report exceeds its preflight estimate")
            budget.require_output(path, planned_bytes=estimate)
            path.write_text(text, encoding="utf-8")
            allocated = allocated_bytes(directory)
            if (self.metadata["storage"]["actual_output_allocated_bytes"] == allocated
                    and "after" in self.metadata["storage"]):
                break
            self.metadata["storage"]["actual_output_allocated_bytes"] = allocated
            self.metadata["storage"]["after"] = budget.check()
        else:
            raise RuntimeError("compact output accounting did not stabilise")
        return path


def execute_ensemble(jobs: Sequence[JobSpec], *, workers: int = 1, timeout_seconds: float = 60.,
                     output_directory: str | Path | None = None, budget=None,
                     limits: EnsembleLimits = EnsembleLimits()) -> EnsembleReport:
    """Execute independent jobs with a global queue-inclusive wall deadline.

    workers=1 is a serial spawn-process reference; workers>1 enables local process
    concurrency. Results are sorted by run_id, independent of completion order.
    Unfinished jobs remain explicit timeout records; no retry or reassignment.
    """
    jobs = tuple(jobs)
    if any(not isinstance(job, JobSpec) for job in jobs):
        raise ValueError("jobs must contain immutable JobSpec instances")
    estimates = _estimates(jobs, limits)
    if len({job.run_id for job in jobs}) != len(jobs):
        raise ValueError("run_id must be unique; no silent overwrite or deduplication")
    _integer(workers, "workers")
    if workers > min(limits.max_workers, os.cpu_count() or 1):
        raise ValueError("workers exceed advertised worker/CPU bounds")
    if (isinstance(timeout_seconds, bool) or not math.isfinite(timeout_seconds)
            or not 0 < timeout_seconds <= limits.max_timeout_seconds):
        raise ValueError("timeout_seconds must be finite, positive and within max_timeout_seconds")
    if output_directory is not None:
        if budget is None:
            raise ValueError("saved output requires an explicit fixed DataBudget")
        budget.require_output(output_directory, planned_bytes=estimates["estimated_output_bytes"])
        if Path(output_directory).exists():
            raise FileExistsError(output_directory)
    from . import _core
    core_hash = hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest()
    started = time.perf_counter()
    deadline = time.monotonic() + float(timeout_seconds)
    results = {}
    pool = ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
                               initializer=_worker_initialise)
    pending = {}
    terminated = False

    def receive(future, job):
        try:
            return future.result()
        except Exception as error:
            result = _base_result(job.to_dict(), core_hash)
            result.update(status="failed", error={"type": type(error).__name__, "message": str(error)[:4096]})
            return result

    try:
        for job in jobs:
            payload = json.dumps({"job": job.to_dict(), "expected_hash": core_hash, "deadline": deadline},
                                 allow_nan=False)
            pending[pool.submit(_execute_payload, payload)] = job
        while pending:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            finished, _ = wait(pending, timeout=remaining, return_when=FIRST_COMPLETED)
            if not finished:
                break
            for future in finished:
                job = pending.pop(future)
                results[job.run_id] = receive(future, job)
        if pending:
            for future, job in pending.items():
                if future.done():
                    results[job.run_id] = receive(future, job)
                    continue
                future.cancel()
                result = _base_result(job.to_dict(), core_hash)
                result.update(status="timeout", error={"type": "EnsembleDeadline",
                              "message": "global queue-inclusive deadline; unfinished worker state is not a checkpoint"})
                results[job.run_id] = result
            _terminate_pool(pool)
            terminated = True
    except Exception as error:
        # Dispatch/infrastructure failures do not discard completed jobs or
        # silently resubmit the remaining work to a different seed/worker.
        for future, job in pending.items():
            if future.done():
                results[job.run_id] = receive(future, job)
        for job in jobs:
            if job.run_id not in results:
                result = _base_result(job.to_dict(), core_hash)
                result.update(status="failed", error={"type": type(error).__name__, "message": str(error)[:4096]})
                results[job.run_id] = result
        _terminate_pool(pool)
        terminated = True
    except BaseException:
        _terminate_pool(pool)
        terminated = True
        raise
    finally:
        if not terminated:
            pool.shutdown(wait=True, cancel_futures=True)
    ordered = [results[key] for key in sorted(results)]
    # Retain every failed/timeout identity instead of silently drawing another
    # seed. A successful-only scientific dataset requires an explicit exclusion
    # decision; execution completion cannot certify unbiased physical sampling.
    metadata = {"backend": "local_cpu_spawn_processes", "workers": workers,
                "native_threads_per_worker": 1, "timeout_seconds": float(timeout_seconds),
                "timeout_policy": "one global wall deadline including queue/startup; terminate unfinished workers; no retry",
                "elapsed_seconds": time.perf_counter() - started, "core_binary_sha256": core_hash,
                "estimates": estimates, "limits": dict(limits.__dict__),
                "status_counts": {status: sum(result["status"] == status for result in ordered)
                                  for status in ("success", "failed", "timeout")},
                "result_order": "lexicographic run_id; independent of worker completion order",
                "identity_policy": "condition/lineage supplied by caller; lineage label is not restart state",
                "scientific_scope": "independent Brownian/dumbbell CPU reference runs; no GPU/MPI or modern network MC implementation"}
    report = EnsembleReport(metadata, ordered)
    if output_directory is not None:
        report.save(output_directory, budget=budget)
    return report
