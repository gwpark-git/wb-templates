"""Manifest, bounded launcher, MPI/reference equivalence and ML provenance."""
from dataclasses import asdict, replace
import hashlib
import json
import os
from pathlib import Path
import struct
import sys

import numpy as np
import pytest

from bdmc.mpi import (MPIExecutionError, MPIJob, MANIFEST_SCHEMA,
                      mpi_dataset, run_mpi_ensemble, serialize_jobs)
from bdmc.storage import DataBudget


@pytest.fixture
def budget(tmp_path):
    baseline = tmp_path/"baseline.json"
    root = tmp_path/"runs"
    root.mkdir()
    baseline.write_text(json.dumps({"schema": "bdmc.additional-data-budget.v1",
                                   "created_at": "unit-test-only", "baseline_allocated_bytes": {str(root): 0},
                                   "additional_limit_bytes": 10_000_000_000, "soft_limit_bytes": 9_000_000_000,
                                   "measurement": "scoped unit-test allocation"}))
    return DataBudget(baseline), root


def job(**changes):
    return MPIJob(**{"run_id": "run-0", "condition_id": "condition-0", "model": "hookean",
                     "seed": 713, "parameter": 1., **changes})


@pytest.mark.parametrize("changes", [{"seed": True}, {"seed": -1}, {"seed": 2**64},
                                    {"run_id": "bad\tname"}, {"condition_id": ""}, {"model": "network"},
                                    {"steps": 0}, {"count": 4097}, {"dt": float("nan")},
                                    {"diffusion": -1}, {"parameter": 0}, {"q0": 0},
                                    {"model": "repulsive", "shear_rate": .1},
                                    {"model": "repulsive", "box_length": 1.}, {"dt": 1.}])
def test_invalid_job_is_rejected(changes):
    with pytest.raises(ValueError):
        job(**changes)


def test_explicit_manifest_is_stable_and_bounded():
    jobs = [job(seed=0), job(run_id="run-1", model="repulsive", seed=2**64-1, parameter=0.)]
    text = serialize_jobs(jobs)
    assert text.startswith(MANIFEST_SCHEMA+"\nrun_id\tcondition_id\tlineage_id\tmodel\tseed\t")
    assert text == serialize_jobs(jobs)
    assert text.splitlines()[2].split("\t")[:5] == ["run-0", "condition-0", "run-0", "hookean", "0"]
    assert str(2**64-1) in text
    with pytest.raises(ValueError, match="unique"):
        serialize_jobs([job(), job()])
    with pytest.raises(ValueError, match="seeds must be unique"):
        serialize_jobs([job(), job(run_id="different-run", condition_id="different-condition", lineage_id="different-lineage")])
    with pytest.raises(ValueError, match="256"):
        serialize_jobs([])
    with pytest.raises(ValueError, match="update-work"):
        serialize_jobs([job(count=4096, steps=100000)])


def test_launcher_requires_fresh_registered_output_and_existing_budget(budget, monkeypatch):
    data_budget, root = budget
    monkeypatch.setattr("bdmc.mpi.shutil.which", lambda value: "/bin/false")
    with pytest.raises(TypeError, match="DataBudget"):
        run_mpi_ensemble([job()], executable="/bin/false", output=root/"a", budget=None)
    with pytest.raises(ValueError, match="outside"):
        run_mpi_ensemble([job()], executable="/bin/false", output=root.parent/"outside", budget=data_budget)
    with pytest.raises(FileExistsError):
        run_mpi_ensemble([job()], executable="/bin/false", output=root, budget=data_budget)
    for changes in ({"ranks": 0}, {"ranks": 64, "threads": 5}, {"timeout": 0},
                    {"backend": "pretend_gpu"}, {"launcher_args": "--oversubscribe"}):
        with pytest.raises(ValueError):
            run_mpi_ensemble([job()], executable="/bin/false", output=root/"unused", budget=data_budget, **changes)


def _successful_report(jobs):
    results = []
    for item in jobs:
        results.append({"run_id": item.run_id, "condition_id": item.condition_id, "lineage_id": item.lineage_id,
                        "model": item.model, "seed": item.seed, "status": "success",
                        "spec": {key: value for key, value in asdict(item).items() if key in
                                 ("parameter", "count", "steps", "dt", "diffusion", "mobility", "box_length", "sigma", "q0", "shear_rate", "warmup")},
                        "step": item.steps+item.warmup, "time": (item.steps+item.warmup)*item.dt,
                        "observation_lag": item.steps*item.dt,
                        "summary": {"connector_second_moment": .1, "polymer_shear_stress_per_dumbbell": -.02},
                        "initial_state_fnv1a64": "0000000000000000", "origin_state_fnv1a64": "0000000000000001",
                        "final_state_fnv1a64": "0000000000000002"})
    return {"schema": "bdmc.mpi-ensemble.v1", "status": "success", "jobs": [asdict(item) for item in jobs],
            "manifest_sha256": hashlib.sha256(serialize_jobs(jobs).encode()).hexdigest(),
            "executable_sha256": "unit-test", "native": {"schema": "bdmc.mpi-ensemble.v1", "status": "success",
                "execution": {"backend": "serial", "threads_per_rank": 1, "ranks": 1}, "results": results}}


def test_dataset_preserves_condition_lineage_and_manifest_without_seed_features():
    jobs = [job(), job(run_id="run-1", seed=714, lineage_id="continuation-1"),
            job(run_id="run-2", seed=715, condition_id="condition-1", parameter=2.)]
    report = _successful_report(jobs)
    data = mpi_dataset(report, model="hookean", target_names=("connector_second_moment",))
    np.testing.assert_array_equal(data.groups, [item.run_id for item in jobs])
    np.testing.assert_array_equal(data.group_keys["condition"], [item.condition_id for item in jobs])
    np.testing.assert_array_equal(data.group_keys["continuation"], [item.lineage_id for item in jobs])
    assert "seed" not in data.feature_names
    assert data.metadata["runs"][0]["seed"] == 713
    assert data.metadata["execution"]["backend"] == "serial"
    assert data.metadata["execution"]["arithmetic_precision"] == "float64"
    report["native"]["execution"]["ranks"] = 3
    assert mpi_dataset(report, model="hookean", target_names=("connector_second_moment",)).fingerprint == data.fingerprint
    report["native"]["execution"]["backend"] = "metal"
    metal = mpi_dataset(report, model="hookean", target_names=("connector_second_moment",))
    assert metal.metadata["execution"]["arithmetic_precision"] == "float32"
    assert metal.fingerprint != data.fingerprint
    for features in (("seed",), ("run_id",), ("count", "count")):
        with pytest.raises(ValueError, match="features"):
            mpi_dataset(report, model="hookean", feature_names=features, target_names=("connector_second_moment",))
    report["native"]["results"][0]["seed"] += 1
    with pytest.raises(ValueError, match="identity"):
        mpi_dataset(report, model="hookean", target_names=("connector_second_moment",))


@pytest.mark.parametrize("mutation", ["parameter", "time", "manifest"])
def test_dataset_rejects_changed_physical_specification_or_manifest(mutation):
    report = _successful_report([job()])
    if mutation == "parameter":
        report["native"]["results"][0]["spec"]["parameter"] += .1
    elif mutation == "time":
        report["native"]["results"][0]["time"] += .1
    else:
        report["jobs"][0]["parameter"] += .1
    with pytest.raises(ValueError, match="identity|SHA256"):
        mpi_dataset(report, model="hookean", target_names=("connector_second_moment",))


def test_launcher_retains_failures_without_retry_or_reseed(budget, monkeypatch):
    data_budget, root = budget
    calls = []

    def failed(command, **kwargs):
        calls.append(command)
        return 3, b'{"schema":"bdmc.mpi-ensemble.v1","status":"failed","error":"FENE"}', b"diagnostic", None

    monkeypatch.setattr("bdmc.mpi.shutil.which", lambda value: "/bin/false")
    monkeypatch.setattr("bdmc.mpi._bounded_process", failed)
    with pytest.raises(MPIExecutionError) as caught:
        run_mpi_ensemble([job()], executable="/bin/false", output=root/"failed", budget=data_budget)
    assert len(calls) == 1
    assert caught.value.report["status"] == "failed"
    assert caught.value.report["jobs"][0]["seed"] == 713
    saved = json.loads((root/"failed"/"report.json").read_text())
    assert saved["native"]["error"] == "FENE"
    assert (root/"failed"/"stderr.log").read_text() == "diagnostic"


def test_timeout_terminates_real_launcher_group_and_saves_failure(budget, monkeypatch):
    data_budget, root = budget
    monkeypatch.setattr("bdmc.mpi.shutil.which", lambda value: "/bin/sh")
    with pytest.raises(MPIExecutionError, match="timeout") as caught:
        run_mpi_ensemble([job()], executable="/bin/true", output=root/"timed-out", budget=data_budget,
                         launcher_args=("-c", "sleep 30 & wait", "timeout-test"), timeout=.1)
    assert caught.value.report["elapsed_seconds"] < 5
    assert caught.value.report["returncode"] != 0
    assert json.loads((root/"timed-out"/"report.json").read_text())["status"] == "failed"


def test_process_output_bound_terminates_actual_process(monkeypatch, tmp_path):
    from bdmc import mpi
    monkeypatch.setattr(mpi, "MAX_STDOUT_BYTES", 1024)
    code, stdout, _, problem = mpi._bounded_process(
        [sys.executable, "-c", "import sys,time; sys.stdout.buffer.write(b'x'*8192); sys.stdout.flush(); time.sleep(30)"],
        directory=tmp_path, environment=os.environ.copy(), timeout=5.)
    assert code != 0
    assert len(stdout) == 1024
    assert "output exceeded" in problem


def _fingerprint(values):
    value = 14695981039346656037
    for number in np.asarray(values).ravel():
        for byte in struct.pack("<d", float(number)):
            value = ((value ^ byte)*1099511628211) % 2**64
    return f"{value:016x}"


@pytest.mark.skipif(not os.environ.get("BDMC_MPI_EXECUTABLE"), reason="requires explicitly built native MPI executable")
def test_actual_mpi_matches_serial_reference_and_feeds_grouped_learning(budget):
    from bdmc.reference import DumbbellSimulation, DumbbellSpec, RepulsiveSimulation, RepulsiveSpec
    from bdmc.ml import grouped_split, load_artifact, train

    data_budget, root = budget
    jobs = [job(run_id=f"run-{index}", condition_id=f"condition-{index//2}",
                model=model, seed=713+index, parameter=.8+.1*index, warmup=2,
                box_length=2. if model == "repulsive" else 4.,
                shear_rate=0 if model == "repulsive" else (-.75 if index % 2 else .75))
            for index, model in enumerate(("repulsive", "hookean", "fene", "hookean", "fene", "hookean"))]
    reports = []
    for ranks in (1, 2, 3):
        report = run_mpi_ensemble(jobs, executable=os.environ["BDMC_MPI_EXECUTABLE"], output=root/f"ranks-{ranks}",
                                  budget=data_budget, ranks=ranks, timeout=30,
                                  launcher=os.environ.get("BDMC_MPI_LAUNCHER", "mpiexec"),
                                  launcher_args=tuple(json.loads(os.environ.get("BDMC_MPI_LAUNCHER_ARGS", "[]"))))
        reports.append(report)
    assert reports[0]["native"]["results"] == reports[1]["native"]["results"] == reports[2]["native"]["results"]
    for item, result in zip(jobs, reports[0]["native"]["results"]):
        if item.model == "repulsive":
            sim = RepulsiveSimulation(RepulsiveSpec(n_particles=item.count, box_length=item.box_length,
                dt=item.dt, diffusion=item.diffusion, mobility=item.mobility, epsilon=item.parameter,
                sigma=item.sigma, seed=item.seed))
            key = "unwrapped_positions"
        else:
            sim = DumbbellSimulation(DumbbellSpec(n_dumbbells=item.count, dt=item.dt,
                diffusion=item.diffusion, mobility=item.mobility, spring_constant=item.parameter,
                spring_law=item.model, q0=item.q0, protocol="simple_shear", shear_rate=item.shear_rate, seed=item.seed))
            key = "connectors"
        assert _fingerprint(sim.observe()[key]) == result["initial_state_fnv1a64"]
        sim.advance(item.warmup)
        assert _fingerprint(sim.observe()[key]) == result["origin_state_fnv1a64"]
        sim.advance(item.steps)
        assert _fingerprint(sim.observe()[key]) == result["final_state_fnv1a64"]
    data = mpi_dataset(reports[0], model="hookean", target_names=("connector_second_moment",))
    split = grouped_split(data, group_by=("condition",))
    split.validate(data)
    assert len(data.X) == 3
    artifact = train(data, method="ridge", split=split, seed=19)
    with data_budget.writer(root/"trained", planned_bytes=100_000):
        artifact.save(root/"trained")
    np.testing.assert_allclose(load_artifact(root/"trained").predict(data.X), artifact.predict(data.X), rtol=0, atol=0)
