"""Actual spawn-process contracts and bounded compact output, no large runs."""
from dataclasses import FrozenInstanceError, replace
from concurrent.futures import Future
import json
import multiprocessing
import time

import numpy as np
import pytest

from bdmc import Simulation, SimulationSpec
from bdmc.ensemble import (EnsembleLimits, JobSpec, execute_ensemble,
                            _base_result, _execute_payload, _worker_initialise)
from bdmc.storage import DataBudget, DataBudgetExceeded, allocated_bytes, create_data_baseline


def brownian_job(run_id="run-b", seed=17, **changes):
    arguments = dict(run_id=run_id, spec=SimulationSpec(n_particles=8, dimension=2, box_lengths=[10., 10.],
                                                       dt=.01, diffusivity=.5, seed=seed),
                     steps=7, stride=3, observables=("msd", "potential_energy", "mean_unwrapped_position"),
                     condition_id="condition-b", lineage_id=f"lineage-{run_id}")
    arguments.update(changes)
    return JobSpec(**arguments)


def dumbbell_job():
    return JobSpec("run-a", SimulationSpec(model="dumbbell", n_particles=2, dimension=1,
                    box_lengths=[20.], dt=.1, diffusivity=0., seed=8), steps=2, stride=1,
                    observables=("connector_second_moment", "potential_energy", "msd"),
                    condition_id="condition-a", lineage_id="lineage-a",
                    initial_positions=((0.,), (2.,)), particle_ids=(19, 4))


def test_job_snapshot_is_deeply_immutable_and_json_roundtrippable():
    spec = SimulationSpec(n_particles=4, dimension=1, box_lengths=[9.], seed=1)
    positions = [[0.], [1.], [2.], [3.]]
    job = brownian_job(spec=spec, initial_positions=positions)
    spec.seed = 999
    spec.box_lengths = [1.]
    positions[0][0] = 100
    assert job.spec["seed"] == 1 and job.spec["box_lengths"] == (9.,)
    assert job.initial_positions[0] == (0.,)
    with pytest.raises(TypeError):
        job.spec["seed"] = 2
    with pytest.raises(FrozenInstanceError):
        job.run_id = "changed"
    restored = JobSpec(**json.loads(json.dumps(job.to_dict())))
    assert restored.to_dict() == job.to_dict()


def test_actual_serial_and_two_worker_results_match_each_seed_and_identity():
    jobs = [brownian_job(), dumbbell_job(), brownian_job("run-c", seed=18)]
    serial = execute_ensemble(jobs, workers=1, timeout_seconds=20)
    parallel = execute_ensemble(tuple(reversed(jobs)), workers=2, timeout_seconds=20)
    assert [row["run_id"] for row in parallel.results] == ["run-a", "run-b", "run-c"]
    for a, b in zip(serial.results, parallel.results):
        assert a["status"] == b["status"] == "success"
        for key in ("run_id", "condition_id", "lineage_id", "seed", "spec", "core_binary_sha256", "summary"):
            assert a[key] == b[key]
        assert a["core_hash_verified"] and b["core_hash_verified"]
        assert a["native_thread_limit"] == b["native_thread_limit"] == 1
        if b["native_thread_pools"]:
            assert all(pool["num_threads"] == 1 for pool in b["native_thread_pools"])
    assert parallel.results[1]["summary"]["values"]["msd"] != parallel.results[2]["summary"]["values"]["msd"]
    assert parallel.metadata["status_counts"] == {"success": 3, "failed": 0, "timeout": 0}


def test_dumbbell_compact_observables_match_known_deterministic_answers():
    result = execute_ensemble([dumbbell_job()], timeout_seconds=20).results[0]
    assert result["status"] == "success"
    assert result["summary"]["steps"] == [0, 1, 2]
    np.testing.assert_allclose(result["summary"]["values"]["connector_second_moment"],
                               [4., 2.56, 1.6384], rtol=1e-13)
    np.testing.assert_allclose(result["summary"]["values"]["potential_energy"],
                               [2., 1.28, .8192], rtol=1e-13)
    np.testing.assert_allclose(result["summary"]["values"]["msd"], [0., .04, .1296], atol=1e-13)


def test_unwrapped_displacements_and_final_partial_stride_match_direct_core():
    spec = SimulationSpec(n_particles=8, dimension=2, box_lengths=[.1, .1], dt=.01,
                          diffusivity=1., seed=37)
    job = brownian_job(spec=spec)
    direct = Simulation(spec).sample(7, stride=3, include_initial=True)
    expected = np.square(direct["unwrapped_positions"]).sum(axis=2).mean(axis=1)
    result = execute_ensemble([job], timeout_seconds=20).results[0]
    assert result["summary"]["steps"] == [0, 3, 6, 7]
    np.testing.assert_array_equal(result["summary"]["times"], direct["times"])
    np.testing.assert_array_equal(result["summary"]["values"]["msd"], expected)
    assert expected[-1] > 2 * .1 ** 2  # A wrapped-coordinate MSD would be bounded here.
    assert "positions" not in result["summary"] and "forces" not in result["summary"]


def test_worker_failure_is_retained_beside_success_with_no_retry_or_partial_frame():
    failure = JobSpec("bad", SimulationSpec(n_particles=2, dimension=1), steps=1, stride=1,
                      observables=("mean_unwrapped_position",), condition_id="extreme",
                      lineage_id="lineage-bad", initial_positions=((1e308,), (1e308,)))
    report = execute_ensemble([failure, brownian_job("good")], workers=2, timeout_seconds=20)
    failed, passed = report.results
    assert failed["run_id"] == "bad" and failed["status"] == "failed"
    assert failed["error"]["type"] == "FloatingPointError"
    assert failed["summary"]["steps"] == []
    assert failed["summary"]["values"]["mean_unwrapped_position"] == []
    assert passed["status"] == "success" and len(report.results) == 2
    assert report.metadata["status_counts"] == {"success": 1, "failed": 1, "timeout": 0}


def test_global_timeout_preserves_every_job_and_terminates_owned_workers():
    before = {process.pid for process in multiprocessing.active_children()}
    started = time.monotonic()
    report = execute_ensemble([brownian_job("early"), brownian_job("queued")], workers=1,
                              timeout_seconds=1e-6)
    assert time.monotonic() - started < 10
    assert [row["run_id"] for row in report.results] == ["early", "queued"]
    assert all(row["status"] == "timeout" for row in report.results)
    assert all(row["error"]["type"] == "EnsembleDeadline" for row in report.results)
    assert not ({process.pid for process in multiprocessing.active_children()} - before)


def test_binary_change_is_a_retained_worker_failure_not_a_silent_backend_switch():
    payload = json.dumps({"job": brownian_job().to_dict(), "expected_hash": "0" * 64,
                          "deadline": time.monotonic() + 20})
    result = _execute_payload(payload)
    assert result["status"] == "failed" and not result["core_hash_verified"]
    assert "core binary changed" in result["error"]["message"]


def test_saved_compact_report_requires_budget_preflight_and_records_actual_growth(tmp_path):
    root = tmp_path / "runs"
    root.mkdir()
    baseline = tmp_path / "budget.json"
    create_data_baseline([root], baseline)
    budget = DataBudget(baseline)
    output = root / "ensemble"
    report = execute_ensemble([dumbbell_job()], timeout_seconds=20, output_directory=output, budget=budget)
    saved = output / "ensemble.json"
    record = json.loads(saved.read_text())
    assert set(output.iterdir()) == {saved}
    assert record["metadata"]["storage"]["actual_report_file_bytes"] == saved.stat().st_size
    assert record["metadata"]["storage"]["actual_output_allocated_bytes"] == allocated_bytes(output)
    assert record == report.to_dict()
    assert saved.stat().st_size < report.metadata["estimates"]["estimated_output_bytes"]
    assert saved.stat().st_size < 2_000_000
    assert record["results"][0]["status"] == "success"
    with pytest.raises(FileExistsError):
        report.save(output, budget=budget)
    with pytest.raises(ValueError, match="explicit fixed DataBudget"):
        report.save(root / "unbudgeted", budget=None)


def test_output_scope_or_budget_rejection_happens_before_worker_launch(tmp_path, monkeypatch):
    root = tmp_path / "runs"
    root.mkdir()
    baseline = tmp_path / "budget.json"
    create_data_baseline([root], baseline)
    budget = DataBudget(baseline)
    from bdmc import ensemble

    def unexpected(*args, **kwargs):
        raise AssertionError("workers must not launch before output preflight succeeds")

    monkeypatch.setattr(ensemble, "ProcessPoolExecutor", unexpected)
    with pytest.raises(ValueError, match="explicit fixed DataBudget"):
        execute_ensemble([brownian_job()], output_directory=root / "missing")
    with pytest.raises(ValueError, match="outside"):
        execute_ensemble([brownian_job()], output_directory=tmp_path / "outside", budget=budget)
    record = json.loads(baseline.read_text())
    record["soft_limit_bytes"] = 1
    baseline.write_text(json.dumps(record))
    with pytest.raises(DataBudgetExceeded):
        execute_ensemble([brownian_job()], output_directory=root / "over-budget", budget=budget)
    assert not (root / "over-budget").exists()


@pytest.mark.parametrize("change", [dict(steps=-1), dict(stride=0), dict(run_id=""),
                                    dict(condition_id=""), dict(lineage_id=""),
                                    dict(observables=("unknown",)), dict(observables=("msd", "msd")),
                                    dict(observables=("connector_second_moment",)),
                                    dict(initial_positions=((0.,),)), dict(particle_ids=(1, 1))])
def test_invalid_job_contract_is_rejected(change):
    with pytest.raises((ValueError, TypeError)):
        brownian_job(**change)


@pytest.mark.parametrize("change", [dict(workers=0), dict(workers=999), dict(timeout_seconds=0),
                                    dict(timeout_seconds=float("nan")),
                                    dict(timeout_seconds=301),
                                    dict(limits=EnsembleLimits(max_frames_per_job=1)),
                                    dict(limits=EnsembleLimits(max_particles=1)),
                                    dict(limits=EnsembleLimits(max_particle_steps=1)),
                                    dict(limits=EnsembleLimits(max_estimated_output_bytes=1))])
def test_excessive_or_invalid_execution_is_rejected(change):
    with pytest.raises(ValueError):
        execute_ensemble([brownian_job()], **change)


def test_duplicate_ids_and_empty_jobs_are_rejected():
    with pytest.raises(ValueError, match="unique"):
        execute_ensemble([brownian_job(), brownian_job()])
    with pytest.raises(ValueError, match="1..max_jobs"):
        execute_ensemble([])


def test_dispatch_failure_retains_completed_and_unscheduled_jobs_without_retry(monkeypatch):
    from bdmc import ensemble
    calls = []

    class BrokenDispatcher:
        _processes = {}

        def __init__(self, **kwargs):
            pass

        def submit(self, function, payload):
            envelope = json.loads(payload)
            calls.append(envelope["job"]["run_id"])
            if len(calls) == 2:
                raise RuntimeError("dispatch infrastructure failed")
            result = _base_result(envelope["job"], envelope["expected_hash"])
            result["status"] = "success"
            future = Future()
            future.set_result(result)
            return future

        def shutdown(self, **kwargs):
            assert kwargs["wait"] is False

    monkeypatch.setattr(ensemble, "ProcessPoolExecutor", BrokenDispatcher)
    report = execute_ensemble([brownian_job("a"), brownian_job("b"), brownian_job("c")])
    assert calls == ["a", "b"]  # No retry, and c was never submitted.
    assert [row["status"] for row in report.results] == ["success", "failed", "failed"]
    assert report.results[1]["error"]["message"] == "dispatch infrastructure failed"
    assert report.results[2]["error"]["message"] == "dispatch infrastructure failed"
