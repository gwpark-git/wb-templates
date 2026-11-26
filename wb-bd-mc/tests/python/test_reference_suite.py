"""Negative controls for the frozen-reference checker, without simulations."""

import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "reference_suite.py"
SPEC = importlib.util.spec_from_file_location("reference_suite", SCRIPT)
suite = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(suite)


def write_table(path, rows):
    path.write_text("".join("\t".join(str(value) for value in row) + "\n" for row in rows))


@pytest.fixture
def network(tmp_path):
    """Two particles, two conserved chains: loops become two bridges."""
    run = tmp_path / "network"
    run.mkdir()
    parameters = {**suite.COMMON, "Method": "NAPLE_ASSOCIATION", "Step": "NO_EQUILIBRATION",
                  "Np": "2", "N_chains_per_particle": "1", "Nt": "3",
                  "N_skip_file": "1", "N_skip_ener": "1", "cell_list": "FALSE"}
    positions = [[0., 1., 2., 3., 0., 0., 0., 4., 5., 1., 0., 0., 0.],
                 [.001, 1.01, 2., 3., 0., 0., 0., 4., 5.01, 1., 0., 0., 0.]]
    energy = [[0.] * 31 for _ in range(2)]
    energy[0][5] = .02
    energy[1][0], energy[1][4], energy[1][5], energy[1][30] = .001, 2., .03, 1.
    write_table(run / "case.traj", positions)
    write_table(run / "case.ener", energy)
    write_table(run / "case.chain", [[0, 1, 0, 1], [0, 1, 1, 0]])
    write_table(run / "case.hash", [[0], [1], [0, 1], [1, 0]])
    write_table(run / "case.weight", [[2], [2], [0, 2], [0, 2]])
    (run / "stdout.log").write_text("START SIMULATION\nTotal simulation time = 0.1\n")
    return run, parameters


def candidate(network, tmp_path):
    output = tmp_path / "candidate"
    shutil.copytree(network[0], output)
    return output


def mutate(path, row, column, value):
    rows = suite.table(path)
    rows[row][column] = value
    write_table(path, rows)


def test_checker_accepts_independently_constructed_conserved_topology(network):
    summary = suite.validate_run(*network)
    assert summary["max_bridges"] == 2
    assert summary["chain_end_count"] == 4
    assert summary["distinct_topology_frames"] == 2


def test_mutated_scientific_field_fails_exact_and_tolerant_comparison(network, tmp_path):
    output = candidate(network, tmp_path)
    mutate(output / "case.ener", 1, 9, .125)
    for exact in (False, True):
        with pytest.raises(AssertionError, match="ener: scientific difference.*column 9"):
            suite.compare_runs(network[0], output, exact=exact)


def test_only_energy_timing_column_is_excluded(network, tmp_path):
    output = candidate(network, tmp_path)
    mutate(output / "case.ener", 1, 5, 9000.)
    suite.validate_run(output, network[1])
    suite.compare_runs(network[0], output, exact=True)
    mutate(output / "case.ener", 1, 6, .001)
    with pytest.raises(AssertionError, match="column 6"):
        suite.compare_runs(network[0], output)


@pytest.mark.parametrize("suffix,column", [("traj", 1), ("ener", 9), ("ener", 5)])
def test_nonfinite_rejected_even_in_ignored_timing(network, tmp_path, suffix, column):
    output = candidate(network, tmp_path)
    mutate(output / f"case.{suffix}", 1, column, float("nan"))
    with pytest.raises(AssertionError, match="nonfinite"):
        suite.validate_run(output, network[1])
    with pytest.raises(AssertionError, match="nonfinite"):
        suite.compare_runs(network[0], output)


@pytest.mark.parametrize("suffix", suite.SUFFIXES)
def test_missing_table_rejected(network, tmp_path, suffix):
    output = candidate(network, tmp_path)
    (output / f"case.{suffix}").unlink()
    with pytest.raises(AssertionError, match="missing table"):
        suite.validate_run(output, network[1])


@pytest.mark.parametrize("mutation", ["short_row", "extra_row", "blank_row", "nonnumeric", "timestamp"])
def test_table_schema_and_frame_controls(network, tmp_path, mutation):
    output = candidate(network, tmp_path)
    path = output / "case.traj"
    rows = suite.table(path)
    if mutation == "short_row":
        rows[1].pop()
        write_table(path, rows)
    elif mutation == "extra_row":
        write_table(path, rows + rows[-1:])
    elif mutation == "blank_row":
        path.write_text(path.read_text() + "\n")
    elif mutation == "nonnumeric":
        path.write_text(path.read_text().replace("1.01", "corrupt"))
    else:
        mutate(path, 1, 0, .002)
    with pytest.raises(AssertionError):
        suite.validate_run(output, network[1])


def test_reconstructed_chain_topology_catches_corrupt_adjacency(network, tmp_path):
    output = candidate(network, tmp_path)
    mutate(output / "case.weight", 3, 1, 1)
    with pytest.raises(AssertionError, match="adjacency/labelled-chain"):
        suite.validate_run(output, network[1])


def test_topology_integers_compare_exactly(network, tmp_path):
    output = candidate(network, tmp_path)
    mutate(output / "case.chain", 1, 0, 0.00000001)
    with pytest.raises(AssertionError, match="nonintegral topology"):
        suite.validate_run(output, network[1])
    with pytest.raises(AssertionError, match="chain: scientific difference"):
        suite.compare_runs(network[0], output)


@pytest.mark.parametrize("diagnostic", ["Bad condition missing", "ERR: unsupported", "ERROR: failure",
                                       "runtime error: out of bounds", "AddressSanitizer: overflow"])
def test_runtime_diagnostics_rejected(network, tmp_path, diagnostic):
    output = candidate(network, tmp_path)
    with (output / "stdout.log").open("a") as stream:
        stream.write(diagnostic + "\n")
    with pytest.raises(AssertionError, match="input/runtime diagnostic"):
        suite.validate_run(output, network[1])


def checksum_fixture(tmp_path, manifest):
    dataset = tmp_path / "frozen"
    dataset.mkdir()
    (dataset / "manifest.json").write_text(json.dumps(manifest))
    suite.write_json(dataset / "CHECKSUMS.json", {"algorithm": "sha256", "files": {
        "manifest.json": suite.sha256(dataset / "manifest.json")}})
    return dataset


def test_tampered_fixture_rejected_before_any_process_or_output(tmp_path):
    dataset = checksum_fixture(tmp_path, {"schema": suite.SCHEMA})
    (dataset / "manifest.json").write_text('{"schema":"tampered"}')
    run_dir = tmp_path / "never-created"
    with pytest.raises(AssertionError, match="checksum mismatch"):
        suite.verify(Path("/nonexistent-executable"), dataset, run_dir, None)
    assert not run_dir.exists()


def test_rehashed_wrong_manifest_schema_is_still_rejected(tmp_path):
    dataset = checksum_fixture(tmp_path, {"schema": "unknown-schema"})
    with pytest.raises(AssertionError, match="manifest schema/scope"):
        suite.load_dataset(dataset)


def test_missing_or_unlisted_fixture_is_rejected(tmp_path):
    dataset = checksum_fixture(tmp_path, {"schema": suite.SCHEMA})
    (dataset / "unlisted.txt").write_text("not registered")
    with pytest.raises(AssertionError, match="inventory mismatch"):
        suite.load_dataset(dataset)


def test_create_refuses_to_replace_even_an_empty_existing_dataset(tmp_path):
    dataset = tmp_path / "existing"
    dataset.mkdir()
    with pytest.raises(FileExistsError, match="already exists"):
        suite.create(Path("/nonexistent-executable"), dataset, tmp_path / "runs", None)
    assert list(dataset.iterdir()) == []


def test_frozen_corpus_integrity_and_protocol_activity():
    dataset = SCRIPT.parent / "reference/v1"
    assert dataset.is_dir(), "shipped frozen reference corpus is missing"
    manifest = suite.load_dataset(dataset)
    assert len(manifest["cases"]) == 9
    assert len(manifest["protocol_checks"]) == 4


def test_brownian_clock_revisions_reject_each_others_semantics(tmp_path):
    source = SCRIPT.parent / "reference/v1/brownian"
    run = tmp_path / "brownian"
    shutil.copytree(source, run)
    parameters = suite.parse_input((run / "input.inp").read_text())
    suite.validate_run(run, parameters, "legacy-v1")
    with pytest.raises(AssertionError, match="energy time-field mismatch"):
        suite.validate_run(run, parameters, "repaired-v2")
    rows = suite.table(run / "case.ener")
    trajectory = suite.table(run / "case.traj")
    for row, frame in zip(rows, trajectory):
        row[0] = frame[0]
    write_table(run / "case.ener", rows)
    summary = suite.validate_run(run, parameters, "repaired-v2")
    assert summary["energy_time_field"] == "pre-step force/frame timestamp"
    with pytest.raises(AssertionError, match="energy time-field mismatch"):
        suite.validate_run(run, parameters, "legacy-v1")


def test_unknown_behavioral_revision_cannot_reinterpret_an_authenticated_corpus(tmp_path):
    source = SCRIPT.parent / "reference/v1"
    dataset = tmp_path / "frozen"
    shutil.copytree(source, dataset)
    manifest_path = dataset / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["behavioral_revision"] = "unknown-semantics"
    manifest_path.write_text(json.dumps(manifest))
    checksums_path = dataset / "CHECKSUMS.json"
    checksums = json.loads(checksums_path.read_text())
    checksums["files"]["manifest.json"] = suite.sha256(manifest_path)
    checksums_path.write_text(json.dumps(checksums))
    with pytest.raises(AssertionError, match="unsupported behavioral revision"):
        suite.load_dataset(dataset)


def test_helper_import_does_not_change_python_package_path():
    original_path = sys.path[:]
    spec = importlib.util.spec_from_file_location("reference_suite_path_control", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert sys.path == original_path


@pytest.mark.parametrize("parent_exits", [False, True])
def test_process_group_cleanup_blocks_delayed_descendant_writes(tmp_path, parent_exits):
    marker = tmp_path / "escaped-write.txt"
    ready = tmp_path / "descendant-ready.txt"
    script = tmp_path / "fork_writer.py"
    script.write_text(
        "import os, pathlib, signal, sys, time\n"
        "child = os.fork()\n"
        "if child == 0:\n"
        "    signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
        "    pathlib.Path(sys.argv[2]).write_text(str(os.getpid()))\n"
        "    time.sleep(0.5)\n"
        "    pathlib.Path(sys.argv[1]).write_text('descendant survived')\n"
        "    os._exit(0)\n"
        "while not pathlib.Path(sys.argv[2]).exists(): time.sleep(0.002)\n"
        "if sys.argv[3] == 'exit': sys.exit(0)\n"
        "time.sleep(30)\n")
    command = [sys.executable, str(script), str(marker), str(ready), "exit" if parent_exits else "wait"]
    previous_handler = signal.getsignal(signal.SIGTERM)
    with (tmp_path / "process.log").open("w") as log:
        if parent_exits:
            process = suite.run_process(command, cwd=tmp_path, env=os.environ.copy(), log=log, timeout=1)
            assert process.returncode == 0
        else:
            with pytest.raises(subprocess.TimeoutExpired):
                suite.run_process(command, cwd=tmp_path, env=os.environ.copy(), log=log, timeout=.15)
    assert signal.getsignal(signal.SIGTERM) == previous_handler
    assert ready.exists(), "fake executable must actually create a descendant before cleanup"
    time.sleep(.6)
    assert not marker.exists(), "ordinary descendant escaped process-group cleanup"


def test_managed_mode_rejects_a_wrapper_that_is_not_group_leader(tmp_path, monkeypatch):
    monkeypatch.setattr(suite.os, "getpgrp", lambda: os.getpid() + 1)
    with (tmp_path / "process.log").open("w") as log:
        with pytest.raises(AssertionError, match="must be process-group leader"):
            suite.run_process([sys.executable, "-c", "raise SystemExit(99)"], cwd=tmp_path,
                              env={**os.environ, "BDMC_REFERENCE_MANAGED_GROUP": "1"}, log=log, timeout=.1)


@pytest.mark.parametrize("mode", ["normal_parent", "helper_timeout", "outer_timeout"])
def test_managed_outer_group_contains_simulator_descendants(tmp_path, mode):
    marker = tmp_path / "escaped-write.txt"
    ready = tmp_path / "descendant-ready.txt"
    child = tmp_path / "fork_writer.py"
    child.write_text(
        "import os, pathlib, signal, sys, time\n"
        "pid = os.fork()\n"
        "if pid == 0:\n"
        "    signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
        "    pathlib.Path(sys.argv[2]).write_text(str(os.getpgrp()))\n"
        "    time.sleep(0.6)\n"
        "    pathlib.Path(sys.argv[1]).write_text('escaped managed group')\n"
        "    os._exit(0)\n"
        "while not pathlib.Path(sys.argv[2]).exists(): time.sleep(0.002)\n"
        "if sys.argv[3] == 'normal_parent': sys.exit(0)\n"
        "time.sleep(30)\n")
    driver = tmp_path / "managed_driver.py"
    driver.write_text(
        "import importlib.util, os, pathlib, subprocess, sys\n"
        f"spec=importlib.util.spec_from_file_location('reference', {str(SCRIPT)!r})\n"
        "suite=importlib.util.module_from_spec(spec); spec.loader.exec_module(suite)\n"
        "with pathlib.Path(sys.argv[5]).open('w') as log:\n"
        "    try:\n"
        "        suite.run_process([sys.executable,*sys.argv[1:5]],cwd=pathlib.Path.cwd(),\n"
        "                          env=os.environ.copy(),log=log,\n"
        "                          timeout=.15 if sys.argv[4]=='helper_timeout' else 30)\n"
        "    except subprocess.TimeoutExpired: sys.exit(3)\n")
    env = {**os.environ, "BDMC_REFERENCE_MANAGED_GROUP": "1"}
    process = subprocess.Popen([sys.executable, str(driver), str(child), str(marker), str(ready), mode,
                                str(tmp_path / "process.log")], cwd=tmp_path, env=env, start_new_session=True)
    try:
        if mode == "outer_timeout":
            with pytest.raises(subprocess.TimeoutExpired):
                process.wait(timeout=.25)
        else:
            assert process.wait(timeout=2) == (3 if mode == "helper_timeout" else 0)
    finally:
        # This is the outer unified runner's ownership boundary: every managed
        # simulator child must share this group even when the wrapper has exited.
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()
    assert ready.exists(), "fake simulator must create a real descendant"
    assert int(ready.read_text()) == process.pid, "simulator escaped the outer managed process group"
    time.sleep(.7)
    assert not marker.exists(), "outer cleanup failed to stop managed descendant"
