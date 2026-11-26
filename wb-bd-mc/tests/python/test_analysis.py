"""Coordination integration and restored explicit-column text tools."""
import os
from pathlib import Path
import subprocess

import numpy as np
import pytest

from bdmc.analysis import autocorrelation, coordination_number, pair_correlation


@pytest.mark.parametrize("normalization", ["finite_N", "bulk"])
def test_coordination_exact_counts_and_shell_integral(normalization):
    r = pair_correlation([[0, 0, 0], [1, 0, 0], [2.5, 0, 0]], [10] * 3,
                         [0, .5, 1.2, 2, 3], normalization=normalization)
    np.testing.assert_allclose(coordination_number(r), [0, 2/3, 4/3, 2])
    np.testing.assert_allclose(coordination_number(r, source="g"), coordination_number(r))
    assert coordination_number(r, 3) == 2
    assert coordination_number(r, 0) == 0
    with pytest.raises(ValueError, match="stored radial edge"):
        coordination_number(r, 1.1)
    with pytest.raises(ValueError, match="source"):
        coordination_number(r, source="unknown")
    with pytest.raises(ValueError):
        coordination_number({**r, "kind": "directional_3d"})


def test_missing_inner_shell_requires_explicit_truncation():
    r = pair_correlation([[0], [1]], [10], [.5, 1.5])
    with pytest.raises(ValueError, match="above zero"):
        coordination_number(r)
    np.testing.assert_array_equal(coordination_number(r, allow_truncated=True), [1])


@pytest.fixture
def executables():
    directory = os.environ.get("BDMC_ANALYSIS_TOOL_DIR")
    if not directory:
        pytest.skip("BDMC_ANALYSIS_TOOL_DIR not set; standalone tools not exercised")
    paths = tuple(Path(directory) / name for name in ["measure_stress_autocorrelation", "correlation_function_individual_chain"])
    for path in paths:
        assert path.is_file(), f"declared tool missing: {path}"
    return paths


def run_tool(executable, source, output, *arguments):
    return subprocess.run([str(executable), str(source), str(output), *arguments],
                          text=True, capture_output=True, timeout=30)


def test_stress_tool_named_and_index_columns_both_estimators(executables, tmp_path):
    stress, _ = executables
    source = tmp_path / "stress.dat"
    values = np.array([[10., 1, 5], [10.5, 2, 4], [11., 3, 1], [11.5, 4, 0], [12., 7, -1]])
    np.savetxt(source, values, header="time sigma_xy sigma_xz")
    for window, lag in [("all_origins", 4), ("legacy_fixed_window", 1)]:
        for columns in [("time", "sigma_xy,sigma_xz"), ("0", "1,2")]:
            output = tmp_path / f"{window}_{columns[0]}.dat"
            proc = run_tool(stress, source, output, "--time-column", columns[0], "--columns", columns[1],
                            "--max-lag", str(lag), "--window", window, "--center", "--threads", "2")
            assert proc.returncode == 0, proc.stderr
            actual = np.loadtxt(output)
            reference = autocorrelation(values[:, 1:], lag, center=True, window=window)
            np.testing.assert_allclose(actual[:, 0], np.arange(lag + 1) * .5)
            np.testing.assert_array_equal(actual[:, 1], reference["origin_counts"])
            np.testing.assert_allclose(actual[:, 2:], reference["values"])


def test_connector_tool_explicit_vector_contract(executables, tmp_path):
    _, chain = executables
    source, output = tmp_path / "q.dat", tmp_path / "q_corr.dat"
    q = np.arange(30.).reshape(5, 6) / 10
    np.savetxt(source, np.column_stack([np.arange(5.), q]))
    proc = run_tool(chain, source, output, "--time-column", "0", "--columns", "1,2,3,4,5,6", "--max-lag", "4")
    assert proc.returncode == 0, proc.stderr
    actual = np.loadtxt(output)
    expected = np.array([
        (q[:5-lag].reshape(-1, 2, 3) * q[lag:].reshape(-1, 2, 3)).sum(axis=2).mean(axis=0)
        for lag in range(5)
    ])
    np.testing.assert_allclose(actual[:, 3:], expected)
    np.testing.assert_allclose(actual[:, 2], expected.mean(axis=1))


def test_tools_reject_malformed_input_options_and_overwrite(executables, tmp_path):
    stress, chain = executables
    source = tmp_path / "input.dat"
    source.write_text("0 1 2 3\n1 2 3 4\n2 3 4 5\n3 4 5 6\n")
    for arguments in [(), ("--time-column", "0"), ("--time-column", "0", "--columns", "4"),
                      ("--time-column", "0", "--columns", "1,1"),
                      ("--time-column", "0", "--columns", "1", "--threads", "0"),
                      ("--time-column", "0", "--columns", "1", "--max-lag", "99")]:
        proc = run_tool(stress, source, tmp_path / "invalid.dat", *arguments)
        assert proc.returncode != 0
        assert not (tmp_path / "invalid.dat").exists()
    assert run_tool(chain, source, tmp_path / "invalid.dat", "--time-column", "0", "--columns", "1,2").returncode != 0
    for content in ["0 1\n1 2 3\n", "0 1\n2 2\n3 3\n", "0 1\n0 2\n", "0 1\n1 nan\n"]:
        source.write_text(content)
        assert run_tool(stress, source, tmp_path / "invalid.dat", "--time-column", "0", "--columns", "1").returncode != 0
    output = tmp_path / "keep.dat"
    output.write_text("original evidence\n")
    assert run_tool(stress, source, output, "--time-column", "0", "--columns", "1").returncode != 0
    assert output.read_text() == "original evidence\n"
