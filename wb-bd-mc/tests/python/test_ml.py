"""Independent fixtures for the simulation/ML contract, not scientific regression."""
import hashlib
import importlib.util
import json
from dataclasses import replace

import numpy as np
import pytest

from bdmc.ml import (Dataset, DataSplit, RegressionTask, Standardizer,
                     available_methods, compare, grouped_split, load_artifact,
                     register_method, regression_metrics, train)


@pytest.fixture
def linear_data():
    X = np.linspace(-3, 3, 72)[:, None]
    y = np.column_stack((3 * X[:, 0] - 2, -X[:, 0] + 0.25))
    return Dataset(X, y, np.repeat(np.arange(24), 3), feature_names=("control",),
                   target_names=("first", "second"), metadata={"source": "analytic_fixture"},
                   task=RegressionTask("linear_fixture", "ensemble_mean", ("1", "1")))


def test_grouped_split_is_disjoint_complete_and_reproducible(linear_data):
    a = grouped_split(linear_data, seed=27)
    b = grouped_split(linear_data, seed=27)
    a.validate(linear_data)
    assert np.array_equal(a.train, b.train)
    assert np.array_equal(a.validation, b.validation)
    assert np.array_equal(a.test, b.test)
    assert not np.array_equal(a.train, grouped_split(linear_data, seed=28).train)
    groups = [set(linear_data.groups[getattr(a, key)]) for key in ("train", "validation", "test")]
    assert groups[0].isdisjoint(groups[1]) and groups[0].isdisjoint(groups[2])
    assert groups[1].isdisjoint(groups[2])


def test_overlapping_continuations_are_joined_transitively():
    data = Dataset(np.arange(8)[:, None], np.arange(8),
                   np.array(["a", "a", "b", "c", "d", "e", "f", "g"]),
                   group_keys={"continuation": np.array(["x", "y", "y", "z", "u", "v", "w", "q"])})
    split = grouped_split(data, group_by=("groups", "continuation"))
    for indices in (split.train, split.validation, split.test):
        if 0 in indices:
            assert 1 in indices and 2 in indices
    split.validate(data)


def test_condition_holdout_keeps_replicates_together():
    data = Dataset(np.arange(24)[:, None], np.arange(24), np.arange(24),
                   group_keys={"condition": np.repeat(np.arange(6), 4)})
    split = grouped_split(data, group_by=("groups", "condition"))
    split.validate(data)
    assert sum(len(set(data.labels("condition")[part]))
               for part in (split.train, split.validation, split.test)) == 6


def test_condition_only_request_cannot_split_one_run_across_conditions():
    data = Dataset(np.arange(8)[:, None], np.arange(8),
                   np.array(["a", "a", "b", "b", "c", "d", "e", "f"]),
                   group_keys={"condition": np.array(["x", "y", "y", "z", "u", "v", "w", "q"])})
    split = grouped_split(data, group_by=("condition",))
    assert "groups" in split.group_by
    split.validate(data)
    for indices in (split.train, split.validation, split.test):
        if 0 in indices:
            assert {0, 1, 2, 3}.issubset(set(indices))
    manual = DataSplit(split.train, split.validation, split.test, data.fingerprint, ("condition",))
    with pytest.raises(ValueError, match="source groups"):
        manual.validate(data)


def test_bad_splits_are_rejected(linear_data):
    with pytest.raises(ValueError, match="three nonoverlapping"):
        grouped_split(Dataset([[1], [2], [3]], [1, 2, 3], [0, 0, 1]))
    with pytest.raises(ValueError, match="fractions"):
        grouped_split(linear_data, fractions=(0.8, 0.2, 0.2))
    with pytest.raises(ValueError, match="unknown group"):
        grouped_split(linear_data, group_by=("missing",))
    split = grouped_split(linear_data)
    changed = Dataset(linear_data.X, linear_data.y + 1, linear_data.groups)
    with pytest.raises(ValueError, match="different dataset"):
        split.validate(changed)
    # Complete index coverage is insufficient if adjacent windows share a run.
    leaky = DataSplit(np.arange(0, 72, 3), np.arange(1, 72, 3), np.arange(2, 72, 3),
                      linear_data.fingerprint)
    with pytest.raises(ValueError, match="leakage"):
        leaky.validate(linear_data)


def test_train_only_standardization_with_extreme_heldout_values():
    X = np.array([0, 2, 100, 200, -100, -200])[:, None]
    data = Dataset(X, 2 * X[:, 0] + 1, np.arange(6))
    split = DataSplit(np.array([0, 1]), np.array([2, 3]), np.array([4, 5]), data.fingerprint)
    artifact = train(data, split=split, config={"alpha": 0})
    np.testing.assert_allclose(artifact.feature_scaler.mean, [1])
    np.testing.assert_allclose(artifact.feature_scaler.scale, [1])
    np.testing.assert_allclose(artifact.target_scaler.mean, [3])
    np.testing.assert_allclose(artifact.target_scaler.scale, [2])
    np.testing.assert_allclose(artifact.predict([[3], [1000]]), [[7], [2001]], atol=1e-10)


def test_ridge_matches_exact_independent_linear_solution(linear_data):
    artifact = train(linear_data, config={"alpha": 0})
    np.testing.assert_allclose(artifact.predict([[-1], [0], [2]]),
                               [[-5, 1.25], [-2, 0.25], [4, -1.75]], atol=1e-12)
    assert artifact.metrics["test"]["rmse"] < 1e-12
    assert artifact.uncertainty_kind is None
    with pytest.raises(ValueError, match="uncertainty"):
        artifact.predict([[0]], return_std=True)


def test_ridge_handles_rank_deficiency_and_constant_features():
    X = np.column_stack((np.arange(12), np.arange(12), np.ones(12)))
    data = Dataset(X, 2 * np.arange(12) + 3, np.arange(12))
    artifact = train(data, config={"alpha": 0})
    np.testing.assert_allclose(artifact.predict([[4, 4, 1]]), [[11]], atol=1e-12)
    assert artifact.feature_scaler.scale[2] == 1


def test_metrics_have_known_answers_and_constant_r2_is_undefined():
    metrics = regression_metrics(np.array([[1.], [3.]]), np.array([[2.], [1.]]), ("value",))
    assert metrics["rmse"] == pytest.approx(np.sqrt(2.5))
    assert metrics["mae"] == 1.5
    assert metrics["per_target"]["value"]["r2"] == -1.5
    constant = regression_metrics(np.ones((2, 1)), np.ones((2, 1)), ("value",))
    assert constant["per_target"]["value"]["r2"] is None


def test_artifact_roundtrip_captures_provenance(linear_data, tmp_path):
    artifact = train(linear_data, seed=14)
    destination = artifact.save(tmp_path / "model")
    loaded = load_artifact(destination)
    np.testing.assert_array_equal(loaded.predict(linear_data.X), artifact.predict(linear_data.X))
    assert loaded.provenance["dataset_fingerprint"] == linear_data.fingerprint
    assert loaded.provenance["dataset_metadata"]["source"] == "analytic_fixture"
    assert loaded.provenance["preprocessing_fit_partition"] == "train"
    assert loaded.metrics == artifact.metrics
    np.testing.assert_array_equal(loaded.split.train, artifact.split.train)
    manifest = json.loads((destination / "manifest.json").read_text())
    assert manifest["uncertainty_calibrated"] is False
    with pytest.raises(FileExistsError):
        artifact.save(destination)


def test_artifact_rejects_corruption_and_unknown_schema(linear_data, tmp_path):
    artifact = train(linear_data)
    destination = artifact.save(tmp_path / "corrupt")
    with (destination / "weights.npz").open("ab") as stream:
        stream.write(b"corruption")
    with pytest.raises(ValueError, match="checksum"):
        load_artifact(destination)
    other = artifact.save(tmp_path / "future")
    manifest_path = other / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["schema_version"] = 999
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="schema version"):
        load_artifact(other)


def test_artifact_never_loads_object_arrays(linear_data, tmp_path):
    destination = train(linear_data).save(tmp_path / "object")
    weights = destination / "weights.npz"
    np.savez(weights, model__coefficients=np.array([{"payload": "must not unpickle"}], dtype=object))
    manifest_path = destination / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["weights_sha256"] = hashlib.sha256(weights.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Object arrays cannot be loaded"):
        load_artifact(destination)


@pytest.mark.parametrize("X,y,groups", [([[np.nan], [1], [2]], [1, 2, 3], [0, 1, 2]),
                                        ([[1], [2]], [1], [0, 1]),
                                        ([[1], [2]], [1, 2], [0]),
                                        ([[1], [2]], [1, 2], [None, "a"])])
def test_dataset_rejects_invalid_inputs(X, y, groups):
    with pytest.raises((ValueError, TypeError)):
        Dataset(X, y, groups)


def test_schema_configuration_and_provenance_validation(linear_data):
    with pytest.raises(ValueError, match="feature_names"):
        Dataset([[1], [2], [3]], [1, 2, 3], [0, 1, 2], feature_names=("a", "b"))
    with pytest.raises((ValueError, TypeError)):
        Dataset([[1], [2], [3]], [1, 2, 3], [0, 1, 2], metadata={"bad": float("nan")})
    with pytest.raises(ValueError, match="target_kind"):
        RegressionTask(target_kind="undefined")
    with pytest.raises(ValueError, match="unknown method"):
        train(linear_data, method="not_implemented_gnn")
    with pytest.raises(ValueError, match="configuration"):
        train(linear_data, config={"typo": 1})
    with pytest.raises(ValueError, match="nonnegative"):
        train(linear_data, config={"alpha": -1})
    with pytest.raises(ValueError, match="width"):
        train(linear_data).predict([[1, 2]])


def test_compare_reuses_one_split_and_registry_supports_extensions(linear_data):
    def fit(X, y, config, seed):
        return {"mean": y.mean(axis=0)}, {}

    def predict(X, state, config):
        return np.repeat(state["mean"][None, :], len(X), axis=0), None

    register_method("test_constant", fit=fit, predict=predict)
    artifacts = compare(linear_data, {"ridge": {"alpha": 0}, "test_constant": {}}, seed=3)
    assert available_methods()["test_constant"]["available"]
    np.testing.assert_array_equal(artifacts["ridge"].split.test, artifacts["test_constant"].split.test)
    assert artifacts["ridge"].metrics["test"]["rmse"] < artifacts["test_constant"].metrics["test"]["rmse"]
    with pytest.raises(ValueError, match="already registered"):
        register_method("test_constant", fit=fit, predict=predict)


def test_optional_dependency_error_is_explicit(linear_data, monkeypatch):
    original = importlib.util.find_spec
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None if name == "sklearn" else original(name))
    assert not available_methods()["gaussian_process"]["available"]
    with pytest.raises(ImportError, match="sklearn"):
        train(linear_data, method="gaussian_process")


def test_gaussian_process_matches_independent_kernel_conditioning(linear_data, tmp_path):
    pytest.importorskip("sklearn")
    artifact = train(linear_data, method="gaussian_process", config={"noise_variance": 0.05})
    points = np.array([[-1.25], [0.75]])
    prediction, std = artifact.predict(points, return_std=True)
    X_train = artifact.feature_scaler.transform(linear_data.X[artifact.split.train])
    y_train = artifact.target_scaler.transform(linear_data.y[artifact.split.train])
    query = artifact.feature_scaler.transform(points)
    K = np.exp(-0.5 * ((X_train[:, None] - X_train[None, :]) ** 2).sum(axis=2))
    K += np.eye(len(K)) * 0.05
    cross = np.exp(-0.5 * ((query[:, None] - X_train[None, :]) ** 2).sum(axis=2))
    mean_expected = artifact.target_scaler.inverse(cross @ np.linalg.solve(K, y_train))
    variance = 1 - np.sum(cross * np.linalg.solve(K, cross.T).T, axis=1)
    std_expected = np.sqrt(variance)[:, None] * artifact.target_scaler.scale
    np.testing.assert_allclose(prediction, mean_expected, atol=1e-12)
    np.testing.assert_allclose(std, std_expected, atol=1e-12)
    assert "uncalibrated" in artifact.uncertainty_kind
    loaded = load_artifact(artifact.save(tmp_path / "gp"))
    mean_loaded, std_loaded = loaded.predict(points, return_std=True)
    np.testing.assert_array_equal(mean_loaded, prediction)
    np.testing.assert_array_equal(std_loaded, std)


def test_gp_mean_only_skips_variance_and_std_uses_no_generic_solve(linear_data, monkeypatch):
    pytest.importorskip("sklearn")
    from bdmc.ml import methods

    artifact = train(linear_data, method="gaussian_process")
    definition = methods.get_method("gaussian_process")

    def unexpected(*args, **kwargs):
        raise AssertionError("unnecessary generic solve or variance calculation")

    monkeypatch.setattr(np.linalg, "solve", unexpected)
    monkeypatch.setitem(methods._METHODS, "gaussian_process", replace(definition, predict_std=unexpected))
    assert artifact.predict([[0.5]]).shape == (1, 2)
    monkeypatch.setitem(methods._METHODS, "gaussian_process", definition)
    mean, std = artifact.predict([[0.5]], return_std=True)
    assert mean.shape == std.shape == (1, 2)


def test_torch_mlp_learns_small_problem_and_exports_numpy_inference(tmp_path):
    torch = pytest.importorskip("torch")
    X = np.linspace(-1, 1, 60)[:, None]
    data = Dataset(X, 0.7 * X[:, 0] + 0.2, np.arange(60))
    rng_before = torch.random.get_rng_state().clone()
    artifact = train(data, method="mlp", config={"hidden_layers": [12], "epochs": 200}, seed=8)
    assert torch.equal(rng_before, torch.random.get_rng_state())
    assert artifact.metrics["test"]["rmse"] < 0.04
    points = np.array([[-0.5], [0.5]])
    tensor = torch.tensor(artifact.feature_scaler.transform(points), dtype=torch.float64)
    state = artifact.state
    expected = torch.tanh(tensor @ torch.tensor(state["weight_0"].T) + torch.tensor(state["bias_0"]))
    expected = expected @ torch.tensor(state["weight_1"].T) + torch.tensor(state["bias_1"])
    expected = artifact.target_scaler.inverse(expected.numpy())
    np.testing.assert_allclose(artifact.predict(points), expected, atol=1e-14)
    loaded = load_artifact(artifact.save(tmp_path / "mlp"))
    np.testing.assert_array_equal(loaded.predict(points), artifact.predict(points))
