"""Learning budgets count independent components, not correlated frames."""
import json

import numpy as np
import pytest

from bdmc.ml import (Dataset, DataSplit, TrainedArtifact, evaluate_learning_efficiency,
                     Standardizer, grouped_split, load_artifact, register_method, train)
from bdmc.ml.data import independent_components
from bdmc.ml import efficiency, platform


@pytest.fixture
def experiment():
    X, y, runs, conditions, continuations = [], [], [], [], []
    for condition in range(12):
        for replicate in range(2):
            for window in range(1 + condition % 3):
                X.append([condition / 11, window / 3])
                y.append([2 * condition / 11 + 3 * window / 3 + 1])
                runs.append(f"run-{condition}-{replicate}")
                conditions.append(f"condition-{condition}")
                continuations.append(f"lineage-{condition}")
    data = Dataset(X, y, np.array(runs), feature_names=("control", "lag"), target_names=("response",),
                   group_keys={"condition": np.array(conditions), "continuation": np.array(continuations)})
    split = grouped_split(data, seed=10, group_by=("groups", "condition", "continuation"))
    return data, split


def test_nested_whole_component_budgets_pair_methods_and_freeze_evaluation(experiment):
    data, split = experiment
    report = evaluate_learning_efficiency(data, {"ridge": {"alpha": 0}}, budgets=(1, 2, 4),
                                          split=split, repeats=3, seed=8, retain_artifacts=True)
    assert len(report.records) == 9
    components = independent_components(data, group_by=split.group_by)
    for repeat in range(3):
        previous = set()
        for budget in (1, 2, 4):
            record = next(r for r in report.records if r["repeat"] == repeat and r["budget_components"] == budget)
            selected = set(record["train_source_rows"])
            assert previous.issubset(selected)
            assert selected.issubset(set(split.train))
            assert selected.isdisjoint(split.validation) and selected.isdisjoint(split.test)
            assert sum(set(component).issubset(selected) for component in components) == budget
            assert record["validation_source_rows"] == split.validation.tolist()
            assert record["test_source_rows"] == split.test.tolist()
            assert record["counts"]["train"]["rows"] == len(selected)
            assert record["counts"]["train"]["unique_labels"]["condition"] == budget
            assert record["counts"]["train"]["unique_labels"]["groups"] == 2 * budget
            artifact = report.artifacts[("ridge", repeat, budget)]
            np.testing.assert_allclose(artifact.feature_scaler.mean, data.X[sorted(selected)].mean(axis=0))
            assert "test" not in artifact.metrics
            previous = selected


def test_same_seed_reproduces_subsets_predictions_and_paired_differences(experiment):
    data, split = experiment

    def fit(X, y, config, seed):
        return {"mean": y.mean(axis=0)}, {}

    def predict(X, state, config):
        return np.tile(state["mean"], (len(X), 1)), None

    register_method("efficiency_constant", fit=fit, predict=predict)
    args = dict(budgets=(2, 4), split=split, repeats=2, seed=19)
    methods = {"ridge": {"alpha": 0}, "efficiency_constant": {}}
    a = evaluate_learning_efficiency(data, methods, **args)
    b = evaluate_learning_efficiency(data, methods, **args)
    for ra, rb in zip(a.records, b.records):
        assert ra["train_source_rows"] == rb["train_source_rows"]
        assert ra["validation_metrics"] == rb["validation_metrics"]
        assert ra["training_seed"] == rb["training_seed"]
    for repeat in range(2):
        for budget in (2, 4):
            records = [r for r in a.records if r["repeat"] == repeat and r["budget_components"] == budget]
            assert records[0]["subset_id"] == records[1]["subset_id"]
            assert records[0]["train_source_rows"] == records[1]["train_source_rows"]
    for difference in a.paired_differences:
        expected = []
        for repeat in range(2):
            current = [r for r in a.records if r["repeat"] == repeat
                       and r["budget_components"] == difference["budget_components"]]
            expected.append(current[1]["validation_metrics"]["per_target"]["response"]["rmse"]
                            - current[0]["validation_metrics"]["per_target"]["response"]["rmse"])
        stats = difference["validation_per_target_rmse_difference"]["response"]
        assert stats["mean"] == pytest.approx(np.mean(expected))
        assert stats["sample_std"] == pytest.approx(np.std(expected, ddof=1))
    assert not a.artifacts  # Numeric models are not retained by default.


def test_transitive_source_and_condition_component_is_never_partially_sampled():
    data = Dataset(np.arange(12)[:, None], np.arange(12),
                   np.array(["a", "a", "b", "b", "c", "d", "e", "f", "g", "h", "i", "j"]),
                   group_keys={"condition": np.array(["x", "y", "y", "z", "c", "d", "e", "f", "g", "h", "i", "j"])})
    split = DataSplit(np.arange(8), np.array([8, 9]), np.array([10, 11]), data.fingerprint,
                      ("groups", "condition"))
    report = evaluate_learning_efficiency(data, {"ridge": {}}, budgets=(1, 3, 5), split=split, repeats=2)
    for record in report.records:
        selected = set(record["train_source_rows"])
        assert not selected.intersection(range(4)) or set(range(4)).issubset(selected)
    assert report.metadata["available_train_components"] == 5


def test_test_scoring_is_disabled_and_optional_final_is_one_predesignated_fit(experiment, monkeypatch):
    data, split = experiment
    test_features = data.X[split.test]
    original = TrainedArtifact.predict
    test_calls = []

    def record(self, X, **kwargs):
        if np.array_equal(X, test_features):
            test_calls.append((self.method, self.provenance["dataset_metadata"]["efficiency_source"]))
        return original(self, X, **kwargs)

    monkeypatch.setattr(TrainedArtifact, "predict", record)
    a = evaluate_learning_efficiency(data, {"ridge": {}}, budgets=(1, 3), split=split, repeats=2)
    assert not test_calls
    assert a.final_test is None
    b = evaluate_learning_efficiency(data, {"ridge": {}}, budgets=(1, 3), split=split, repeats=2,
                                     final_test_method="ridge", final_test_repeat=1,
                                     inference_repeats=2, retain_artifacts=True)
    assert len(test_calls) == 3  # One warmup + two timed calls; exactly one fitted model.
    assert all(method == "ridge" and source["repeat"] == 1 and source["budget_components"] == 3
               for method, source in test_calls)
    assert b.final_test["method"] == "ridge" and b.final_test["repeat"] == 1
    for key, artifact in b.artifacts.items():
        assert ("test" in artifact.metrics) == (key == ("ridge", 1, 3))


def test_cost_accounting_times_fit_and_inference_separately_and_charges_source_once(experiment, monkeypatch):
    data, split = experiment
    tick = [0.0]
    monkeypatch.setattr(platform, "perf_counter", lambda: tick[0])
    monkeypatch.setattr(efficiency, "perf_counter", lambda: tick[0])

    def fit(X, y, config, seed):
        tick[0] += 7.0
        return {"mean": y.mean(axis=0)}, {}

    def predict(X, state, config):
        tick[0] += 3.0
        return np.tile(state["mean"], (len(X), 1)), None

    def baseline(X):
        tick[0] += 2.0
        return 2 * X[:, :1] + 3 * X[:, 1:2] + 1

    register_method("efficiency_clock", fit=fit, predict=predict)
    source_costs = {group: i + 1.0 for i, group in enumerate(sorted(set(data.groups)))}
    report = evaluate_learning_efficiency(data, {"efficiency_clock": {}}, budgets=(1, 2), split=split,
                                          repeats=1, baseline=baseline, source_cost_seconds=source_costs)
    for record in report.records:
        assert record["fit_seconds"] == 7.0
        assert record["validation_inference_seconds"] == 3.0
        selected = record["train_source_rows"]
        expected = sum(source_costs[group] for group in set(data.groups[selected]))
        assert record["counts"]["train"]["simulation_seconds"] == expected
    assert report.baseline["fit_seconds"] == 0
    assert report.baseline["validation_inference_seconds"] == 2.0
    assert report.baseline["validation_metrics"]["rmse"] < 1e-14
    assert report.summary[0]["fit_seconds"]["sample_std"] is None


def test_external_prediction_baseline_metrics_and_saved_report_roundtrip(experiment, tmp_path):
    data, split = experiment
    reference = 2 * data.X[:, :1] + 3 * data.X[:, 1:2] + 1
    report = evaluate_learning_efficiency(data, {"ridge": {"alpha": 0}}, budgets=(2,), split=split,
                                          repeats=1, baseline=reference, baseline_name="analytic_linear",
                                          retain_artifacts=True, final_test_method="ridge")
    assert report.baseline["validation_inference_seconds"] is None
    assert report.baseline["validation_metrics"]["rmse"] == 0
    assert report.baseline["final_test_metrics"]["rmse"] == 0
    saved = report.save(tmp_path / "report", save_artifacts=True)
    payload = json.loads(saved.read_text())
    assert payload["metadata"]["frozen_test_source_rows"] == split.test.tolist()
    assert len(payload["artifact_paths"]) == 1
    artifact = load_artifact(saved.parent / payload["artifact_paths"][0]["path"])
    np.testing.assert_allclose(artifact.predict(data.X[split.validation]),
                               report.artifacts[("ridge", 0, 2)].predict(data.X[split.validation]))
    assert artifact.timings["fit_seconds"] >= 0
    with pytest.raises(FileExistsError):
        report.save(saved.parent)


def test_train_can_seal_test_without_changing_default_contract(experiment):
    data, split = experiment
    assert "test" in train(data, split=split).metrics
    sealed = train(data, split=split, evaluate_test=False)
    assert "test" not in sealed.metrics
    assert sealed.provenance["test_evaluated"] is False


def test_identical_decimal_feature_is_constant_despite_variance_roundoff():
    values = np.full((45, 1), 0.3)
    # NumPy's floating mean can give a nonzero numerical std for identical inputs.
    scaler = Standardizer.fit(values)
    np.testing.assert_array_equal(scaler.mean, [0.3])
    np.testing.assert_array_equal(scaler.scale, [1.0])
    np.testing.assert_array_equal(scaler.transform(values), np.zeros_like(values))
    np.testing.assert_allclose(scaler.transform([[0.6]]), [[0.3]], atol=1e-15)


def test_single_condition_budget_does_not_amplify_constant_parameter_roundoff():
    controls = np.repeat([0.3, 0.6, 0.9], 45)
    time = np.tile(np.linspace(0.01, 1, 45), 3)
    data = Dataset(np.column_stack((controls, time)), 2 * time + 1, np.arange(135),
                   group_keys={"condition": controls.astype(str)})
    split = DataSplit(np.arange(45), np.arange(45, 90), np.arange(90, 135), data.fingerprint,
                      ("groups", "condition"))
    report = evaluate_learning_efficiency(data, {"ridge": {"alpha": 0}}, budgets=(1,),
                                          split=split, repeats=1)
    assert report.records[0]["validation_metrics"]["rmse"] < 1e-12


@pytest.mark.parametrize("option", [dict(budgets=(0,)), dict(budgets=(1.5,)), dict(budgets=(True,)),
                                    dict(budgets=(2, 1)), dict(budgets=(1, 1)), dict(budgets=(100,)),
                                    dict(repeats=0), dict(inference_repeats=0), dict(seed=-1),
                                    dict(final_test_method="not_predeclared"), dict(final_test_repeat=3),
                                    dict(max_fits=1)])
def test_invalid_or_unbounded_design_is_rejected(experiment, option):
    data, split = experiment
    arguments = dict(budgets=(1, 2), split=split, repeats=2)
    arguments.update(option)
    with pytest.raises(ValueError):
        evaluate_learning_efficiency(data, {"ridge": {}}, **arguments)


def test_missing_source_cost_and_wrong_baseline_shape_are_rejected(experiment):
    data, split = experiment
    with pytest.raises(ValueError, match="source costs"):
        evaluate_learning_efficiency(data, {"ridge": {}}, budgets=(1,), split=split,
                                     source_cost_seconds={data.groups[0]: 1})
    with pytest.raises(ValueError, match="baseline predictions"):
        evaluate_learning_efficiency(data, {"ridge": {}}, budgets=(1,), split=split,
                                     baseline=np.ones((2, 1)))


def test_learning_curves_reject_an_existing_leaky_split(experiment):
    data, split = experiment
    leaky = DataSplit(np.arange(0, len(data.X), 3), np.arange(1, len(data.X), 3),
                      np.arange(2, len(data.X), 3), data.fingerprint, split.group_by)
    with pytest.raises(ValueError, match="leakage"):
        evaluate_learning_efficiency(data, {"ridge": {}}, budgets=(1,), split=leaky)
