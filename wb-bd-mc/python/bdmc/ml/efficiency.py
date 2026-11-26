"""Paired learning curves using nested, independent simulation components.

Evaluation partitions are frozen before sampling training budgets. The test
partition is never scored unless one method/repeat is designated in advance.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform as system_platform
from time import perf_counter
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from .data import Dataset, DataSplit, independent_components, json_copy, numeric_array
from .methods import get_method
from .platform import TrainedArtifact, regression_metrics, train

EFFICIENCY_SCHEMA_VERSION = 1


def _positive_integer(value: Any, name: str, *, allow_zero: bool = False) -> int:
    if (not isinstance(value, (int, np.integer)) or isinstance(value, (bool, np.bool_))
            or value < (0 if allow_zero else 1)):
        raise ValueError(f"{name} must be {'a nonnegative' if allow_zero else 'a positive'} integer")
    return int(value)


def _stats(values: Sequence[float]) -> dict[str, float | int | None]:
    return {"mean": float(np.mean(values)),
            "sample_std": float(np.std(values, ddof=1)) if len(values) > 1 else None,
            "n_repeats": len(values)}


def _time_prediction(predict: Callable, X: np.ndarray, repetitions: int):
    # An explicit warmup precedes the timed block for methods and callable baselines.
    predict(X)
    started = perf_counter()
    for _ in range(repetitions):
        prediction = predict(X)
    elapsed = (perf_counter() - started) / repetitions
    return prediction, float(elapsed)


def _subset(dataset: Dataset, selected: np.ndarray, frozen: DataSplit,
            repeat: int, budget: int) -> tuple[Dataset, DataSplit]:
    source_rows = np.concatenate((selected, frozen.validation, frozen.test))
    n_train, n_validation = len(selected), len(frozen.validation)
    metadata = dict(dataset.metadata)
    metadata["efficiency_source"] = {
        "dataset_fingerprint": dataset.fingerprint, "source_row_indices": source_rows.tolist(),
        "repeat": repeat, "budget_components": budget,
        "frozen_split_seed": frozen.seed, "group_by": list(frozen.group_by)}
    subset = Dataset(dataset.X[source_rows], dataset.y[source_rows], dataset.groups[source_rows],
                     dataset.feature_names, dataset.target_names, metadata, dataset.task,
                     {key: labels[source_rows] for key, labels in dataset.group_keys.items()})
    split = DataSplit(np.arange(n_train), np.arange(n_train, n_train + n_validation),
                      np.arange(n_train + n_validation, len(source_rows)), subset.fingerprint,
                      frozen.group_by, frozen.seed)
    split.validate(subset)
    return subset, split


def _summaries(records: list[dict], methods: Sequence[str], budgets: Sequence[int],
               target_names: tuple[str, ...], repeats: int):
    summary = []
    paired = []
    for budget in budgets:
        by_method = {}
        for method in methods:
            current = sorted((record for record in records if record["method"] == method
                              and record["budget_components"] == budget), key=lambda record: record["repeat"])
            by_method[method] = current
            row = {"method": method, "budget_components": budget,
                   "fit_seconds": _stats([r["fit_seconds"] for r in current]),
                   "validation_inference_seconds": _stats([r["validation_inference_seconds"] for r in current]),
                   "validation_per_target_rmse": {
                       name: _stats([r["validation_metrics"]["per_target"][name]["rmse"] for r in current])
                       for name in target_names}}
            if current[0].get("baseline_delta_per_target_rmse") is not None:
                row["baseline_delta_per_target_rmse"] = {
                    name: _stats([r["baseline_delta_per_target_rmse"][name] for r in current])
                    for name in target_names}
            summary.append(row)
        for left, right in combinations(methods, 2):
            # Each pair uses identical training rows for each repeat, not independent resamples.
            assert len(by_method[left]) == len(by_method[right]) == repeats
            paired.append({"budget_components": budget, "left": left, "right": right,
                           "difference": "right_minus_left",
                           "validation_per_target_rmse_difference": {
                               name: _stats([
                                   b["validation_metrics"]["per_target"][name]["rmse"]
                                   - a["validation_metrics"]["per_target"][name]["rmse"]
                                   for a, b in zip(by_method[left], by_method[right])])
                               for name in target_names}})
    return summary, paired


@dataclass
class EfficiencyReport:
    metadata: dict[str, Any]
    records: list[dict[str, Any]]
    summary: list[dict[str, Any]]
    paired_differences: list[dict[str, Any]]
    baseline: dict[str, Any] | None
    final_test: dict[str, Any] | None
    artifacts: dict[tuple[str, int, int], TrainedArtifact] = field(repr=False)

    def to_dict(self) -> dict[str, Any]:
        return json_copy({"schema_version": EFFICIENCY_SCHEMA_VERSION, "metadata": self.metadata,
                          "records": self.records, "summary": self.summary,
                          "paired_differences": self.paired_differences,
                          "baseline": self.baseline, "final_test": self.final_test})

    def save(self, directory: str | Path, *, save_artifacts: bool = False) -> Path:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        report_path = directory / "report.json"
        if report_path.exists():
            raise FileExistsError(f"efficiency report already exists at {directory}")
        payload = self.to_dict()
        if save_artifacts:
            payload["artifact_paths"] = []
            for i, ((method, repeat, budget), artifact) in enumerate(self.artifacts.items()):
                relative = f"artifacts/run-{i:04d}"
                artifact.save(directory / relative)
                payload["artifact_paths"].append({"method": method, "repeat": repeat,
                                                   "budget_components": budget, "path": relative})
        report_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        return report_path


def evaluate_learning_efficiency(
        dataset: Dataset, methods: Mapping[str, Mapping[str, Any]], *,
        budgets: Sequence[int], split: DataSplit, repeats: int = 3, seed: int = 0,
        baseline: Callable[[np.ndarray], np.ndarray] | np.ndarray | None = None,
        baseline_name: str = "external_reference", inference_repeats: int = 3,
        source_cost_seconds: Mapping[str, float] | None = None,
        final_test_method: str | None = None, final_test_repeat: int = 0,
        retain_artifacts: bool = False, max_fits: int = 100) -> EfficiencyReport:
    """Compare fixed methods on nested whole-component training subsets.

    Budgets count connected components under split.group_by, always including
    source lineages. Methods share the same subset and seed within each repeat.
    Test predictions are disabled by default. If final_test_method is supplied,
    score only that predesignated method/repeat at the largest supplied budget.
    No method, hyperparameter or repeat is chosen from test performance.

    An external baseline is a feature-only callable or fixed [N,T] predictions
    in source-dataset row order. Source costs, if supplied, are measured total
    seconds per source group, not per frame, and are charged once per lineage.
    Reports retain metrics/source indices only by default, not every fitted model.
    """
    split.validate(dataset)
    repeats = _positive_integer(repeats, "repeats")
    seed = _positive_integer(seed, "seed", allow_zero=True)
    inference_repeats = _positive_integer(inference_repeats, "inference_repeats")
    max_fits = _positive_integer(max_fits, "max_fits")
    final_test_repeat = _positive_integer(final_test_repeat, "final_test_repeat", allow_zero=True)
    if final_test_repeat >= repeats:
        raise ValueError("final_test_repeat must be less than repeats")
    budgets = tuple(_positive_integer(budget, "budget") for budget in budgets)
    if not budgets or tuple(sorted(set(budgets))) != budgets:
        raise ValueError("budgets must be nonempty, unique and strictly increasing")
    if not methods:
        raise ValueError("at least one fixed method is required")
    if len(methods) * len(budgets) * repeats > max_fits:
        raise ValueError("requested fits exceed max_fits; reduce the bounded experiment")
    methods = json_copy(dict(methods))
    for name, config in methods.items():
        get_method(name, training=True)
        if not isinstance(config, dict):
            raise ValueError("each method configuration must be a mapping")
    if final_test_method is not None and final_test_method not in methods:
        raise ValueError("final_test_method must name a method fixed before evaluation")
    components = independent_components(dataset, group_by=split.group_by)
    # A budget purchases whole independent source components, not arbitrarily
    # many correlated time rows. Paired methods must see the same components;
    # otherwise a method difference also confounds different simulation data.
    train_rows = set(split.train.tolist())
    train_components = [component for component in components if int(component[0]) in train_rows]
    if budgets[-1] > len(train_components):
        raise ValueError("a budget exceeds the frozen training pool's independent components")
    source_costs = None
    if source_cost_seconds is not None:
        source_costs = {str(key): float(value) for key, value in source_cost_seconds.items()}
        if (not set(dataset.groups).issubset(source_costs)
                or any(not np.isfinite(v) or v < 0 for v in source_costs.values())):
            raise ValueError("source costs must cover every source group with finite nonnegative seconds")

    def counts(indices: np.ndarray) -> dict[str, Any]:
        return {"rows": len(indices),
                "unique_labels": {key: len(set(dataset.labels(key)[indices])) for key in split.group_by},
                "simulation_seconds": None if source_costs is None else
                float(sum(source_costs[group] for group in set(dataset.groups[indices])))}

    baseline_array = None
    baseline_report = None
    if baseline is not None:
        if not isinstance(baseline_name, str) or not baseline_name:
            raise ValueError("baseline_name must be nonempty")
        if not callable(baseline):
            baseline_array = numeric_array(baseline, ndim=2, name="baseline predictions")
            if baseline_array.shape != dataset.y.shape:
                raise ValueError("baseline predictions must match the full source target shape")
            baseline_prediction = baseline_array[split.validation]
            baseline_seconds = None
        else:
            baseline_prediction, baseline_seconds = _time_prediction(
                baseline, dataset.X[split.validation], inference_repeats)
        baseline_report = {"name": baseline_name, "fitted_by_evaluator": False,
                           "fit_seconds": 0.0, "validation_inference_seconds": baseline_seconds,
                           "validation_metrics": regression_metrics(dataset.y[split.validation],
                                                                       baseline_prediction, dataset.target_names)}

    records, artifacts = [], {}
    final_test = None
    for repeat in range(repeats):
        order = np.random.default_rng(np.random.SeedSequence([seed, repeat])).permutation(len(train_components))
        training_seed = int(np.random.SeedSequence([seed, repeat, 1]).generate_state(1)[0])
        for budget in budgets:
            selected = np.sort(np.concatenate([train_components[i] for i in order[:budget]]))
            subset, subset_split = _subset(dataset, selected, split, repeat, budget)
            subset_id = hashlib.sha256(selected.astype("<i8").tobytes()).hexdigest()
            for method, config in methods.items():
                artifact = train(subset, method=method, split=subset_split, config=config,
                                 seed=training_seed, evaluate_test=False)
                prediction, inference_seconds = _time_prediction(
                    artifact.predict, dataset.X[split.validation], inference_repeats)
                metrics = regression_metrics(dataset.y[split.validation], prediction, dataset.target_names)
                if retain_artifacts:
                    artifacts[(method, repeat, budget)] = artifact
                record = {"method": method, "repeat": repeat, "training_seed": training_seed,
                          "effective_config": json_copy(artifact.config),
                          "budget_components": budget, "subset_id": subset_id,
                          "train_source_rows": selected.tolist(),
                          "validation_source_rows": split.validation.tolist(),
                          "test_source_rows": split.test.tolist(),
                          "counts": {"train": counts(selected), "validation": counts(split.validation),
                                     "test": counts(split.test)},
                          "fit_seconds": artifact.timings["fit_seconds"],
                          "validation_inference_seconds": inference_seconds,
                          "validation_metrics": metrics,
                          "baseline_delta_per_target_rmse": None if baseline_report is None else {
                              name: metrics["per_target"][name]["rmse"]
                              - baseline_report["validation_metrics"]["per_target"][name]["rmse"]
                              for name in dataset.target_names}}
                records.append(record)
                if method == final_test_method and repeat == final_test_repeat and budget == budgets[-1]:
                    prediction, elapsed = _time_prediction(artifact.predict, dataset.X[split.test], inference_repeats)
                    test_metrics = regression_metrics(dataset.y[split.test], prediction, dataset.target_names)
                    artifact.metrics["test"] = test_metrics
                    artifact.provenance["test_evaluated"] = True
                    final_test = {"method": method, "repeat": repeat, "budget_components": budget,
                                  "selection": "externally predesignated; not selected from test metrics",
                                  "metrics": test_metrics, "inference_seconds": elapsed}
    if final_test is not None and baseline_report is not None:
        if baseline_array is not None:
            prediction, elapsed = baseline_array[split.test], None
        else:
            prediction, elapsed = _time_prediction(baseline, dataset.X[split.test], inference_repeats)
        baseline_report["final_test_metrics"] = regression_metrics(dataset.y[split.test], prediction, dataset.target_names)
        baseline_report["final_test_inference_seconds"] = elapsed
    summary, paired = _summaries(records, tuple(methods), budgets, dataset.target_names, repeats)
    metadata = {"dataset_fingerprint": dataset.fingerprint, "dataset_metadata": dataset.metadata,
                "feature_names": list(dataset.feature_names), "target_names": list(dataset.target_names),
                "task": dataset.task.__dict__, "group_by": list(split.group_by),
                "frozen_split_seed": split.seed, "sampling_seed": seed, "repeats": repeats,
                "budgets_components": list(budgets), "available_train_components": len(train_components),
                "methods": methods, "inference_repeats": inference_repeats,
                "retained_artifacts": bool(retain_artifacts), "max_fits": max_fits,
                "versions": artifact.provenance["versions"],
                "hardware": {"machine": system_platform.machine(), "processor": system_platform.processor()},
                "timing_policy": "fit includes train-only preprocessing, estimator fit and numeric export; inference is mean warm-batch time",
                "simulation_cost_policy": "caller-supplied total seconds charged once per selected source lineage; null if unknown",
                "evaluation_policy": "validation for curves; test only externally predesignated method/repeat at largest budget",
                "frozen_validation_source_rows": split.validation.tolist(),
                "frozen_test_source_rows": split.test.tolist(),
                "final_test_requested": final_test_method is not None}
    return EfficiencyReport(json_copy(metadata), records, summary, paired, baseline_report, final_test, artifacts)
