"""Method-independent training, evaluation and versioned inference artifacts."""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
from time import perf_counter
from typing import Any, Mapping

import numpy as np

from .data import Dataset, DataSplit, Standardizer, grouped_split, json_copy, numeric_array
from .methods import get_method

ARTIFACT_SCHEMA_VERSION = 1


def regression_metrics(y_true: np.ndarray, prediction: np.ndarray,
                       target_names: tuple[str, ...]) -> dict[str, Any]:
    truth = numeric_array(y_true, ndim=2, name="truth")
    prediction = numeric_array(prediction, ndim=2, name="prediction")
    if truth.shape != prediction.shape or truth.shape[1] != len(target_names):
        raise ValueError("prediction, truth and target names must match")
    residual = prediction - truth
    mse = np.mean(residual ** 2, axis=0)
    denominator = ((truth - truth.mean(axis=0)) ** 2).sum(axis=0)
    per_target = {}
    for i, name in enumerate(target_names):
        per_target[name] = {
            "mse": float(mse[i]), "rmse": float(np.sqrt(mse[i])),
            "mae": float(np.mean(np.abs(residual[:, i]))),
            "r2": None if denominator[i] == 0 else
            float(1 - (residual[:, i] ** 2).sum() / denominator[i])}
    return {"n_rows": len(truth), "rmse": float(np.sqrt(np.mean(residual ** 2))),
            "mae": float(np.mean(np.abs(residual))), "per_target": per_target}


def _numeric_state(state: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    if not state:
        raise ValueError("inference state cannot be empty")
    result = {}
    for name, value in state.items():
        if not isinstance(name, str) or not name or not name.replace("_", "").isalnum():
            raise ValueError("inference array names must be alphanumeric with optional underscores")
        array = np.array(value, dtype=np.float64, copy=True)
        if not array.size or not np.isfinite(array).all():
            raise ValueError("inference state must contain nonempty finite numeric arrays")
        array.setflags(write=False)
        result[name] = array
    return result


def _versions() -> dict[str, str]:
    versions = {"python": platform.python_version(), "numpy": np.__version__}
    for package in ("scikit-learn", "torch", "bdmc"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            continue
    return versions


@dataclass
class TrainedArtifact:
    method: str
    config: dict[str, Any]
    state: dict[str, np.ndarray]
    feature_scaler: Standardizer
    target_scaler: Standardizer
    feature_names: tuple[str, ...]
    target_names: tuple[str, ...]
    provenance: dict[str, Any]
    split: DataSplit
    metrics: dict[str, Any]
    timings: dict[str, float] = field(default_factory=dict)

    @property
    def uncertainty_kind(self) -> str | None:
        return get_method(self.method).uncertainty

    def predict(self, X: np.ndarray, *, return_std: bool = False):
        """Predict target values; optional GP std is latent, not physical noise."""
        X = numeric_array(X, ndim=2, name="prediction features")
        if X.shape[1] != len(self.feature_names):
            raise ValueError("prediction feature width does not match the trained schema")
        method = get_method(self.method)
        if return_std and not method.uncertainty:
            raise ValueError(f"{self.method} does not provide an uncertainty estimate")
        mean, std = method.predict(self.feature_scaler.transform(X), self.state, self.config)
        mean = np.asarray(mean)
        if mean.shape != (len(X), len(self.target_names)) or not np.isfinite(mean).all():
            raise ValueError("method returned invalid prediction shape or values")
        result = self.target_scaler.inverse(mean)
        if not np.isfinite(result).all():
            raise ValueError("prediction overflowed in target units")
        if return_std:
            if method.predict_std is not None:
                std = method.predict_std(self.feature_scaler.transform(X), self.state, self.config)
            if std is None or np.shape(std) != mean.shape or not np.isfinite(std).all() or np.any(std < 0):
                raise ValueError("method returned an invalid standard deviation")
            return result, np.asarray(std) * self.target_scaler.scale
        return result

    def save(self, directory: str | Path) -> Path:
        """Write JSON + NPZ, refusing to overwrite an existing artifact."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        manifest_path = directory / "manifest.json"
        weights_path = directory / "weights.npz"
        if manifest_path.exists() or weights_path.exists():
            raise FileExistsError(f"artifact already exists at {directory}")
        state = _numeric_state(self.state)
        arrays = {f"model__{name}": value for name, value in state.items()}
        arrays.update(feature_mean=self.feature_scaler.mean, feature_scale=self.feature_scaler.scale,
                      target_mean=self.target_scaler.mean, target_scale=self.target_scaler.scale,
                      split_train=self.split.train, split_validation=self.split.validation,
                      split_test=self.split.test)
        np.savez_compressed(weights_path, **arrays)
        # The digest detects changed numeric bytes; it does not authenticate a
        # publisher or prove that the recorded targets are scientifically valid.
        checksum = hashlib.sha256(weights_path.read_bytes()).hexdigest()
        manifest = {"schema_version": ARTIFACT_SCHEMA_VERSION, "method": self.method,
                    "config": self.config, "feature_names": list(self.feature_names),
                    "target_names": list(self.target_names), "provenance": self.provenance,
                    "metrics": self.metrics, "timings": self.timings,
                    "uncertainty_kind": self.uncertainty_kind,
                    "uncertainty_calibrated": False, "weights_sha256": checksum,
                    "split": {"dataset_fingerprint": self.split.dataset_fingerprint,
                              "group_by": list(self.split.group_by), "seed": self.split.seed}}
        manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n",
                                 encoding="utf-8")
        return directory


def train(dataset: Dataset, *, method: str = "ridge", split: DataSplit | None = None,
          config: Mapping[str, Any] | None = None, seed: int = 0,
          evaluate_test: bool = True) -> TrainedArtifact:
    """Fit preprocessing and estimator exclusively on the training partition.

    Validation and optional test metrics are reported separately. Set
    evaluate_test=False while comparing methods or simulation budgets. This function performs
    no hyperparameter search or automatic test-set selection.
    """
    split = grouped_split(dataset, seed=seed) if split is None else split
    split.validate(dataset)
    definition = get_method(method, training=True)
    fit_started = perf_counter()
    # Held-out feature extremes must not influence centering/scales either:
    # preprocessing is learned state and can leak information just like weights.
    feature_scaler = Standardizer.fit(dataset.X[split.train])
    target_scaler = Standardizer.fit(dataset.y[split.train])
    state, effective_config = definition.fit(
        feature_scaler.transform(dataset.X[split.train]),
        target_scaler.transform(dataset.y[split.train]), json_copy(dict(config or {})), seed)
    numeric_state = _numeric_state(state)
    fit_seconds = perf_counter() - fit_started
    artifact = TrainedArtifact(
        method, json_copy(effective_config), numeric_state, feature_scaler, target_scaler,
        dataset.feature_names, dataset.target_names,
        {"dataset_fingerprint": dataset.fingerprint, "dataset_metadata": json_copy(dataset.metadata),
         "task": json_copy(dataset.task.__dict__), "training_seed": int(seed),
         "preprocessing_fit_partition": "train", "versions": _versions(),
         "test_evaluated": bool(evaluate_test),
         "platform": sys.platform,
         "target_uncertainty_note": "Between-realization variation is not estimated by squared-error regression."},
        split, {}, {"fit_seconds": fit_seconds})
    partitions = [("train", split.train), ("validation", split.validation)]
    if evaluate_test:
        partitions.append(("test", split.test))
    artifact.metrics = {name: regression_metrics(dataset.y[indices], artifact.predict(dataset.X[indices]),
                                                 dataset.target_names)
                        for name, indices in partitions}
    return artifact


def compare(dataset: Dataset, methods: Mapping[str, Mapping[str, Any]], *,
            split: DataSplit | None = None, seed: int = 0) -> dict[str, TrainedArtifact]:
    """Fit explicitly named methods on the same frozen split and simulation budget."""
    # This convenience reports each method's test metrics by train's default.
    # For model/budget selection with a sealed test partition use train(...,
    # evaluate_test=False) or evaluate_learning_efficiency's predeclared policy.
    split = grouped_split(dataset, seed=seed) if split is None else split
    return {name: train(dataset, method=name, config=config, split=split, seed=seed)
            for name, config in methods.items()}


def load_artifact(directory: str | Path) -> TrainedArtifact:
    """Load numeric-only state. Unknown versions/methods and corrupt arrays fail."""
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"),
                          parse_constant=lambda value: (_ for _ in ()).throw(
                              ValueError(f"nonfinite JSON value: {value}")))
    if manifest.get("schema_version") != ARTIFACT_SCHEMA_VERSION:
        raise ValueError("unsupported artifact schema version")
    method = get_method(manifest["method"])
    if manifest.get("uncertainty_kind") != method.uncertainty:
        raise ValueError("artifact uncertainty contract differs from the registered method")
    weights_path = directory / "weights.npz"
    if hashlib.sha256(weights_path.read_bytes()).hexdigest() != manifest["weights_sha256"]:
        raise ValueError("artifact weight checksum mismatch")
    with np.load(weights_path, allow_pickle=False) as archive:
        arrays = {name: archive[name].copy() for name in archive.files}
    state = _numeric_state({name.removeprefix("model__"): value for name, value in arrays.items()
                            if name.startswith("model__")})
    feature_names = tuple(manifest["feature_names"])
    target_names = tuple(manifest["target_names"])
    for names in (feature_names, target_names):
        if not names or len(set(names)) != len(names) or any(not isinstance(n, str) or not n for n in names):
            raise ValueError("invalid artifact column names")
    scalers = []
    for prefix, width in (("feature", len(feature_names)), ("target", len(target_names))):
        mean = numeric_array(arrays[f"{prefix}_mean"], ndim=1, name=f"{prefix} mean")
        scale = numeric_array(arrays[f"{prefix}_scale"], ndim=1, name=f"{prefix} scale")
        if mean.shape != (width,) or scale.shape != (width,) or np.any(scale <= 0):
            raise ValueError("invalid artifact scaler shape or values")
        scalers.append(Standardizer(mean, scale))
    split_manifest = manifest["split"]
    partition_arrays = [arrays[f"split_{name}"] for name in ("train", "validation", "test")]
    for indices in partition_arrays:
        if indices.ndim != 1 or indices.dtype.kind not in "iu" or not len(indices):
            raise ValueError("invalid saved split indices")
    all_indices = np.concatenate(partition_arrays)
    if not np.array_equal(np.sort(all_indices), np.arange(len(all_indices))):
        raise ValueError("saved split does not partition its recorded rows")
    split = DataSplit(*partition_arrays, split_manifest["dataset_fingerprint"],
                      tuple(split_manifest["group_by"]), split_manifest["seed"])
    if "groups" not in split.group_by or len(set(split.group_by)) != len(split.group_by):
        raise ValueError("saved split must include distinct grouping keys and source groups")
    artifact = TrainedArtifact(manifest["method"], json_copy(manifest["config"]), state,
                               *scalers, feature_names, target_names, manifest["provenance"],
                               split, manifest["metrics"], manifest.get("timings", {}))
    artifact.predict(np.zeros((1, len(feature_names))))
    return artifact
