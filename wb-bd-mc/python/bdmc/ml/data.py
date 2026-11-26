"""Validated regression data and leakage-aware partitions."""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from typing import Any, Mapping

import numpy as np


def numeric_array(value: Any, *, ndim: int, name: str) -> np.ndarray:
    # Detach numeric storage from mutable simulation/analysis buffers. A frozen
    # dataclass alone would not prevent a caller from changing an aliased array
    # after the dataset fingerprint or train/test partition was recorded.
    array = np.array(value, dtype=np.float64, copy=True)
    if array.ndim != ndim or not array.size or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a nonempty, finite {ndim}-D numeric array")
    array.setflags(write=False)
    return array


def json_copy(value: Any) -> Any:
    """Reject arbitrary objects, NaNs and unserialisable provenance."""
    return json.loads(json.dumps(value, allow_nan=False))


@dataclass(frozen=True)
class RegressionTask:
    name: str = "regression"
    target_kind: str = "response"
    target_units: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("task name must be nonempty")
        if self.target_kind not in {"response", "ensemble_mean", "realization"}:
            raise ValueError("target_kind must be response, ensemble_mean or realization")
        if any(not isinstance(unit, str) for unit in self.target_units):
            raise ValueError("target_units must contain strings")


@dataclass(frozen=True)
class Dataset:
    """Rows may be correlated; groups must identify independent source lineages.

    Additional group_keys (e.g. condition or continuation) support stricter splits.
    Seeds belong in metadata, not in physical predictive features.
    """

    X: np.ndarray
    y: np.ndarray
    groups: np.ndarray
    feature_names: tuple[str, ...] = ()
    target_names: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)
    task: RegressionTask = field(default_factory=RegressionTask)
    group_keys: Mapping[str, np.ndarray] = field(default_factory=dict)

    def __post_init__(self) -> None:
        X = numeric_array(self.X, ndim=2, name="X")
        y_source = np.asarray(self.y)
        y = numeric_array(y_source[:, None] if y_source.ndim == 1 else y_source,
                          ndim=2, name="y")
        if len(X) != len(y):
            raise ValueError("X and y row counts differ")
        object.__setattr__(self, "X", X)
        object.__setattr__(self, "y", y)
        for attribute, width, prefix in (("feature_names", X.shape[1], "x"),
                                          ("target_names", y.shape[1], "y")):
            names = tuple(getattr(self, attribute)) or tuple(f"{prefix}{i}" for i in range(width))
            if (len(names) != width or len(set(names)) != width
                    or any(not isinstance(n, str) or not n for n in names)):
                raise ValueError(f"{attribute} must have one unique nonempty name per column")
            object.__setattr__(self, attribute, names)
        keys = dict(self.group_keys)
        if "groups" in keys:
            raise ValueError("'groups' is reserved; use the groups argument")
        keys["groups"] = self.groups
        clean_keys = {}
        for name, labels in keys.items():
            if not isinstance(name, str) or not name:
                raise ValueError("group key names must be nonempty strings")
            raw = np.asarray(labels)
            if raw.ndim != 1 or len(raw) != len(X):
                raise ValueError(f"group key {name} must have one label per row")
            if raw.dtype.kind not in "USiu":
                raise ValueError("group labels must be strings or integers, without missing values")
            labels = np.array(raw, dtype=str, copy=True)
            if np.any(labels == ""):
                raise ValueError("group labels cannot be empty")
            labels.setflags(write=False)
            clean_keys[name] = labels
        object.__setattr__(self, "groups", clean_keys.pop("groups"))
        object.__setattr__(self, "group_keys", clean_keys)
        metadata = json_copy(dict(self.metadata))
        object.__setattr__(self, "metadata", metadata)
        if not isinstance(self.task, RegressionTask):
            raise ValueError("task must be a RegressionTask")
        if self.task.target_units and len(self.task.target_units) != y.shape[1]:
            raise ValueError("target_units must match target columns")

    def labels(self, key: str) -> np.ndarray:
        if key == "groups":
            return self.groups
        if key not in self.group_keys:
            raise ValueError(f"unknown group key: {key}")
        return self.group_keys[key]

    @property
    def fingerprint(self) -> str:
        digest = hashlib.sha256()
        for array in (self.X, self.y):
            digest.update(str(array.shape).encode())
            digest.update(array.tobytes())
        for key in sorted(("groups", *self.group_keys)):
            digest.update(json.dumps([key, self.labels(key).tolist()], ensure_ascii=True).encode())
        digest.update(json.dumps({"features": self.feature_names, "targets": self.target_names,
                                  "metadata": self.metadata, "task": self.task.__dict__},
                                 sort_keys=True, allow_nan=False).encode())
        return digest.hexdigest()


@dataclass(frozen=True)
class DataSplit:
    train: np.ndarray
    validation: np.ndarray
    test: np.ndarray
    dataset_fingerprint: str
    group_by: tuple[str, ...] = ("groups",)
    seed: int = 0

    def validate(self, dataset: Dataset) -> None:
        if self.dataset_fingerprint != dataset.fingerprint:
            raise ValueError("split belongs to a different dataset")
        if ("groups" not in self.group_by or len(set(self.group_by)) != len(self.group_by)):
            raise ValueError("split must declare distinct grouping keys including source groups")
        indices = [np.asarray(getattr(self, name)) for name in ("train", "validation", "test")]
        for index in indices:
            if index.ndim != 1 or not index.size or index.dtype.kind not in "iu":
                raise ValueError("partitions must be nonempty one-dimensional integer arrays")
        all_indices = np.concatenate(indices)
        if not np.array_equal(np.sort(all_indices), np.arange(len(dataset.X))):
            raise ValueError("partitions must cover every row exactly once")
        for key in self.group_by:
            seen: set[str] = set()
            labels = dataset.labels(key)
            for index in indices:
                current = set(labels[index])
                if seen & current:
                    raise ValueError(f"group leakage across partitions for {key}")
                seen.update(current)


def independent_components(dataset: Dataset, *,
                           group_by: tuple[str, ...] = ("groups",)) -> tuple[np.ndarray, ...]:
    """Whole source-lineage components joined transitively across selected keys."""
    if not group_by or len(set(group_by)) != len(group_by):
        raise ValueError("group_by must contain distinct grouping keys")
    group_by = tuple(dict.fromkeys(("groups", *group_by)))
    parent = list(range(len(dataset.X)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for key in group_by:
        first: dict[str, int] = {}
        for i, label in enumerate(dataset.labels(key)):
            if label in first:
                parent[find(i)] = find(first[label])
            else:
                first[label] = i
    # Independence is transitive: a continuation shared by two runs must keep
    # both runs together, even when their immediate run labels differ. Splitting
    # time frames randomly would let the same trajectory inform held-out scores.
    component_rows: dict[int, list[int]] = {}
    for i in range(len(parent)):
        component_rows.setdefault(find(i), []).append(i)
    return tuple(np.array(rows, dtype=np.int64) for rows in component_rows.values())


def grouped_split(dataset: Dataset, *, seed: int = 0,
                  fractions: tuple[float, float, float] = (0.6, 0.2, 0.2),
                  group_by: tuple[str, ...] = ("groups",)) -> DataSplit:
    """Split connected group components, keeping all selected keys disjoint.

    Fractions count independent components, not rows. Set group_by to include
    'condition' when measuring new-condition generalisation. Source 'groups' are
    always included, even if omitted from group_by. Overlapping run and
    continuation identifiers are joined transitively before random allocation.
    """
    fractions_array = np.asarray(fractions, dtype=float)
    if (fractions_array.shape != (3,) or not np.isfinite(fractions_array).all()
            or np.any(fractions_array <= 0) or not np.isclose(fractions_array.sum(), 1)):
        raise ValueError("fractions must be three positive values summing to one")
    components = independent_components(dataset, group_by=group_by)
    group_by = tuple(dict.fromkeys(("groups", *group_by)))
    if len(components) < 3:
        raise ValueError("three nonoverlapping independent components are required")
    order = np.random.default_rng(seed).permutation(len(components))
    targets = len(components) * fractions_array
    counts = np.maximum(1, np.floor(targets).astype(int))
    while counts.sum() < len(components):
        counts[np.argmax(targets - counts)] += 1
    while counts.sum() > len(components):
        candidates = np.where(counts > 1, counts - targets, -np.inf)
        counts[np.argmax(candidates)] -= 1
    partitions = []
    start = 0
    for count in counts:
        rows = np.array(sorted(i for j in order[start:start + count] for i in components[j]),
                        dtype=np.int64)
        rows.setflags(write=False)
        partitions.append(rows)
        start += count
    result = DataSplit(*partitions, dataset.fingerprint, tuple(group_by), int(seed))
    result.validate(dataset)
    return result


@dataclass(frozen=True)
class Standardizer:
    mean: np.ndarray
    scale: np.ndarray

    @classmethod
    def fit(cls, values: np.ndarray) -> "Standardizer":
        values = numeric_array(values, ndim=2, name="training values")
        # Identical decimal values can have a tiny nonzero np.std because the
        # floating mean is rounded. Never amplify that roundoff into a feature.
        constant = np.all(values == values[0], axis=0)
        mean = np.where(constant, values[0], values.mean(axis=0))
        scale = values.std(axis=0)
        return cls(mean, np.where(~constant & (scale > 0), scale, 1.0))

    def transform(self, values: np.ndarray) -> np.ndarray:
        return (values - self.mean) / self.scale

    def inverse(self, values: np.ndarray) -> np.ndarray:
        return values * self.scale + self.mean
