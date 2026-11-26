"""Exact same-backend restart, RNG-cache preservation and strict state validation."""

import copy
import json

import numpy as np
import pytest

from bdmc import _core


def assert_identical(left, right):
    a, b = left.observe(), right.observe()
    for name in ("positions", "unwrapped_positions", "forces", "particle_ids", "fixed_bonds", "dynamic_bonds"):
        np.testing.assert_array_equal(a[name], b[name])
    assert a["potential_energy"] == b["potential_energy"]
    assert a["time"] == b["time"]
    assert a["step"] == b["step"]
    assert left.metadata == right.metadata
    assert left.checkpoint()["rng"] == right.checkpoint()["rng"]


@pytest.mark.parametrize("model,policy,n,dimension,prefix,total", [
    ("brownian", "legacy_uniform", 7, 1, 17, 93),
    ("brownian", "gaussian", 3, 1, 0, 60),
    ("brownian", "gaussian", 3, 1, 1, 60),
    ("brownian", "gaussian", 3, 1, 2, 60),
    ("brownian", "gaussian", 3, 1, 7, 60),
    ("dumbbell", "legacy_uniform", 4, 2, 11, 70),
    ("dumbbell", "gaussian", 4, 3, 13, 70),
])
@pytest.mark.parametrize("explicit_positions", [False, True])
def test_exact_continuation_including_gaussian_spare(
    model, policy, n, dimension, prefix, total, explicit_positions,
):
    spec = _core.SimulationSpec(
        model=model, noise_policy=policy, n_particles=n, dimension=dimension,
        box_lengths=list(np.arange(10, 10 + dimension, dtype=float)),
        dt=0.003, mobility=0.9, diffusivity=0.7, spring_alpha=1.1, seed=1997,
    )
    initial = (np.arange(n * dimension).reshape(n, dimension) * 0.13 - 0.2
               if explicit_positions else None)
    ids = np.arange(n, dtype=np.uint64) * 19 + 100
    continuous = _core.Simulation(spec, initial, ids)
    split = _core.Simulation(spec, initial, ids)
    continuous.advance(total)
    split.advance(prefix)
    snapshot = split.checkpoint()
    assert snapshot["schema"] == "bdmc.checkpoint"
    assert snapshot["version"] == 1
    assert snapshot["geometry"]["box_lengths"] == spec.box_lengths
    assert snapshot["initialization"]["explicit_initial_positions"] == explicit_positions
    assert snapshot["rng"]["has_spare"] == (policy == "gaussian" and (prefix * n * dimension) % 2 == 1)
    assert len(snapshot["rng"]["engine_state"].split()) == 313
    restored = _core.Simulation.from_checkpoint(snapshot)
    assert_identical(split, restored)
    split.advance(total - prefix)
    restored.sample(total - prefix, stride=9, include_initial=False)
    assert_identical(continuous, split)  # snapshot did not draw noise
    assert_identical(continuous, restored)


def test_snapshot_owns_copied_arrays_and_metadata_is_json_compatible():
    spec = _core.SimulationSpec(n_particles=3, dimension=1, noise_policy="gaussian", seed=17)
    original = _core.Simulation(spec)
    original.advance(3)
    snapshot = original.checkpoint()
    restored = _core.Simulation.from_checkpoint(snapshot)
    for name, dtype in (("particle_ids", np.uint64), ("positions", np.float64), ("unwrapped_positions", np.float64)):
        array = snapshot["state"][name]
        assert array.dtype == np.dtype(dtype)
        assert array.flags.owndata and array.flags.c_contiguous
    json_metadata = copy.deepcopy(snapshot)
    for section, names in (("state", ["particle_ids", "positions", "unwrapped_positions"]),
                            ("topology", ["fixed_bonds", "dynamic_bonds"])):
        for name in names:
            json_metadata[section].pop(name)
    assert json.loads(json.dumps(json_metadata))["version"] == 1
    snapshot["state"]["positions"][:] = 123
    snapshot["state"]["particle_ids"][:] = 42
    snapshot["spec"]["dt"] = 9
    assert_identical(original, restored)
    original.advance(7)
    restored.advance(7)
    assert_identical(original, restored)


def test_simple_dynamic_graph_state_is_preserved_without_a_new_update_law():
    spec = _core.SimulationSpec(model="dumbbell", n_particles=4, dimension=1)
    source = _core.Simulation(spec, np.array([[0.1], [0.2], [0.4], [0.8]]), [10, 20, 30, 40])
    snapshot = source.checkpoint()
    snapshot["topology"]["dynamic_bonds"] = np.array([[10, 20], [10, 40]], dtype=np.uint64)
    restored = _core.Simulation.from_checkpoint(snapshot)
    restored.advance(10)
    np.testing.assert_array_equal(restored.dynamic_bonds, [[10, 20], [10, 40]])
    np.testing.assert_array_equal(restored.fixed_bonds, [[10, 30], [20, 40]])
    np.testing.assert_array_equal(restored.checkpoint()["topology"]["dynamic_bonds"], [[10, 20], [10, 40]])


def corruption_snapshot():
    spec = _core.SimulationSpec(model="dumbbell", n_particles=4, dimension=1, seed=1991)
    original = _core.Simulation(spec, np.array([[0.1], [0.2], [0.4], [0.8]]), [10, 20, 30, 40])
    original.advance(7)
    return original, original.checkpoint()


CORRUPT_FIELDS = [
    (("schema",), "different"),
    (("version",), 999),
    (("version",), True),
    (("core_version",), "future"),
    (("backend",), "gpu"),
    (("geometry", "boundary"), "shear"),
    (("geometry", "box_lengths"), [20]),
    (("spec", "n_particles"), 3),
    (("spec", "dimension"), 2),
    (("spec", "dt"), -1),
    (("spec", "noise_policy"), "gaussian"),
    (("spec", "model"), "brownian"),
    (("state", "step"), 8),
    (("state", "step"), -1),
    (("state", "step"), True),
    (("state", "time"), 0.5),
    (("state", "time"), np.nan),
    (("state", "positions"), np.full((4, 1), np.nan)),
    (("state", "positions"), np.full((4, 1), 10.0)),
    (("state", "positions"), np.zeros((4,))),
    (("state", "positions"), np.zeros((4, 1), dtype=np.float32)),
    (("state", "unwrapped_positions"), np.full((4, 1), np.inf)),
    (("state", "unwrapped_positions"), np.full((4, 1), 0.12345)),
    (("state", "particle_ids"), np.array([10, 10, 30, 40], dtype=np.uint64)),
    (("state", "particle_ids"), np.array([10, 20, 30, 40], dtype=np.int64)),
    (("state", "particle_ids"), [-1, 20, 30, 40]),
    (("topology", "fixed_bonds"), np.empty((0, 2), dtype=np.uint64)),
    (("topology", "fixed_bonds"), np.array([[10, 20], [30, 40]], dtype=np.uint64)),
    (("topology", "dynamic_bonds"), np.array([[10, 999]], dtype=np.uint64)),
    (("topology", "dynamic_bonds"), np.array([[10, 10]], dtype=np.uint64)),
    (("topology", "dynamic_bonds"), np.array([[10, 20], [10, 20]], dtype=np.uint64)),
    (("topology", "dynamic_bonds"), np.array([[20, 10]], dtype=np.uint64)),
    (("topology", "dynamic_bonds"), np.array([[10, 20]], dtype=np.int64)),
    (("rng", "engine"), "unknown"),
    (("rng", "engine_state"), "not an engine"),
    (("rng", "engine_state"), "1 2 3"),
    (("rng", "engine_state"), "0 " * 313),
    (("rng", "noise_policy"), "gaussian"),
    (("rng", "has_spare"), True),
    (("rng", "has_spare"), 1),
    (("rng", "spare"), np.nan),
    (("initialization", "explicit_initial_positions"), 1),
]


@pytest.mark.parametrize("path,value", CORRUPT_FIELDS)
def test_corrupted_snapshot_rejected_without_mutating_existing_simulation(path, value):
    original, snapshot = corruption_snapshot()
    before = original.observe()
    destination = snapshot
    for key in path[:-1]:
        destination = destination[key]
    destination[path[-1]] = value
    with pytest.raises(ValueError):
        _core.Simulation.from_checkpoint(snapshot)
    after = original.observe()
    np.testing.assert_array_equal(before["positions"], after["positions"])
    assert after["step"] == 7


def test_missing_unknown_fields_and_rng_cursor_errors_rejected():
    _, good = corruption_snapshot()
    mutations = []
    missing = copy.deepcopy(good)
    missing["state"].pop("time")
    mutations.append(missing)
    custom = copy.deepcopy(good)
    custom["components"] = {"learned_force": "cannot silently restore"}
    mutations.append(custom)
    extra_rng = copy.deepcopy(good)
    extra_rng["rng"]["engine_state"] += " 0"
    mutations.append(extra_rng)
    invalid_cursor = copy.deepcopy(good)
    words = invalid_cursor["rng"]["engine_state"].split()
    words[-1] = "313"
    invalid_cursor["rng"]["engine_state"] = " ".join(words)
    mutations.append(invalid_cursor)
    for snapshot in mutations:
        with pytest.raises(ValueError):
            _core.Simulation.from_checkpoint(snapshot)


def test_json_list_numeric_arrays_restore_with_the_same_contract():
    original, snapshot = corruption_snapshot()
    for section, names in (("state", ["particle_ids", "positions", "unwrapped_positions"]),
                            ("topology", ["fixed_bonds", "dynamic_bonds"])):
        for name in names:
            snapshot[section][name] = snapshot[section][name].tolist()
    restored = _core.Simulation.from_checkpoint(json.loads(json.dumps(snapshot)))
    original.advance(5)
    restored.advance(5)
    assert_identical(original, restored)
