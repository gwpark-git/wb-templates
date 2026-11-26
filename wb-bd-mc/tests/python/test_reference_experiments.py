"""Reference experiment semantics, frozen identities and actual native labels."""
import json

import numpy as np
import pytest

from bdmc.ml import grouped_split
from bdmc.reference_experiments import (DEFAULT_SEED_MANIFEST, ReferencePlan,
    contact_coordination, load_seed_manifest, reference_dataset,
    uncorrelated_repulsive_baseline)


def test_fixed_seed_manifest_distinguishes_conditions_replicas_and_physics():
    manifest = load_seed_manifest()
    seeds = [seed for name in ("repulsive", "hookean", "fene")
             for row in manifest[name]["seeds"] for seed in row]
    assert len(seeds) == 81
    assert len(set(seeds)) == 81
    assert not set(seeds) & set(manifest["reserved_physics_seeds"])
    assert manifest["manifest_sha256"]


def test_duplicate_seed_is_rejected_without_replacement(monkeypatch):
    manifest = json.loads(DEFAULT_SEED_MANIFEST.read_text())
    manifest["fene"]["seeds"][0][0] = manifest["repulsive"]["seeds"][0][0]
    monkeypatch.setattr(type(DEFAULT_SEED_MANIFEST), "read_text", lambda path: json.dumps(manifest))
    with pytest.raises(ValueError, match="distinct identities"):
        load_seed_manifest(DEFAULT_SEED_MANIFEST)


@pytest.mark.parametrize("override", [{"particles": 9}, {"steps": True},
    {"replicas": 4}, {"conditions": 2}, {"dt": float("nan")},
    {"steps": 100, "stride": 1}, {"particles": 512, "steps": 5000}])
def test_unbounded_or_ambiguous_plan_rejected(override):
    with pytest.raises(ValueError):
        ReferencePlan(**override)


def test_contact_coordination_uses_periodic_images_and_excludes_self():
    positions = [[.1, 0., 0.], [3.9, 0., 0.], [2., 2., 2.]]
    assert contact_coordination(positions, box_length=4., sigma=.5) == pytest.approx(2 / 3)
    # The cutoff is strict, matching the native repulsive support.
    assert contact_coordination([[0., 0., 0.], [1., 0., 0.]],
                                box_length=4., sigma=1.) == 0.


def test_uncorrelated_energy_baseline_matches_independent_radial_integration():
    n, length, epsilon, sigma = 64, 4., 2.5, 1.
    radius = np.linspace(0., sigma, 10001)
    x = radius / sigma
    potential = epsilon / 3 * (1 - x) ** 2 * (2 + x)
    integral = np.trapezoid(4 * np.pi * radius ** 2 * potential, radius)
    expected_energy = (n - 1) / (2 * length ** 3) * integral
    baseline = uncorrelated_repulsive_baseline(particles=n, box_length=length,
        epsilon=epsilon, sigma=sigma, diffusion=.5, lag=.25)
    assert baseline[0] == pytest.approx(expected_energy, rel=1e-7)
    assert baseline[1] == pytest.approx((n - 1) / length ** 3 * 4 * np.pi / 3)
    assert baseline[2] == pytest.approx(.75)


@pytest.mark.parametrize("model", ["repulsive", "hookean", "fene"])
def test_actual_reference_data_keep_condition_and_replica_lineages(model):
    pytest.importorskip("bdmc._reference")
    plan = ReferencePlan(particles=8, dumbbells=32, conditions=3, replicas=2,
                         steps=25, stride=10, repulsive_warmup=5)
    data, prior, counts = reference_dataset(model, plan=plan)
    assert counts["n_runs"] == 6
    assert counts["n_rows"] == 18  # final five-step partial stride is retained
    assert data.y.shape == prior.shape == (18, 3)
    assert np.all(data.X[:, 1] > 0)
    assert data.feature_names[0] in ("epsilon", "spring_constant")
    assert data.metadata["source"]["native_binary"]["sha256"]
    assert data.metadata["noise_policy"] == "gaussian"
    assert all(run["seed"] for run in data.metadata["runs"])
    assert len(set(data.groups)) == 6
    split = grouped_split(data, seed=143000, group_by=("groups", "condition"))
    labels = [set(data.group_keys["condition"][getattr(split, name)])
              for name in ("train", "validation", "test")]
    assert len(set.union(*labels)) == 3
    assert not labels[0] & labels[1] and not labels[0] & labels[2] and not labels[1] & labels[2]
    if model == "repulsive":
        assert np.all(data.y >= 0)
        assert np.all(data.y[:, 1] <= plan.particles - 1)
    else:
        assert np.all(data.y[:, 0] > 0)
        assert all(run["spec"]["protocol"] == "simple_shear" for run in data.metadata["runs"])
        if model == "hookean":
            assert "exact finite-Euler" in data.metadata["baseline_kind"]
        else:
            assert "not an exact FENE" in data.metadata["baseline_kind"]
