"""Fixed, bounded experiments for interacting Brownian and standard dumbbells.

These compare separate response estimators; trained models never enter the
integrator. A finite-horizon response is not an equilibrium/transport result.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np

from .ml import (Dataset, RegressionTask, available_methods, compare,
                 grouped_split, load_artifact, regression_metrics)


DEFAULT_SEED_MANIFEST = Path(__file__).with_name("two_reference_seeds.json")
MAX_OUTPUT_BYTES = 128_000_000
MODEL_NAMES = ("repulsive", "hookean", "fene")


@dataclass(frozen=True)
class ReferencePlan:
    """Numerical design frozen before labels or held-out results are observed.

    Smaller subsets support interface tests; the versioned conditions/seeds do
    not change. Bounds prevent this example becoming an unattended large sweep.
    """

    particles: int = 64
    dumbbells: int = 512
    conditions: int = 9
    replicas: int = 3
    steps: int = 1000
    stride: int = 100
    repulsive_warmup: int = 500
    dt: float = .002

    def __post_init__(self):
        limits = {"particles": (8, 512), "dumbbells": (1, 2048),
                  "conditions": (3, 9), "replicas": (1, 3),
                  "steps": (1, 5000), "stride": (1, 5000),
                  "repulsive_warmup": (0, 5000)}
        for name, (low, high) in limits.items():
            value = getattr(self, name)
            if (isinstance(value, (bool, np.bool_))
                    or not isinstance(value, (int, np.integer))
                    or not low <= value <= high):
                raise ValueError(f"{name} must be an integer in [{low}, {high}]")
        side = round(self.particles ** (1 / 3))
        if side ** 3 != self.particles:
            raise ValueError("particles must form a complete cubic initial lattice")
        if not np.isfinite(self.dt) or not 1e-5 <= self.dt <= .005:
            raise ValueError("dt must be finite and in [1e-5, .005]")
        frames = (self.steps + self.stride - 1) // self.stride
        if frames > 50:
            raise ValueError("at most 50 retained response times per run")
        runs = self.conditions * self.replicas
        pair_work = runs * self.particles * (self.particles - 1) // 2 * (
            self.steps + self.repulsive_warmup)
        bead_work = 2 * runs * self.dumbbells * self.steps
        if pair_work + bead_work > 300_000_000:
            raise ValueError("bounded example exceeds its 300 million update-work estimate")


def load_seed_manifest(path=DEFAULT_SEED_MANIFEST) -> dict:
    """Read the versioned identities; randomness is never a predictive input."""
    path = Path(path)
    record = json.loads(path.read_text())
    if record.get("schema") != "bdmc.two-reference-seeds.v1":
        raise ValueError("unsupported two-reference seed manifest")
    all_seeds = []
    for name in MODEL_NAMES:
        parameter = "epsilon" if name == "repulsive" else "spring_constant"
        entry = record[name]
        values = np.asarray(entry[parameter], dtype=float)
        seeds = np.asarray(entry["seeds"])
        if (values.shape != (9,) or not np.isfinite(values).all()
                or np.any(values <= 0) or np.any(np.diff(values) <= 0)
                or seeds.shape != (9, 3) or seeds.dtype.kind not in "iu"
                or np.any(seeds < 0)):
            raise ValueError("manifest requires nine increasing positive conditions and 9x3 integer seeds")
        all_seeds.extend(int(value) for value in seeds.flat)
    reserved = record["reserved_physics_seeds"]
    if any(isinstance(value, bool) or not isinstance(value, int) or value < 0
           for value in reserved):
        raise ValueError("physics seed identities must be nonnegative integers")
    all_seeds.extend(reserved)
    if len(set(all_seeds)) != len(all_seeds):
        raise ValueError("simulation and physics seeds must have distinct identities")
    for key in ("split_seed", "training_seed"):
        if isinstance(record[key], bool) or not isinstance(record[key], int) or record[key] < 0:
            raise ValueError(f"{key} must be a nonnegative integer")
    return {**record, "manifest_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def contact_coordination(positions, *, box_length: float, sigma: float) -> float:
    """Mean number of neighbors with minimum-image distance strictly < sigma."""
    positions = np.asarray(positions, dtype=float)
    if (positions.ndim != 2 or positions.shape[1] != 3 or len(positions) < 2
            or not np.isfinite(positions).all() or not np.isfinite([box_length, sigma]).all()
            or not 0 < sigma <= box_length / 2):
        raise ValueError("finite Nx3 positions and 0 < sigma <= box_length/2 required")
    i, j = np.triu_indices(len(positions), k=1)
    delta = positions[j] - positions[i]
    delta -= box_length * np.rint(delta / box_length)
    count = np.count_nonzero(np.square(delta).sum(axis=1) < sigma ** 2)
    return float(2 * count / len(positions))


def uncorrelated_repulsive_baseline(*, particles: int, box_length: float,
                                   epsilon: float, sigma: float,
                                   diffusion: float, lag: float) -> np.ndarray:
    """Prior comparator: uniform independent positions and ideal diffusion.

    For U=epsilon/3*(1-r/sigma)^2*(2+r/sigma), the ball integral
    is pi*epsilon*sigma^3/9. Pair counting uses (N-1)/V. This is
    deliberately NOT interacting equilibrium theory or an exact MSD law.
    """
    if (isinstance(particles, (bool, np.bool_)) or not isinstance(particles, (int, np.integer))
            or particles < 2 or not np.isfinite([box_length, epsilon, sigma, diffusion, lag]).all()
            or not 0 < sigma <= box_length / 2 or min(epsilon, diffusion, lag) < 0):
        raise ValueError("invalid uncorrelated baseline parameters")
    other_density = (particles - 1) / box_length ** 3
    return np.array([other_density * np.pi * epsilon * sigma ** 3 / 18,
                     other_density * 4 * np.pi * sigma ** 3 / 3,
                     6 * diffusion * lag])


def _source_provenance() -> dict:
    from . import _reference
    from . import reference, reference_theory
    files = {"native_binary": Path(_reference.__file__),
             "python_interface": Path(reference.__file__),
             "theory_oracles": Path(reference_theory.__file__),
             "experiment": Path(__file__)}
    return {name: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for name, path in files.items()}


def _sample_case(simulation, plan, run_id, *, warmup=0):
    try:
        if warmup:
            simulation.advance(warmup)
        return simulation.sample(plan.steps, stride=plan.stride, include_initial=True)
    except Exception as error:
        raise RuntimeError(f"{run_id}: native sampling failed; no replacement seed: {error}") from error


def reference_dataset(model: str, *, plan: ReferencePlan = ReferencePlan(),
                      seed_manifest=DEFAULT_SEED_MANIFEST) -> tuple[Dataset, np.ndarray, dict]:
    """Generate one compact response dataset, preserving all run identities.

    No raw trajectories are written. Times within a run are correlated and
    replicas share one physical condition: both identities must stay together.
    """
    from .reference import (DumbbellSimulation, DumbbellSpec,
                            RepulsiveSimulation, RepulsiveSpec)
    from .reference_theory import hookean_em_moment

    if model not in MODEL_NAMES:
        raise ValueError(f"model must be one of {MODEL_NAMES}")
    manifest = load_seed_manifest(seed_manifest)
    X, y, baseline, groups, conditions, runs = [], [], [], [], [], []
    parameter_name = "epsilon" if model == "repulsive" else "spring_constant"
    started = time.perf_counter()
    total_simulation_seconds = 0.
    for c, parameter in enumerate(manifest[model][parameter_name][:plan.conditions]):
        condition_id = f"{model}-condition-{c:02d}"
        for replica, seed in enumerate(manifest[model]["seeds"][c][:plan.replicas]):
            run_id = f"{condition_id}-replica-{replica:02d}-seed-{seed}"
            simulation_started = time.perf_counter()
            if model == "repulsive":
                # Density remains one for every smaller cubic interface test.
                box_length = float(round(plan.particles ** (1 / 3)))
                spec = RepulsiveSpec(n_particles=plan.particles, box_length=box_length,
                    dt=plan.dt, mobility=1., diffusion=1., epsilon=float(parameter),
                    sigma=1., seed=seed, noise_policy="gaussian", force_search="cell_list")
                simulation = RepulsiveSimulation(spec)
                batch = _sample_case(simulation, plan, run_id, warmup=plan.repulsive_warmup)
                origin = batch["unwrapped_positions"][0]
            else:
                spec = DumbbellSpec(n_dumbbells=plan.dumbbells, dt=plan.dt,
                    mobility=1., diffusion=.5, spring_constant=float(parameter),
                    spring_law=model, q0=4., protocol="simple_shear", shear_rate=.75,
                    seed=seed, noise_policy="gaussian")
                simulation = DumbbellSimulation(spec, np.zeros((plan.dumbbells, 3)))
                batch = _sample_case(simulation, plan, run_id)
            simulation_seconds = time.perf_counter() - simulation_started
            first_step = int(batch["steps"][0])
            first_time = float(batch["times"][0])
            for frame in range(1, len(batch["times"])):
                n = int(batch["steps"][frame]) - first_step
                lag = float(batch["times"][frame]) - first_time
                if model == "repulsive":
                    observed = [float(batch["potential_energy"][frame]) / plan.particles,
                        contact_coordination(batch["positions"][frame], box_length=box_length, sigma=1.),
                        float(np.square(batch["unwrapped_positions"][frame] - origin).sum(axis=1).mean())]
                    prior = uncorrelated_repulsive_baseline(particles=plan.particles,
                        box_length=box_length, epsilon=float(parameter), sigma=1., diffusion=1., lag=lag)
                else:
                    moment, stress = batch["conformation"][frame], batch["polymer_stress"][frame]
                    observed = [float(np.trace(moment)), float(stress[0, 1]),
                                float(stress[0, 0] - stress[1, 1])]
                    gradient = np.zeros((3, 3)); gradient[0, 1] = .75
                    exact_moment = hookean_em_moment(np.zeros((3, 3)), n, plan.dt,
                        velocity_gradient=gradient, diffusion=.5, mobility=1.,
                        spring_constant=float(parameter))
                    prior = [float(np.trace(exact_moment)), parameter * exact_moment[0, 1],
                             parameter * (exact_moment[0, 0] - exact_moment[1, 1])]
                X.append([float(parameter), lag]); y.append(observed); baseline.append(prior)
                groups.append(run_id); conditions.append(condition_id)
            runs.append({"run_id": run_id, "condition_id": condition_id,
                         "replica": replica, "seed": seed, "spec": spec.to_dict(),
                         "metadata": dict(simulation.metadata),
                         "recorded_steps": np.asarray(batch["steps"]).tolist(),
                         "recorded_times": np.asarray(batch["times"]).tolist()})
            # Wall times are kept outside the deterministic dataset identity.
            total_simulation_seconds += simulation_seconds
    if model == "repulsive":
        target_names = ("potential_energy_per_particle", "contact_coordination", "finite_lag_msd")
        units = ("energy", "dimensionless", "length_squared")
        baseline_kind = "uncorrelated uniform-position static moments plus free-diffusion MSD; approximate comparator"
        preparation = f"regular cubic lattice; {plan.repulsive_warmup} warmup Euler steps; stationarity unestablished"
    else:
        target_names = ("connector_second_moment", "polymer_shear_stress_per_dumbbell", "first_normal_stress_per_dumbbell")
        units = ("length_squared", "energy", "energy")
        baseline_kind = ("exact finite-Euler Hookean ensemble moment" if model == "hookean" else
                         "finite-Euler Hookean comparator; not an exact FENE constitutive solution")
        preparation = "zero connectors, simple shear rate .75 from time zero; transient response"
    metadata = {"schema": "bdmc.two-reference-response.v1", "model": model,
        "ensemble_id": manifest["ensemble_id"], "seed_manifest_sha256": manifest["manifest_sha256"],
        "seed_policy": manifest["policy"], "plan": asdict(plan), "runs": runs,
        "source": _source_provenance(), "noise_policy": "gaussian",
        "preparation": preparation, "baseline_kind": baseline_kind,
        "features": "physical parameter and lag from fixed preparation; protocol and all other parameters fixed",
        "label": "finite-ensemble realization; correlated times are not independent samples",
        "observable_conventions": ("energy/N; mean contact coordination=2*pairs(r<sigma)/N; unwrapped displacement from post-warmup origin"
            if model == "repulsive" else "positive Cauchy spring stress <Q outer F>-kT I per dumbbell; multiply by number density for bulk stress"),
        "sampling": "include time zero only as MSD origin; exclude it from training; retain final partial stride",
        "purpose": "bounded response-surrogate methodology; no final network learning target selected"}
    data = Dataset(np.asarray(X), np.asarray(y), np.asarray(groups),
        feature_names=(parameter_name, "observation_lag"), target_names=target_names,
        group_keys={"condition": np.asarray(conditions)}, metadata=metadata,
        task=RegressionTask(name=f"{model}_reference_response", target_kind="realization", target_units=units))
    return data, np.asarray(baseline), {"n_runs": len(runs), "n_rows": len(X),
        "simulation_seconds": total_simulation_seconds,
        "generation_seconds": time.perf_counter() - started}


def _replica_noise(data: Dataset) -> dict:
    """Describe finite-replica dispersion, without counting times as replicas."""
    standard_errors = []
    condition = data.group_keys["condition"]
    for label in np.unique(condition):
        for lag in np.unique(data.X[condition == label, 1]):
            values = data.y[(condition == label) & (data.X[:, 1] == lag)]
            if len(values) > 1:
                standard_errors.append(values.std(axis=0, ddof=1) / np.sqrt(len(values)))
    return {"n_condition_time_cells": len(standard_errors),
        "mean_replica_standard_error": (None if not standard_errors else
            dict(zip(data.target_names, np.mean(standard_errors, axis=0).tolist()))),
        "interpretation": "descriptive SEM from fixed independent replicas at matched condition/time; few replicas, no coverage claim"}


def run_two_reference_experiments(output, *, budget, plan: ReferencePlan = ReferencePlan(),
                                  seed_manifest=DEFAULT_SEED_MANIFEST,
                                  include_torch: bool = True) -> dict:
    """Simulate, split by whole condition, train fixed estimators and reload.

    The caller supplies the existing DataBudget; no baseline is created/reset.
    Failed simulations stop and remain recorded. No reseeding or test tuning.
    """
    from threadpoolctl import threadpool_limits

    # Apply a runtime cap too: a notebook may have imported NumPy before its
    # environment variables were set. threadpoolctl comes with scikit-learn.
    with threadpool_limits(limits=1):
        return _run_two_reference_experiments(output, budget=budget, plan=plan,
            seed_manifest=seed_manifest, include_torch=include_torch)


def _run_two_reference_experiments(output, *, budget, plan, seed_manifest, include_torch):
    from .io import save_dataset
    from .storage import allocated_bytes
    from threadpoolctl import threadpool_info

    output = Path(output).resolve()
    project_root = budget.path.parent
    if output.parent != (project_root / "runs").resolve() or not output.name.startswith("two-reference-"):
        raise ValueError("output must be <budget-root>/runs/two-reference-<new-run-name>")
    manifest = load_seed_manifest(seed_manifest)
    methods = {"ridge": {"alpha": 1e-6},
               "gaussian_process": {"length_scale": 1.5, "noise_variance": .02}}
    optional_status = {}
    if include_torch:
        if available_methods()["mlp"]["available"]:
            import torch
            torch.set_num_threads(1)
            methods["mlp"] = {"hidden_layers": [16, 16], "epochs": 200, "learning_rate": .01}
            optional_status["mlp"] = "included with fixed predeclared configuration"
        else:
            optional_status["mlp"] = "not run: optional torch dependency unavailable"
    else:
        optional_status["mlp"] = "not requested"
    frozen = {"schema": "bdmc.two-reference-experiment.v1", "created_at": datetime.now(timezone.utc).isoformat(),
        "plan": asdict(plan), "seed_manifest": manifest, "methods": methods,
        "source": _source_provenance(),
        "optional_methods": optional_status, "models": MODEL_NAMES,
        "evaluation": "all replicas and times at a parameter condition stay in one partition; fixed split, no tuning",
        "output_cap_bytes": MAX_OUTPUT_BYTES, "platform": platform.platform(),
        "threadpools": threadpool_info(), "native_execution": "sequential CPU reference kernels",
        "physics_validation": "separate fixed-seed scientific checks; ML labels are excluded from acceptance checks"}
    with budget.writer(output, planned_bytes=MAX_OUTPUT_BYTES):
        output.mkdir(parents=True, exist_ok=False)
        (output / "frozen_plan.json").write_text(json.dumps(frozen, indent=2, allow_nan=False) + "\n")
    report = {"status": "running", **frozen, "results": {}}
    started = time.perf_counter()
    try:
        for model in MODEL_NAMES:
            data, baseline, counts = reference_dataset(model, plan=plan, seed_manifest=seed_manifest)
            split = grouped_split(data, seed=manifest["split_seed"], group_by=("groups", "condition"))
            directory = output / model
            save_dataset(directory / "dataset", data, budget=budget)
            artifacts = compare(data, methods, split=split, seed=manifest["training_seed"])
            truth = data.y[split.test]
            prior_metrics = regression_metrics(truth, baseline[split.test], data.target_names)
            mean = np.tile(data.y[split.train].mean(axis=0), (len(split.test), 1))
            inference = {"features": data.X[split.test], "truth": truth,
                         "baseline": baseline[split.test], "test_indices": split.test}
            records = {}
            for method, artifact in artifacts.items():
                with budget.writer(directory / method, planned_bytes=8_000_000):
                    artifact.save(directory / method)
                expected = artifact.predict(data.X[split.test])
                reload_started = time.perf_counter()
                reloaded = load_artifact(directory / method)
                reload_seconds = time.perf_counter() - reload_started
                inference_started = time.perf_counter()
                restored = reloaded.predict(data.X[split.test])
                inference_seconds = time.perf_counter() - inference_started
                np.testing.assert_allclose(restored, expected, rtol=1e-12, atol=1e-12)
                inference[method] = restored
                comparison = {name: artifact.metrics["test"]["per_target"][name]["rmse"]
                              < prior_metrics["per_target"][name]["rmse"] for name in data.target_names}
                records[method] = {"metrics": artifact.metrics, "timings": artifact.timings,
                    "reload_seconds": reload_seconds, "held_out_inference_seconds": inference_seconds,
                    "save_reload_inference_match": True, "lower_test_rmse_than_prior_by_target": comparison}
            result = {**counts, "fingerprint": data.fingerprint, "baseline_kind": data.metadata["baseline_kind"],
                "baseline_metrics": prior_metrics,
                "training_mean_metrics": regression_metrics(truth, mean, data.target_names),
                "split_rows": {name: len(getattr(split, name)) for name in ("train", "validation", "test")},
                "split_conditions": {name: sorted(set(data.group_keys["condition"][getattr(split, name)]))
                                     for name in ("train", "validation", "test")},
                "replica_noise": _replica_noise(data), "methods": records,
                "metric_note": "targets have different units; assess per-target RMSE, not mixed-unit aggregate RMSE"}
            with budget.writer(directory, planned_bytes=1_000_000):
                np.savez_compressed(directory / "held_out_inference.npz", **inference)
                (directory / "report.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
            report["results"][model] = result
        report["status"] = "completed"
    except Exception as error:
        report["status"] = "failed"
        report["failure"] = {"type": type(error).__name__, "message": str(error),
                             "policy": "stopped; no seed replacement or silent dropped run"}
        raise
    finally:
        report["elapsed_seconds"] = time.perf_counter() - started
        report["output_allocated_bytes_before_final_report"] = allocated_bytes(output)
        report["storage_after"] = budget.check()
        report["limits"] = ["small fixed-domain transient response estimators, not a general constitutive closure",
            "no long-time diffusion, equilibrium, finite-size or time-step convergence claim follows from ML labels",
            "Hookean oracle is exact for its Euler second moment only; repulsive/FENE comparators are approximate",
            "no hyperparameter search, held-out-result tuning or adaptive reseeding",
            "saved numeric inference is not a native restart checkpoint; no GPU/MPI backend",
            "latent GP uncertainty is not calibrated physical sampling uncertainty"]
        text = json.dumps(report, indent=2, allow_nan=False) + "\n"
        with budget.writer(output / "report.json", planned_bytes=len(text.encode()) + 4096):
            (output / "report.json").write_text(text)
        if allocated_bytes(output) > MAX_OUTPUT_BYTES:
            raise RuntimeError("reference output exceeded its predeclared 128MB cap")
    return report
