"""Run fixed interacting-Brownian and Hookean/FENE response experiments.

Example (from the package root after installation):
  python examples/two_reference_models.py \
    --output runs/two-reference-v1 --data-budget DATA_BUDGET.json
"""
import argparse
import json
import os
from pathlib import Path

# Set before NumPy/PyTorch imports, even when the launching shell permits more.
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--data-budget", type=Path, required=True)
    parser.add_argument("--seed-manifest", type=Path)
    parser.add_argument("--without-torch", action="store_true")
    args = parser.parse_args()

    from bdmc.reference_experiments import DEFAULT_SEED_MANIFEST, run_two_reference_experiments
    from bdmc.storage import DataBudget

    report = run_two_reference_experiments(args.output, budget=DataBudget(args.data_budget),
        seed_manifest=args.seed_manifest or DEFAULT_SEED_MANIFEST, include_torch=not args.without_torch)
    print(json.dumps({"status": report["status"], "output": str(args.output),
        "elapsed_seconds": report["elapsed_seconds"],
        "models": {name: {"runs": result["n_runs"], "rows": result["n_rows"],
            "test_rmse_by_target": {method: {target: metrics["rmse"]
                for target, metrics in record["metrics"]["test"]["per_target"].items()}
                for method, record in result["methods"].items()}}
            for name, result in report["results"].items()}}, indent=2))


if __name__ == "__main__":
    main()
