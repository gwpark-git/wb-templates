"""Actual native MPI execution checks, usable directly by CTest.

No timing/speedup claim follows from these short rank-count comparisons.
"""
from __future__ import annotations

import argparse
from contextlib import nullcontext
from dataclasses import replace
import hashlib
import importlib.machinery
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

# This is a standalone CTest/unittest launcher with required CLI arguments.
# pytest covers the Python launcher in tests/python/test_mpi.py; collecting
# this command-line harness again would run it without its parsed ARGS.
__test__ = False

SOURCE = Path(__file__).resolve().parents[2]
# Editable-install meta finders can precede sys.path and select another checkout.
# Pin this launcher's package only; leave third-party imports to Python normally.
class _SourcePackageFinder:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "bdmc":
            return importlib.machinery.PathFinder.find_spec(fullname, [str(SOURCE/"python")])
        if fullname.startswith("bdmc."):
            return importlib.machinery.PathFinder.find_spec(fullname, path)


def _import_source_helpers():
    # pytest imports this file even with __test__=False. Keep collection inert:
    # only the standalone launcher may select this source package. A broader
    # validator can then retain its own staged package and compiled extensions.
    sys.meta_path.insert(0, _SourcePackageFinder())
    from bdmc.mpi import MPIJob, _bounded_process, serialize_jobs
    from bdmc.storage import DataBudget
    for name in ("bdmc", "bdmc.mpi", "bdmc.storage"):
        if not Path(sys.modules[name].__file__).resolve().is_relative_to(SOURCE/"python/bdmc"):
            raise RuntimeError(f"MPI test imported a different checkout: {name}")
    return MPIJob, _bounded_process, serialize_jobs, DataBudget


class NativeMPI(unittest.TestCase):
    @classmethod
    def writer(cls, path, planned_bytes):
        # The validation controller can already own this same cross-process
        # lock. Its explicit marker avoids taking a non-reentrant child flock.
        if ARGS.budget and os.environ.get("BDMC_OUTER_BUDGET_MANAGED") == "1":
            return nullcontext(cls.budget.require_output(path, planned_bytes))
        return cls.budget.writer(path, planned_bytes=planned_bytes)

    @classmethod
    def setUpClass(cls):
        cls.root = Path(ARGS.scratch_root or tempfile.gettempdir()).resolve()
        if ARGS.budget:
            # A caller can use an existing shared allowance, never reset it.
            cls.budget = DataBudget(ARGS.budget)
            with cls.writer(cls.root, planned_bytes=1_000_000):
                cls.root.mkdir(parents=True, exist_ok=True)
                cls.temporary = tempfile.TemporaryDirectory(prefix="mpi-native-test-", dir=cls.root)
        else:
            # CTest needs no research workspace or production data allowance.
            # This disposable fixture covers only tiny test manifests/logs.
            cls.root.mkdir(parents=True, exist_ok=True)
            cls.temporary = tempfile.TemporaryDirectory(prefix="mpi-native-test-", dir=cls.root)
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory = Path(cls.temporary.name)
        if not ARGS.budget:
            baseline = cls.directory/"test-only-budget.json"
            baseline.write_text(json.dumps({
                "schema": "bdmc.additional-data-budget.v1", "created_at": "disposable-mpi-test-only",
                "baseline_allocated_bytes": {str(cls.directory): 0},
                "additional_limit_bytes": 2_000_000, "soft_limit_bytes": 1_500_000,
                "measurement": "temporary MPI test fixture; no production allowance"}))
            cls.budget = DataBudget(baseline)
        cls.budget.require_output(cls.directory, planned_bytes=1_000_000)

    def run_native(self, text, ranks=2, backend="serial", threads=1):
        path = self.directory/(hashlib.sha256(text.encode()).hexdigest()+".tsv")
        with self.writer(path, planned_bytes=len(text.encode())+4096):
            path.write_text(text)
        env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE"}
        if hasattr(os, "geteuid") and os.geteuid() == 0:
            env.update(OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        command = [ARGS.mpiexec, *ARGS.launcher_arg, "-n", str(ranks), ARGS.executable,
                   "--jobs-file", str(path), "--backend", backend, "--threads", str(threads)]
        code, stdout, stderr, error = _bounded_process(command, directory=self.directory, environment=env, timeout=25.)
        self.assertIsNone(error, (error, stderr.decode(errors="replace")))
        report = json.loads(stdout)
        self.assertEqual(report["schema"], "bdmc.mpi-ensemble.v1")
        return code, report

    def jobs(self):
        return [MPIJob(run_id=f"job-{index}", condition_id=f"condition-{index//2}",
                       model=model, seed=713+index, parameter=.7+index*.1, count=8,
                       steps=5, warmup=2, diffusion=.5,
                       box_length=2. if model == "repulsive" else 4.,
                       shear_rate=0 if model == "repulsive" else (-.75 if index % 2 else .75))
                for index, model in enumerate(("repulsive", "hookean", "fene", "hookean", "fene"))]

    def test_complete_results_are_identical_for_one_two_and_three_ranks(self):
        text = serialize_jobs(self.jobs())
        reference = None
        for ranks in (1, 2, 3):
            code, report = self.run_native(text, ranks=ranks)
            self.assertEqual(code, 0, report)
            self.assertEqual(report["status"], "success")
            self.assertEqual(report["execution"]["ranks"], ranks)
            self.assertEqual([row["run_id"] for row in report["results"]], [job.run_id for job in self.jobs()])
            if reference is None:
                reference = report["results"]
            else:
                self.assertEqual(report["results"], reference)

    def test_rank_count_exceeding_job_count_keeps_order(self):
        text = serialize_jobs(self.jobs()[:1])
        code, report = self.run_native(text, ranks=3)
        self.assertEqual(code, 0, report)
        self.assertEqual(len(report["results"]), 1)

    def test_mixed_models_openmp_parity_and_rank_count_invariance(self):
        if not ARGS.check_openmp:
            self.skipTest("OpenMP was not explicitly enabled for this native harness")
        text = serialize_jobs(self.jobs())
        code, serial = self.run_native(text, ranks=1)
        self.assertEqual(code, 0, serial)
        self.assertGreater(serial["results"][0]["summary"]["potential_energy_per_particle"], 0)
        code, one = self.run_native(text, ranks=1, backend="openmp", threads=2)
        self.assertEqual(code, 0, one)
        code, two = self.run_native(text, ranks=2, backend="openmp", threads=2)
        self.assertEqual(code, 0, two)
        self.assertEqual(one["results"], two["results"])
        self.assertEqual(two["results"][0]["force_search"], "brute_force")
        for actual, expected in zip(two["results"], serial["results"]):
            self.assertEqual(actual["initial_state_fnv1a64"], expected["initial_state_fnv1a64"])
            self.assertEqual(actual["spec"], expected["spec"])
            for name, value in actual["summary"].items():
                expected_value = expected["summary"][name]
                if isinstance(value, list):
                    for observed, reference in zip(value, expected_value):
                        self.assertAlmostEqual(observed, reference, delta=1e-12)
                else:
                    self.assertAlmostEqual(value, expected_value, delta=1e-12)

    def test_malformed_and_duplicate_input_fail_collectively(self):
        valid = serialize_jobs(self.jobs())
        malformed = ["unknown-schema\n", valid+"\n", valid.replace("713\t", "-1\t", 1),
                     valid.replace("job-1\t", "job-0\t", 1), valid.replace("714\t", "713\t", 1),
                     valid.replace("\t0.001\t", "\tnan\t", 1)]
        for text in malformed:
            code, report = self.run_native(text)
            self.assertNotEqual(code, 0)
            self.assertEqual(report["status"], "failed")
            self.assertEqual(report["phase"], "input")
            self.assertIn("error", report)

    def test_one_rank_physical_failure_does_not_hang_or_reseed(self):
        jobs = self.jobs()
        # Index one runs only on rank one. Tiny q0 provokes a real FENE Euler
        # domain error after valid zero initialization, without a test-only flag.
        jobs[1] = replace(jobs[1], model="fene", q0=1e-12)
        code, report = self.run_native(serialize_jobs(jobs), ranks=2)
        self.assertNotEqual(code, 0)
        self.assertEqual(report["status"], "failed")
        failed = report["results"][1]
        self.assertEqual((failed["run_id"], failed["seed"]), (jobs[1].run_id, jobs[1].seed))
        self.assertEqual(failed["status"], "failed")
        self.assertIn("FENE", failed["error"])
        self.assertTrue(all(row["status"] == "success" for i, row in enumerate(report["results"]) if i != 1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True, type=lambda value: str(Path(value).resolve(strict=True)))
    parser.add_argument("--mpiexec", default=shutil.which("mpiexec") or "mpiexec")
    parser.add_argument("--launcher-arg", action="append", default=[])
    parser.add_argument("--scratch-root")
    parser.add_argument("--budget", type=Path, help="Existing shared budget; default is a tiny disposable test fixture")
    parser.add_argument("--check-openmp", action="store_true")
    ARGS, remaining = parser.parse_known_args()
    MPIJob, _bounded_process, serialize_jobs, DataBudget = _import_source_helpers()
    unittest.main(argv=[sys.argv[0], *remaining])
