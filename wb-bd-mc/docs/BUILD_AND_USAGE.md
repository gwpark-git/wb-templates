# Build and usage

[Attribution and limits](../README.md).

## Python package

The repository/distribution is `AI_modulated_bd_mc`, version `0.1.0`; Python imports and the CLI remain `bdmc`. Work from the repository root. Install Python 3.10+, a C++17 compiler, CMake 3.18+, Python development headers and GSL development files for the default HI component. Pip obtains the declared scikit-build-core/pybind11 build dependencies and NumPy runtime dependency.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install '.[test,reference]'
python -m bdmc run --spec examples/quickstart.json --steps 20 --stride 5 --output runs/quickstart
```

The smoke example writes a fresh `run.json` and compressed trajectory with initial/post-step frames. Do not reuse an existing output directory. Optional extras are `reference` (SciPy), `ml` (scikit-learn), `torch` (PyTorch) and `test` (pytest). Ridge needs NumPy; GP requires scikit-learn and MLP requires PyTorch. To build without the bounded HI component:

```sh
python -m pip install . -Ccmake.define.BDMC_BUILD_HYDRODYNAMICS=OFF
```

HI-specific tests then require that optional native module and are not a suitable minimal-build test selection. Default builds use portable CPU kernels. Select OpenMP with `-Ccmake.define.BDMC_ENABLE_OPENMP=ON`; MPI additionally needs its development/runtime toolchain. CUDA and Apple options require their actual target tools/hardware; they do not silently fall back. `python tools/build_backend.py --probe` reports tool discovery, not successful execution.

For full-job CPU/Metal timing and accuracy checks, use the [fixed-volume transport benchmark](METAL_BENCHMARK.org). It also provides an explicit CPU-only mode.

## Native and legacy builds

```sh
cmake -S modern -B build/native -DBDMC_BUILD_PYTHON=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build/native --parallel 2
ctest --test-dir build/native --output-on-failure
make -j2
make test-unit
```

Make builds the modified legacy simulator and two explicit-column analysis tools under `build/<compiler-name>/`. Legacy Make/CMake modes require GSL and OpenMP. CMake legacy mode uses `-DBDMC_BUILD_LEGACY=ON`; it is not enabled by default. `short_test.inp` is the original example input, not a proof of a fully validated production setup. Prefer a fresh, small declared input for initial legacy checks; the Python `LegacyRunner` adapter enforces its documented 3D/quiescent/fresh-start limits. Inspect its supported parameters before using it.

For periodic hard-sphere equilibrium sampling, enable `BDMC_BUILD_HARD_SPHERE_MC`
in CMake. The [equilibrium guide](EQUILIBRIUM_EOS.org) gives the standalone
producer build, fixed-seed sample and detached NumPy consumer.

With `-DBDMC_ENABLE_MPI=ON`, CTest also launches actual short independent-run MPI comparisons. Its Python launcher pins this checkout's pure Python helpers, even if another checkout has an editable installation. CTest scratch stays in the build directory and uses a disposable 2 MB test-only allowance (1.5 MB stop point); it grants no production data budget. Direct invocation of `modern/tests/test_mpi_ensemble.py` can instead select an existing `--budget` and a registered `--scratch-root`. Importing this CLI harness during pytest collection does not change package selection. MPI launch configuration remains platform-specific; these cases make no speedup claim.

## Data budget and optional broader checks

Larger validation tools require a fixed local storage budget. Initialize the baseline once per checkout before generating run or build data:

```sh
python tools/check_storage.py --baseline DATA_BUDGET.json --initialize-baseline runs build
python tools/check_storage.py --baseline DATA_BUDGET.json
python tools/validate.py --profile unit --compiler g++ --jobs 2
python tools/validate.py --profile regression --compiler g++ --jobs 2
```

The tool retains its 10 GB additional-data policy and 9 GB stop threshold; initialization does not grant another allowance inside an already budgeted shared workspace. Existing baseline files are never silently replaced. Do not reset the baseline between runs. Budget initialization is operator setup, not a bundled automatic installation hook.

`unit`, `regression`, `full` and `sanitizer` are available profiles. `verify` in `tests/reference_suite.py` reads existing fixtures; `create` is an explicitly separate baseline-capture operation and must not be used to make a failing verification pass. Default v2 preserves deliberate repairs; v1 retains the earlier behavior. The historical workspace-specific `run_repair_acceptance.py` batch entry is intentionally inactive in this standalone distribution.

## Direct API and optional learning examples

```python
from bdmc import Simulation, SimulationSpec
spec = SimulationSpec(n_particles=8, dimension=3, box_lengths=[8, 8, 8],
                      dt=0.001, seed=100)
sim = Simulation(spec)
frames = sim.sample(20, stride=5, include_initial=True)
```

The exact constructor contract is exposed by `SimulationSpec`; the checked JSON smoke input is [quickstart.json](../examples/quickstart.json). Optional learning examples use explicit fresh outputs, such as `python examples/train_calibration.py --model brownian --without-torch --output runs/calibration`. They can generate multiple simulations; inspect their settings and use a justified budget before running them. They demonstrate data/training/inference plumbing, not a validated scientific learning target.

Use [post-processing definitions](POSTPROCESSING.md) before interpreting saved observables, and [validation/fallback](VALIDATION.md) before treating a numerical result as scientific evidence.
