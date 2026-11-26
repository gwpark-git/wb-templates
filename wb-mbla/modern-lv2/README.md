# mBLA C++ — Level 2 Modernization

A C++17 template for Gunwoo Park's cross-flow ultrafiltration solver, preserving
its hollow-fibre concentration-polarization and cake-resistance calculations.
This standalone edition includes the C++ source, CMake and Make build entry
points, local CLI installation, unit comparisons and frozen regression cases.
Level 2 keeps the serial case driver, fixed-point update, numerical kernels
and output formats. It clarifies case execution, wall and cake updates, and
termination. The [operation map and trace procedure](docs/LEVEL_2.org) describe
the former and current implementation.
See [license scope](LICENSE), [notices](NOTICE) and
[third-party evidence](THIRD_PARTY_NOTICES.org).

## Python and Numba within Level 2

The [inner cake-resistance search](docs/PYTHON_CAKE.org) from the
[JCP 2026 C++ solver](https://doi.org/10.1063/5.0336005) now has a Python
implementation and an optional Numba execution path. The outer fixed-point
solver and CLI are C++. The Python/Numba replay matches C++ on 42 recorded
searches and 9,712 trials. The separate solver regression covers eight
CASE1/CASE2 TMP conditions at 200 × 400, including full axial profiles in
the admitted domain and the CASE2 reference setting b = 1.
See the [Level 2 evidence and limits](docs/LEVEL_2.org).

The pressure equation uses a second-order FDM. With prescribed inlet
pressure, the solver reconstructs inlet axial velocity with a second-order
forward difference by default. Frozen original-source outputs use
`inlet_velocity_order=1` to select the historical first-order reconstruction. The quickstart input uses the verified 49 × 240
grid at nominal TMP 5 kPa. [Validation and legacy replay](docs/VALIDATION.org)

From this folder, build, test and install into a new directory. The standard
workflow uses a C++17 compiler, CMake 3.16+, GSL 2+, Python 3 and GNU Make.

```sh
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/mbla-lv2.XXXXXX")"
make test install BUILD_DIR="$work_dir/build" PREFIX="$work_dir/install"
mkdir "$work_dir/run"
cp examples/pressure-inlet.inp "$work_dir/run/case.inp"
(cd "$work_dir/run" && "$work_dir/install/bin/mbla" case.inp)
```

`make install` installs the CLI and notices, defaulting to `BUILD_DIR/install`.
Use `make test` to run the standard checks before installation. For a CLI-only
build, CMake can be used directly with `-DBUILD_TESTING=OFF`. That route needs
the C++ compiler, CMake and GSL.

```sh
cmake -S . -B "$work_dir/cli-build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build "$work_dir/cli-build" --parallel 2
cmake --install "$work_dir/cli-build" --prefix "$work_dir/cli-install"
```

## What the tests compare

- `mbla_units` checks matrix ownership and bounds, compares numerical
  components against specified values or independent formulas, and checks
  fixed-point guard, threshold and stopping semantics.
- `mbla_regression` uses historical first-order inlet axial velocity reconstruction to compare the
solver's ten-column outputs for non-cake,
  mild-cake and demanding-cake cases with frozen original-source fixtures
  under the declared common floating-point policy. The documented output
  precision determines the regression allowance.
- `mbla_input_contract` and `mbla_output_contract` check ordered continuation,
  batch termination, input errors and late write failures in the main output
  and sidecars.

The [pre-correction package check receipt](provenance/PACKAGE_CHECKS.json) records the tested
toolchain, intermediate-state comparisons and checks through the installed
executable. The standalone tests use the included examples and frozen fixtures.

An internal application can link the `mbla_case_driver` target and inject
an equilibrium transport provider. It selects `LIVE` for both diffusion and
viscosity. See the [runtime transport contract](docs/SCIENTIFIC_GUIDE.org)
for reduced units, validation, and the standalone CLI behavior.

## Run the cake-resistance search from Python

The [Python kernel](python/mbla_cake.py) accepts the radial fields for one
resistance search. It provides a standard-library implementation and an
optional NumPy/Numba backend. The compiled path reuses the geometry and
permeation-independent integral and evaluates each exponential once per
radial node and trial. The [Python benchmark guide](docs/PYTHON_CAKE.org)
gives the data contract, usage and measured performance of this inner loop.

## Explore saved profiles interactively

Open the [offline profile explorer](tools/profile-explorer.html).
It starts at the recorded CASE2 comparison (`b = 1`, nominal TMP 1.15 kPa).
**Axial profile** and **Outlet radial profile** show 200 stored axial
stations and 400 outlet samples side by side. The radial probe gives `y/R`,
`r/R` and `y/δ_CP` for the selected case. Choose **Site hollow-fibre reference**
to inspect the saved 5 kPa Python example with concentration and velocity in both
directions. **Open profile data** reads a compatible local file.
[Data format and use](docs/PROFILE_EXPLORER.md)

## Inspect a result

`fout=result.dat` is resolved against the run directory and overwritten. A
completed run prints `MBLA_STATUS=converged`; this names the inherited stopping
criterion. Inspect conservation before interpreting its predictions. Run
`tools/inspect_output.py` from the copied template with the resulting table:

```sh
python3 tools/inspect_output.py "$work_dir/run/result.dat" \
  --radius 0.0005 --permeability 6.7e-10 --csv "$work_dir/run/result.csv"
```

The examples distinguish a non-cake case, a mild cake case, and a demanding
cake regression. The latter exposes substantial particle-flux error; the mild
case also lacks a successful refined-grid confirmation. These are useful
engineering and numerical diagnostics, **not grid-converged cake validation**.
The inherited output has mixed normalization in its cake column; read the
[scientific/I/O guide](docs/SCIENTIFIC_GUIDE.org) before using it.

[Fixed-point method](docs/FIXED_POINT_METHOD.org) explains the history-dependent
cake safety rule and the researcher-stated limit for non-monotonic wall profiles.
[Validation](docs/VALIDATION.org) describes the inlet axial velocity checks
and the historical first-order velocity reconstruction used by the frozen fixtures. The 1 kPa
comparison uses a separate semi-analytic Python solver for non-cake
conditions. The Level 2 Python/Numba implementation covers the inner
cake search.
