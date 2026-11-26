# mBLA C++

[Open mBLA interactive validation](modern-lv2/tools/profile-explorer.html)

Select a recorded case to compare the archive, Level 1 and Level 2
at stored axial stations and across the outlet radius. The separate
**Site hollow-fibre reference** tab shows saved concentration and velocity
profiles through a selected point.

The independent [Level 1](modern-lv1/README.md) and
[Level 2](modern-lv2/README.md) packages each have CMake, Make, local CLI
installation and tests. Level 2 refactors the fixed-point controller, wall
updates and cake search within the existing numerical architecture.
The current C++ case driver declares the Level 3 local-equilibrium transport
interface. [Case driver header](modern-lv2/src/case_driver.h)
The [Level 3 package](modern-lv3/README.md) builds a provider example on that
driver and compares it with the Level 2 result. The [Level 3 design record](docs/LEVEL_3_DESIGN.org) describes the wider
scientific-composability workflow.
[Modernization overview](docs/MODERNIZATION.md)

## Modernization levels

| Level | Implemented and checked scope |
|---|---|
| Level 1 — portable baseline | Build and interface modernization supports replay with the historical first-order inlet axial velocity reconstruction and second-order pressure FDM. See the [Level 1 validation guide](modern-lv1/docs/VALIDATION.org). |
| Level 2 — numerical organization | The fixed-point architecture and output contract are retained. The current C++ solver and the Python/Numba inner cake-resistance kernel have separate replay evidence below. |
| Level 3 — scientific composition | `EquilibriumTransport` supplies scalar `D_c/D_0` and `eta/eta_s` to paired `LIVE` closures. The recorded native BD–MC/mBLA benchmark ran a 5 kPa, 10 nm non-cake case with on-demand and pretabulated providers. The live and cached linear paths matched all 250 axial and 1,440 radial printed values in that case. |

`modern-lv3` uses the transport interface in the Level 2 source tree. The pre-correction composition
benchmark applies to a local-equilibrium, low particle Péclet number case
with free-draining particle dynamics and the built-in osmotic-pressure law.
Full suspension viscosity and cake-regime applicability require their own
scientific evidence. [Runtime transport contract](modern-lv2/docs/SCIENTIFIC_GUIDE.org)

## Evidence from the JCP 2026 reference

The [regression data guide](modern-lv2/docs/PROFILE_EXPLORER.md) explains the
stored grids, quantities and case coverage. Open the [profile explorer](modern-lv2/tools/profile-explorer.html)
to compare the full axial and outlet radial sections side by side. The site
reference tab shows axial and radial concentration and flow through a selected
point in the independent hollow-fibre preview.

The reference is the C++ solver used in Park, Dhont and Nägele,
[J. Chem. Phys. 164, 244119 (2026)](https://doi.org/10.1063/5.0336005).
Level 1 modernizes the build and interfaces. Level 2 preserves the numerical
architecture and ports the **inner cake-resistance search** from C++ to
Python, with optional Numba acceleration. The outer solver and CLI are C++.
Eight selected CASE1/CASE2 TMP runs with prescribed inlet velocity at 200 × 400 reproduce Level 1 output
byte-for-byte and match archived mean TMP and permeate velocity at stored
precision. The accepted profiles span the full axial domain, with CASE2 run at
the reference setting b = 1. Mesh convergence and the high-TMP inlet regime
have their own applicability limits. The [Python/Numba replay](modern-lv2/docs/PYTHON_CAKE.org)
agrees with C++ on 42 saved searches and 9,712 trials. These are bounded
modernization checks for the recorded conditions and routines.

## Pressure-inlet velocity reconstruction

The pressure equation uses a second-order FDM. With prescribed inlet
pressure, the solver reconstructs inlet axial velocity from the solved
pressure with a three-point second-order forward difference. Set
`inlet_velocity_order=1` to replay the historical first-order inlet axial
velocity reconstruction. See the measured refinement,
C++–Python comparison and publication scope.
[Inlet velocity validation](docs/INLET_VELOCITY.org)

## Build the root C++ baseline

Run these commands from the mBLA package root to build its baseline CLI.
For a modernization-level build, use the [Level 1 guide](modern-lv1/README.md),
[Level 2 guide](modern-lv2/README.md), or [Level 3 guide](modern-lv3/README.md).
The Level 3 provider example builds against the sibling Level 2 source.

A C++17 template for Gunwoo Park's cross-flow ultrafiltration solver, preserving
its hollow-fibre concentration-polarization and cake-resistance calculations.
It provides a CMake library/CLI, owned matrix storage, validated configuration
and tables, and small immutable reference cases.
See [license scope](LICENSE), [notices](NOTICE) and
[third-party evidence](THIRD_PARTY_NOTICES.org).

Install or otherwise provide a C++17 compiler, CMake 3.16+, GSL 2+ and Python 3
for the standard-library test runner. No network, plugin, original dataset or
Python scientific stack is needed for the C++ build and standard tests.

```sh
cmake -S . -B /tmp/mbla-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/mbla-build -j2
ctest --test-dir /tmp/mbla-build --output-on-failure
mkdir /tmp/mbla-run
cp examples/pressure-inlet.inp /tmp/mbla-run/case.inp
cd /tmp/mbla-run
/tmp/mbla-build/mbla case.inp
```

`fout=result.dat` is resolved against the run directory and overwritten. A
completed run prints `MBLA_STATUS=converged`; this names the inherited stopping
criterion. Inspect conservation before interpreting its predictions. Run
`tools/inspect_output.py` from the copied template with the resulting table:

```sh
python3 tools/inspect_output.py /tmp/mbla-run/result.dat \
  --radius 0.0005 --permeability 6.7e-10 --csv /tmp/mbla-run/result.csv
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
