# mBLA Level 3 modernization

Level 3 supports **Scientific Composability for Multiscale Simulation**. Use the
[case-driver interface](../modern-lv2/src/case_driver.h) to supply reduced
collective diffusivity $D_c/D_0$ and relative viscosity $\eta/\eta_s$ to the
existing continuum solver. A case selects `LIVE` for both properties in its
input file. The osmotic-pressure law remains separately selectable.

The [C++ example](examples/unity_provider.cpp) injects a synthetic constant
provider and runs the existing solver through `mbla::run_cases`. The held
provider is valid for the lifetime of that call. It returns unity for every
finite numerical trial so the test covers negative Runge–Kutta trial
concentrations as well as physical states. The unity example checks the
interface's behavior. It supplies no dispersion-property estimate.

## Build and verify

Keep `modern-lv2` and `modern-lv3` together inside the mBLA template. Build
with C++17, CMake 3.16+, GSL 2+ and Python 3 for tests. From `modern-lv3`:

```sh
make test BUILD_DIR=/tmp/mbla-lv3-build
```

For a single run, copy the [example input](examples/transport-unity.inp) into a
fresh working directory. Run `/tmp/mbla-lv3-build/mbla-lv3-unity` with that input
path from the working directory. The example writes `result.dat`, `result.dat.Fval`
and `result.dat.Rout` there.

The [compatibility test](tests/compare_unity.py) compares all three tables
byte-for-byte with the Level 2 `CT` closures. The CMake project also runs the
Level 2 regression and compatibility tests. The solver, numerical settings and
fixture inputs are shared between the two provider paths.

## Scientific scope

The [Python/Numba and PyTorch worked guide](docs/PYTHON_NUMBA_ML.org) walks
through an inference batch, a versionable table and a lookup compatible with
Numba. Use it to plan the Python-first solver path and validate the property
boundary before using supervised inference in a physical case.

See the [Level 3 design and internal BD–MC/mBLA benchmark](../docs/LEVEL_3_DESIGN.org)
for the contract, benchmark conditions and applicability. A physical provider
must define its estimators, units, state domain and numerical-trial behavior.
Particles and offline interpolated tables can each supply properties through the
same interface, subject to that contract.

Original Level 3 example code, build and tests use the MIT grant described in
the [mBLA license scope](../LICENSE). The collection carries the shared text
in its root LICENSE. A detached mBLA export carries it as LICENSE.NEW.
The linked Level 2 solver and its dependencies have [their own terms](../modern-lv2/LICENSE).
