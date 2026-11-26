# mBLA modernization

[Open the interactive profile explorer](../modern-lv2/tools/profile-explorer.html) to view
the original prescribed-velocity inlet regression at CASE2 (`b = 1`, nominal TMP 1.15 kPa).
The regression view uses the original hollow-fibre selector, 200 stored
axial stations, and 400 samples in the distinct outlet section. It shows
archive, Level 1 and Level 2 values and their signed replay residual.
[Profile data format and use](../modern-lv2/docs/PROFILE_EXPLORER.md#load-the-recorded-cylinder-regression) /
[Current pressure-inlet validation](INLET_VELOCITY.org)

[Source package](../README.md) /
[Level 1 baseline](../modern-lv1/README.md) /
[Level 2 solver](../modern-lv2/README.md) /
[Level 3 provider example and compatibility test](../modern-lv3/README.md)

mBLA calculates concentration polarization, permeate velocity and cake
resistance in cross-flow ultrafiltration. The reference is the C++ solver used
by Park, Dhont and Nägele in [J. Chem. Phys. 164, 244119 (2026)](https://doi.org/10.1063/5.0336005).

| Scope | Recorded evidence |
|---|---|
| Level 1 Modernization | Modernized build and interfaces. The three frozen prescribed-pressure cases use second-order pressure FDM and first-order inlet axial velocity reconstruction. The current second-order inlet axial velocity reconstruction has independent verification. |
| Level 2 Modernization | Retained the fixed-point numerical architecture. Eight selected velocity-inlet CASE1/CASE2 TMP runs on the 200 × 400 grid matched Level 1 files byte-for-byte and archived mean TMP and permeate velocity at stored precision. Profiles were compared over the full axial domain of the admitted cases. CASE2 used the reference b = 1. After the inlet correction, two selected velocity-inlet cases were rerun against the archived replays. |
| Python and Numba | Ported the inner cake-resistance search to Python with optional Numba execution. The two paths matched C++ on 42 saved searches covering 9,712 trials. The outer solver and CLI are C++. |
| Level 3 Modernization | The C++ solver accepts an `EquilibriumTransport` provider for paired scalar `D_c/D_0` and `eta/eta_s` closures. A recorded native BD–MC/mBLA benchmark exercised live and pretabulated providers for a 5 kPa, 10 nm non-cake case. The live and cached linear paths matched all 250 axial and 1,440 radial printed values. |

The Level 3 package builds against the interface in the Level 2 C++ source
tree. The earlier, pre-correction Level 3 benchmark
applies to local equilibrium at low particle Péclet number, uses free-draining
Brownian dynamics and retains the built-in osmotic-pressure law. The benchmark
checks transport-provider composition for that non-cake case. Mesh convergence,
high-TMP inlet applicability and full hydrodynamic suspension viscosity
require separate scientific evidence.

In that internal benchmark, native BD–MC sampled transport inside the mBLA
process. The live run took 338.000 s, including 337.583 s for particle
sampling. The cached linear solve took 0.396 s with a precomputed table,
whose generation took 473.737 s separately. For repeated continuum runs,
the practical path is to interpolate a validated table within its measured
domain. Current development targets live sampling at selected states to
extend coverage. A future ML provider could implement the same interface
once it meets held-out property, continuum accuracy, domain and cost checks.
[Level 3 design and validation plan](LEVEL_3_DESIGN.org)

In a separate 25 × 120 regression case at 5 kPa and 10 nm, the built-in
unity law and an injected unity provider give byte-identical axial,
integral and radial outputs with the current default. The earlier provider audit showed that the solver can query concentrations outside
the physical range during intermediate numerical stages, even for a
non-cake case. A strict-range provider stopped at the first negative trial.
Tabulated providers therefore need an explicit numerical-trial policy and
a stated validity domain for physical results.

**Project Director** Gunwoo Park.
