# Validation scope and reference to the original source

[Project attribution](../README.md)

## Scope

**Build checks, unit tests and selected regressions cover the cases named below. Full production simulations and scientific parameter ranges still need validation.**

Build checks establish that selected configurations can run. Frozen-trajectory comparisons check repeatability against a captured numerical baseline. Independent physical oracles and convergence studies test model and observable assumptions.

## Modernization status

The modernization claim covers selected checks of capabilities already in the
original BD–MC package. Full workflow equivalence remains open. RPY and SD are
additional development with scientific validation pending. Post-processing
also awaits end-to-end scientific validation. See [modernization levels](MODERNIZATION_LEVELS.org)
for the relationship between tested capabilities and each level.

## Three verification steps

| Step | Evidence for BD–MC | Current scope |
| --- | --- | --- |
| Unit verification | [Network unit tests](../modern/tests/test_network.cpp) check labeled topology, proposal and acceptance behavior, and the kinetic clock. The `BDMC_NETWORK_LEGACY_ORACLE=ON` build also runs selected legacy-law comparisons. | Selected modernized operations and their matched legacy laws. |
| Regression against a retained baseline | [The v2 replay harness](../tests/reference_suite.py) runs the legacy executable directly against the [repaired current-port fixtures](../tests/reference/v2/README.md), including both signs of simple and step shear. | Recorded behavior of the repaired legacy implementation. The full modern active-shear workflow still needs a matching run and event trace. |
| Independent physical validation | [Legacy physics checks](../tests/test_legacy_physics.cpp) use independent lattice enumeration and energy/force/virial identities. [Modern reference checks](../modern/tests/test_reference.cpp) compare free Brownian variance and Hookean finite-step moments with analytical expectations. | The named physical limits and observables. Full transient-network rheology and comparisons with published data remain separate validation tasks. |

The remaining Level 2 workflow check needs exact labeled-event ordering and tolerance-based numerical agreement. Invoke the legacy executable directly for active shear. The Python `LegacyRunner` accepts quiescent cases only.

## Recorded bounded evidence

The [equilibrium EOS comparison](EQUILIBRIUM_EOS.org) adds a separate hard-sphere MC unit target, a stored-observation replay and a matched LAMMPS soft-BD comparison against Carnahan–Starling. Pressure agreement has model-, density- and estimator-specific limits. The MC clock covers trials and provides no physical transport measurement.

The [Level 3 equilibrium contract](EQUILIBRIUM_EOS.org#level-3-acceptance) defines the producer–consumer boundary for equilibrium pressure and structure and provides a detached example.

The [post-processing tests](../tests/python/test_postprocessing.py) exercise explicit estimator definitions. When comparing reconstructed energy or stress against rounded data, `rtol=2e-4, atol=2e-6` describes the chosen compatibility allowance. Treat printer-error propagation and frozen-output regression as separate comparisons.

The [Metal transport benchmark](METAL_BENCHMARK.org) defines a paired runtime and numerical comparison. The Linux CPU reference and explicit Metal-unavailable status were exercised. Native Apple throughput and float32 transport accuracy still need device evidence.

The [equilibrium free-space pair benchmark](HI_EQUILIBRIUM.org) uses the
existing RPY/SD mobility and a reversible thermal pair integrator. It compares
matched equilibrium observables against a pure-BD unit-mobility reference.
It is a two-sphere verification case for mobility and stress relaxation.

## Unchecked or incomplete scope

- Exhaustive source, branch, model and parameter coverage; full production rheology, relaxation, diffusion and viscosity convergence; timestep/size/statistical uncertainty across all workloads.
- Cross-platform bitwise equivalence. Modern `std::mt19937_64` and legacy GSL/libc initialization are distinct; the same seed does not make their paths identical. Compiler/libm/thread scheduling and update order can affect results.
- Complete migrated modern transient-network control, complete network/HI coupling, domain-decomposed network execution and a validated learned kinetic controller.
- CUDA/Apple target build/execution/performance on the unavailable hardware/toolchains; complete periodic HI, thermal SD and total suspension-stress theory through the bounded adapters.

## If a check or result looks wrong

1. Preserve the smallest reproducer: input, seed(s), source hashes, compiler/library versions, thread counts, exact command and original output. Preserve failed results and the original fixture.
2. Check units, clock/recording phase, geometry/shear sign, sparse schema, weights, particle/chain identities and estimator normalization. Distinguish all-pair RDF from selected distances and random-force diagnostics from conservative stress.
3. Consult **<https://github.com/gwpark-git/dynamics_of_networks_and_colloids.git>**, beginning with the recorded branch/commit and the corresponding force, stochastic-transition or recorder function. Keep an untouched upstream checkout for comparison. Do not assume the latest upstream branch is the same historical or publication edition.
4. Compare the discrepancy with the governing equation or a separately constructed oracle. Original code is a primary reference for provenance and intended implementation, not an automatic guarantee that either implementation is physically correct.
5. Review a proposed correction against that definition; preserve alternate valid observables, retain the failure, rerun the relevant tests and record any deliberate behavioral change. Do not loosen a tolerance, overwrite a reference or discard failed seeds simply to obtain a pass.

Report conclusions at the checked level. A short successful run does not establish a converged scientific result. License terms and the AI disclosure do not replace the user's scientific validation.
