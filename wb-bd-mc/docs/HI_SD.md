# Bounded imported HI references

[Attribution and license scope](../NOTICE.md).

The self-contained source snapshot under `third_party/hi_native` supports dense free-space Oseen/RPY Brownian references (N<=128) and deterministic athermal Stokesian references (N<=64). The bridge uses GSL; the pure default/reference components keep separate link boundaries. The snapshot has no runtime dependency on another development worktree.

The SD implementation ports selected JFSD invariant-tensor formulas and the nonzero resistance-table portion. See the unchanged [Apache notice](../third_party/hi_native/third_party/hi_sd/NOTICE), [license](../third_party/hi_native/third_party/hi_sd/LICENSE) and [provenance](../third_party/hi_native/third_party/hi_sd/SOURCE.json). A successful bounded free-space check does not establish periodic HI, thermal SD, complete suspension stress, network/HI coupling or converged transport/viscosity. The original native import hashes are retained in [IMPORT_MANIFEST](../third_party/hi_native/IMPORT_MANIFEST.json).

## Implementation status

RPY and SD are additional development with scientific validation pending. The
native C++ SD reference ports selected JFSD formulas and resistance data and
computes deterministic instantaneous free-space responses for up to 64
particles. A separate reduced integrator uses the resulting mobility in a
two-particle thermal benchmark. Periodic thermal RPY/SD, production transport
and native Apple Metal performance require separate evidence.

## Equilibrium pair benchmark and transport observables

[The equilibrium HI benchmark](HI_EQUILIBRIUM.org) compares pure BD
(\(K_\mathrm{sed}=1\)), thermal RPY and thermal evolution using the
deterministic SD pair mobility. It gives the diffusivity, finite
configurational viscosity and pressure definitions, reproducible values,
and the reasoning that limits the comparison to free-space pairs.
