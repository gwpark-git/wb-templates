# Modern native components

Start with the [package README](../README.md), [AI update scope](../docs/AI_UPDATE_SCOPE.md) and [build guide](../docs/BUILD_AND_USAGE.md). C++17 headers are under `include/bdmc`, implementations under `src`, and physical/API tests under `tests`.

`_core` owns the serial Brownian/permanent-dumbbell reference and default restart state. `_reference` supplies repulsive and driven Hookean/FENE models. `_network` supplies labeled topology/transition/bonded-force primitives, not a complete migrated transient-network controller. `_analysis` supplies structural/correlation kernels. Optional `_hydrodynamics` supplies the bounded imported free-space models. Seeded modern streams are distinct from legacy GSL streams. Backend source presence does not establish execution on every platform; [validation limits](../docs/VALIDATION.md) apply.

The optional [equilibrium sampler](../docs/EQUILIBRIUM_EOS.org) builds a
standalone hard-sphere MC executable under `equilibrium/`. Its CSV interface
supports pressure and structure analysis with a configurational trial clock.

The optional [transport benchmark](../docs/METAL_BENCHMARK.org) compares full fixed-volume property estimation across the serial and Metal execution paths.
