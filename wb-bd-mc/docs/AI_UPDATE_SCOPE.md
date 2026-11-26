# Scope of the modernization

Human project attribution is centralized in the [README](../README.md).

This is a derivative of the [original software](https://github.com/gwpark-git/dynamics_of_networks_and_colloids.git). The table describes the available components. [Modernization levels](MODERNIZATION_LEVELS.org) and [validation](VALIDATION.md) document the tested scope.

| AI update area | Included work | Boundary |
|---|---|---|
| Portable legacy builds | GSL/OpenMP build cleanup, C++17 Make/CMake paths and explicit build outputs | Modernized legacy sources contain behavioral repairs; this is not an untouched upstream copy |
| Legacy correctness repairs | Parameter/finite-domain checks, periodic/sheared geometry, recording clocks, topology update safety and analysis drivers | Supported interacting shear is 3D xy; arbitrary settings and every legacy option are not validated |
| Composable core | Owned state, explicit force/mobility/noise/integrator contracts, stable IDs and direct pybind11/NumPy access | Default Brownian/permanent-dumbbell core is a reference; a topology callback alone supplies no physical kinetics |
| Physical references | Repulsive Brownian and Hookean/FENE dumbbell responses, simple/step/oscillatory drive | Finite benchmark cases do not establish all nonlinear rheology or parameter limits |
| Topology | Exact labeled endpoints, loops/multiple bridges, sparse aggregate adjacency, transactional edits and transition primitives | Full modern transient-network BD–MC controller and modern shear-aware bonded-network integration remain unfinished; use the legacy simulator where appropriate |
| Equilibrium hard spheres | Independent periodic NVT hard-core MC, virtual-compression and contact-shell pressure estimates, cross-engine EOS comparison | Equilibrium state and pressure only. MC sweeps provide no Brownian time or transport coefficients. |
| Search and execution | Sparse occupied-cell neighbor queries, explicit OpenMP execution and real MPI independent-run ensembles | MPI is independent replicas, not single-system domain decomposition |
| Accelerators | CUDA float64, Apple Metal float32 and Accelerate source/build options | Full network-controller acceleration and hardware-specific accuracy remain unverified |
| Restart and provenance | Default-core checkpoint/RNG cache, explicit seeds, observations, lineage and independent ensemble records | Modern seeds do not reproduce old GSL paths; restart does not cover every backend/model |
| Learning workflows | Dataset construction, lineage-grouped splits, ridge/GP/MLP, saved/reloaded numeric inference and bounded learning-budget demos | No adopted or validated learned kinetic controller, complete rheology surrogate or guaranteed generalization |
| Hydrodynamics | Frozen, self-contained imported Oseen/RPY BD and deterministic athermal SD references; GSL bridge and retained third-party notices | Dense free-space reference bounds: BD N<=128, SD N<=64; no complete periodic HI, thermal SD, full suspension stress or coupled network/HI solver |
| Spatial analysis | RDF, directional/planar pair density, count-normalized coordination and explicit correlation tools | Geometric, density and estimator domains must match the data |
| Post-processing | Synchronized ragged topology, chain/edge measures, exact images, linear FFT correlation, stress projection, conservative virial, positive histograms/ECDF, S(q) and explicit MSD/diffusion contracts | End-to-end scientific validation of the analysis workflow is pending |
| Verification | Analytical/invariant tests, synthetic cases, v1/v2 current-port fixtures and bounded replay controls | Frozen outputs detect deviations from the captured behavior. Use analytical and physical checks for independent validation |

The preserved valid choices include chain-weighted versus distinct-edge statistics, single-origin versus all-origin estimators and the historical `1/5`, `1/30` isotropic stress coefficients. Repairs target their mathematical definitions rather than assuming every older variant was wrong.

Use `tools/validate.py` for the general validation checks described in [Build and usage](BUILD_AND_USAGE.md).
