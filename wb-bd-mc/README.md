# BD–MC modernization template (`AI_modulated_bd_mc`)

[Modernization scope](docs/AI_UPDATE_SCOPE.md) / [Validation](docs/VALIDATION.md)

[Modernization levels](docs/MODERNIZATION_LEVELS.org) / [Component guide](docs/CODE_GUIDE.org)

[Equilibrium hard-sphere MC and EOS comparison](docs/EQUILIBRIUM_EOS.org)

[Equilibrium RPY and Stokesian pair transport](docs/HI_EQUILIBRIUM.org)

[Self and collective diffusion](docs/DIFFUSION.org)

The equilibrium path is a bounded Level 3 example. A standalone C11 sampler
exports static hard-sphere observations for an independent EOS analysis.
The guide covers the soft-cosine BD comparator, physical time and transport
limits, a detached smoke run, and the evidence for scientific reuse.

## Origin and status

This is an independently copyable child of **wb-templates**, the reusable source-template layer of Agentic Workbench (wb). The compatible distribution name `AI_modulated_bd_mc`, Python import and CLI `bdmc` are preserved. It is a worked transformation of existing research software, not a claim that this model suits every research task.

**Gunwoo Park** is the original software author and directs this modernization. The original README uses the spelling **Gun Woo Park**. The upstream project is explicitly:

**<https://github.com/gwpark-git/dynamics_of_networks_and_colloids.git>**

This standalone source package derives from the `colloidal_dynamics` branch at `105fd3c96cd289b2f502ace04a332a2555e3307e`. The researcher directed the modernization.

**Validation is incomplete. Basic build, unit, physical-identity and selected regression checks have been performed, but exhaustive unit tests and regression tests have not been completed for every simulation mode, parameter range, platform or complete production simulation. This package cannot fully replace the original source code.** If a result or validation check is anomalous, preserve the reproducer and consult the corresponding original source, its recording conventions and the physical equations before relying on the result. See [validation and fallback](docs/VALIDATION.md).

## Modernization and validation

- **Level 1:** the portable legacy simulator passes selected frozen current-port
  regressions, including simple and step shear.
- **Level 2:** tests compare selected topology, transition and kinetic-clock
  operations with legacy laws. The complete active-shear workflow remains open.
- **Level 3:** a standalone hard-sphere MC sampler exports equilibrium data for
  independent pressure and structure analysis. Its sweep count is a trial clock.

The [modernization levels guide](docs/MODERNIZATION_LEVELS.org) links each level
to its test procedure. Two related areas still need scientific validation.

- **Additional development:** RPY and Stokesian dynamics (SD) are under
  development. Scientific validation is pending. The native C++ SD reference
  computes bounded deterministic free-space responses using selected formulas
  and tables ported from JFSD. See [HI implementation and evidence](docs/HI_SD.md).
- **Post-processing:** scientific validation of the analysis workflows is
  pending. Unit and reconstruction checks cover specific estimators and cases.
  See [observable definitions](docs/POSTPROCESSING.md).

## What is included

- Modern C++17 components and direct Python/NumPy bindings for Brownian and dumbbell references; repulsive, Hookean/FENE and imposed-flow reference models.
- The modified legacy transient-network simulator; labeled-topology and transition primitives; sparse neighbor search and guarded periodic/shear geometry.
- Checkpoint/restart, lineage-aware independent ensembles and optional ridge, Gaussian-process and PyTorch learning workflows.
- OpenMP/MPI build paths and provided CUDA/Apple Metal/Accelerate implementations, with explicit platform limits.
- Bounded free-space HI references and physically defined RDF, correlation and post-processing tools.
- Source tests and small frozen current-port fixtures, with upstream attribution and retained license notices.

[The update inventory](docs/AI_UPDATE_SCOPE.md) distinguishes implemented components from unfinished controllers and unvalidated backends. The distribution name is `AI_modulated_bd_mc`; the compatible Python import and CLI remain `bdmc`.

## Install and run

Use Python 3.10+, a C++17 compiler, CMake 3.18+ and development headers. The default HI build also requires GSL. Legacy builds require GSL and OpenMP. For example, on Debian/Ubuntu:

```sh
sudo apt-get install build-essential cmake pkg-config libgsl-dev python3-dev python3-venv
python3 -m venv .venv
. .venv/bin/activate
python -m pip install '.[test,reference]'
python -m bdmc --help
python -m bdmc run --spec examples/quickstart.json --steps 20 --stride 5 --output runs/quickstart
python -m pytest tests/python/test_postprocessing.py -q
```

Run from this package's root. Output directories must be fresh. The example is a tiny fixed-seed Brownian smoke run, not a validated production setting. See [build options, native tests and data budgets](docs/BUILD_AND_USAGE.md) for legacy execution, optional dependencies, backend selection and larger unattended jobs.

## Attribution and license

The original software and subsequent project-owned material retain the [MIT license](LICENSE), including the original copyright notice. [Third-party notices](NOTICE.md) identify Apache-2.0-derived HI material and external dependency terms. Cite the original project and identify this derivative's exact edition using the [citation metadata](CITATION.cff). Scientific or numerical errors may remain. The license's warranty disclaimer applies subject to applicable law.

## Structure and extension

The [component guide](docs/CODE_GUIDE.org) explains owned physical state, replaceable components, Python buffer ownership, exact topology, execution backends, learning and observer boundaries. The guides use Org, the plain-text format of Emacs org-mode. Markdown `#` and `##` headings correspond to Org `*` and `**`. Read the guide directly in any text editor before extending a component. The source, tests and examples run independently of wb-plugins.
