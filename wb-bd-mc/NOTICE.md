# Attribution and license scope

The standard license texts and retained third-party notices are preserved.

Gunwoo Park is the original project author and directs this modernization. The original public project is <https://github.com/gwpark-git/dynamics_of_networks_and_colloids.git>. The upstream MIT copyright line `Copyright (c) 2020 gunwoo` is retained verbatim, with a separate 2026 attribution line in [LICENSE](LICENSE). The standard MIT grant is unchanged. Its scope includes project-owned code, tests, examples, documentation and packaging material, subject to explicit third-party exceptions below. Copyrightability, third-party rights and recorded human review require their own evidence.

## Included third-party material

Selected deterministic Stokesian-dynamics formulas and a resistance table derive from **JFSD by William Torre and contributors**, repository <https://github.com/torrewk/Python-Jax-Fast-Stokesian-Dynamics>, commit `5234c959c0fdb0ba7ba32d6015d0478014a512c7`.

Retained terms and provenance:

- [Apache-2.0 license](third_party/hi_native/third_party/hi_sd/LICENSE).
- [Original imported notice](third_party/hi_native/third_party/hi_sd/NOTICE).
- [Source/formula provenance](third_party/hi_native/third_party/hi_sd/SOURCE.json).
- [Frozen import hashes](third_party/hi_native/IMPORT_MANIFEST.json).

The derived `third_party/hi_native/lib/hi_sd.cpp` implementation and `third_party/hi_native/third_party/hi_sd/resistance_table.h` retain their Apache-2.0 provenance. See the [HI notes](docs/HI_SD.md) for their scientific scope.

The optional standalone hard-sphere MC target compiles the retained unmodified
PCG32 C source and header by Melissa O'Neill, from
[pcg-c-basic](https://github.com/imneme/pcg-c-basic).
The [Apache-2.0 license](modern/equilibrium/pcg32/PCG_LICENSE.txt) and per-file
copyright notices travel with those files. The EOS analysis reads LAMMPS
output as data. LAMMPS retains its own license and notices.

## External dependencies

Dependencies are obtained separately and retain their own terms. In particular, GSL is GPL-licensed and is linked by the legacy and default HI components; the root MIT license does not relicense GSL or establish that every combined binary may be redistributed under MIT alone. Review the actual GSL version and distribution obligations when distributing linked binaries. This folder distributes source, not bundled GSL binaries. NumPy, SciPy, pybind11, scikit-build-core, scikit-learn, PyTorch, compiler/runtime, MPI and accelerator SDK terms likewise remain their owners' terms. No exhaustive dependency-license or professional legal certification is claimed.

## Scientific disclaimer

The modernization and its documentation may contain scientific, numerical or software errors. Basic/selected tests do not cover all simulations or establish complete equivalence to the original source. Preserve a reproducer and compare with the original source and physical definition when validation is anomalous. See [validation](docs/VALIDATION.md). These statements describe validation status and add no use restrictions to the retained licenses; their warranty/limitation clauses apply subject to mandatory applicable law.
