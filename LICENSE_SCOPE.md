# License scope

The root [MIT license](LICENSE) is the shared grant for project-owned material
identified below. Distinct copyright notices, scope statements and third-party
terms remain with the material they cover. A detached export or wheel carries
the license and notices applicable to its contents.

| Material | Grant and attribution | Exceptions and standalone copies |
|---|---|---|
| Collection guides, templates and `wb-python` documentation | Root MIT for authored material directed by Gunwoo Park | `wb-python` is a documentation entry. Its independently distributed runtime has its own package grant. |
| `wb-bd-mc` | [Package MIT](wb-bd-mc/LICENSE) preserves the original 2020 and current attribution | The combined package and its wheels carry this grant. [Notices](wb-bd-mc/NOTICE.md) cover retained Apache-2.0 HI, PCG32 and other sources. |
| `wb-mbla` | [Root MIT](LICENSE) for the new authored adaptation where the [mBLA scope](wb-mbla/LICENSE) applies it | Inherited C++ and fixtures have separate scope. Keep the [third-party notices](wb-mbla/THIRD_PARTY_NOTICES.org) and level-specific scope files. CMake installation and detached export carry the shared MIT text as `LICENSE.NEW`. |
| mBLA profile explorer | Authored viewer and retained renderer notices | React terms and the original reference data are described in the [profile explorer data guide](wb-mbla/modern-lv2/docs/PROFILE_EXPLORER.md). The data are outside the shared MIT grant. |
| `scientific-io` | [CSV MIT scope](scientific-io/LICENSE_CSV) | The wheel carries its two license files and [notice](scientific-io/NOTICE.md). |
| `plot-styles` and `wb-style` | Root MIT for project-owned presets, adapters and examples | The wheel build reads local MIT copies and retains the [plot-style notices](plot-styles/NOTICE.md) and [wb-style notices](wb-style/NOTICE.md). |
| Quarto documentation example | Root MIT for first-party source and synthetic data | The independently copyable kit carries its copy and [renderer and font notices](documentation/examples/quarto-comparison/THIRD_PARTY_NOTICES.md). |
| BD–MC equilibrium EOS example | Root MIT covers the new documentation | The retained CSV and figures have separate [research provenance](examples/bdmc-eos/README.md). |

GSL is separately installed GPL software linked by some C++ targets. This source
candidate does not bundle linked binaries or establish combined binary release
rights. NumPy, SciPy, Matplotlib, cycler and cmcrameri retain their actual terms.
Authentic batlow data come from installed cmcrameri, not vendored data or an AI
substitute. The documentation examples embed font subsets in rendered artifacts.
The cited paper's CC BY notice does not license separate source/data/model files.
