# Style and color references

## Researcher-authored style

Gunwoo Park authored the original Python plotting scripts and style used for:
G. W. Park, J. K. G. Dhont and G. Nägele, “Modeling crossflow filtration:
Effect of shear on particle-enriched polarization and cake layers,”
*The Journal of Chemical Physics* **164**, 244119 (2026),
[doi:10.1063/5.0336005](https://doi.org/10.1063/5.0336005).
The examples use synthetic data. Citation is distinct from software licensing.

Within this collection, this paper is the source reference for the `wb-style`
adapter. The independent `plot-styles` common personal preset also draws on
researcher-authored helper scripts; this citation does not assert that all of
its defaults or layouts reproduce that paper.

## Color source

The independent `plot-styles` color mode uses **batlow** from **Fabio Crameri**,
*Scientific colour maps*, version8.0 (2023),
[doi:10.5281/zenodo.8035877](https://doi.org/10.5281/zenodo.8035877).
It is accessed through **Callum Rollo's cmcrameri** Python package,
[project source](https://github.com/callumrollo/cmcrameri).
Version1.10's actual wheel metadata explicitly identifies this colour-map edition.
The monochrome `wb-style` example does not call a colormap; this reference explains
the color resource used by the collection, not a new dependency of wb-style.

The inspected cmcrameri1.10 license is MIT: colormaps copyright2020 Fabio Crameri;
Python scripts/packaging copyright2020 Callum Rollo. Preserve its copyright and
permission notice if copying/substantially redistributing the software/data.
The collection reads the independently installed batlow resource; it does not
bundle or replace that data. The retained actual license is in the color package's
THIRD_PARTY_NOTICES directory. Scientific colour maps citation and cmcrameri
packaging attribution serve different purposes. Other dependencies retain their
own terms; a palette license does not clear unrelated inherited style code.

## Named dash presets

The [official Matplotlib Linestyles gallery](https://matplotlib.org/devdocs/gallery/lines_bars_and_markers/linestyles.html)
is the reference for named dash tuples, including **dashdotdotted**.
The [Matplotlib v3.10.7 gallery source](https://github.com/matplotlib/matplotlib/blob/v3.10.7/galleries/examples/lines_bars_and_markers/linestyles.py)
was used for the comparison. The packages store on/off lengths in points as
JSON arrays and supply zero offset to Matplotlib. Existing local variants
retain their stored lengths. See the child technical reference for the exact
differences. The [Matplotlib license notice](THIRD_PARTY_NOTICES/Matplotlib-3.10.7-LICENSE.txt)
is included with the distribution.
