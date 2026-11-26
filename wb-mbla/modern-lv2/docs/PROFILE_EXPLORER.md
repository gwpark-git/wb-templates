# Profile explorer

Open the [Profile Explorer](../tools/profile-explorer.html) directly in a browser. It opens the embedded
2026-11-25 reference comparison at its recorded initial case, CASE2
(`b = 1`, nominal TMP 1.15 kPa). The original case selection and stored
values are unchanged. **Open profile data** accepts another saved cylinder
regression or a v1 profile file and reads the selected file in memory.

## Load the recorded cylinder regression

The offline HTML embeds the original `CUSTOM_DATA.json`
(`wb.mbla.cylinder-comparison.v3`) from *mBLA physical observable recheck
under limited cake coverage*. It contains eight full comparisons and
three coverage summaries. All eight full CASE1/CASE2 comparisons
prescribe inlet velocity. The separate SPHS 20 kPa coverage summary
prescribes inlet pressure and uses second-order pressure FDM with the
historical first-order inlet axial velocity reconstruction. The saved data also retain
the original Level 1/2 and inner-kernel evidence editions. The
[current inlet report](../../docs/INLET_VELOCITY.org) records new
pressure-inlet checks and the two prescribed-velocity spot replays. The
reference study gives the scientific context. [JCP 2026 article](https://doi.org/10.1063/5.0336005)

1. Select a case from the measured summary or the **Case** list. The paired
   CASE1 and CASE2 (`b = 1`) rows show stored mean TMP (kPa), along with archive
   and replay mean permeate velocities (µm/s). Nominal TMP identifies the run.
2. In **Axial profile**, the schematic, chart and slider select one of the case's
   **200 stored stations**. They include both endpoints, `z/L = 0` and
   `z/L = 1`. The position readout gives the ordinal, exact stored normalized
   coordinate and derived distance in millimetres. The table shows stored
   column values and signed Level 2 minus archive residuals. **Whole-axis
   regression differences** gives each changed observable's maximum absolute
   difference, changed-station count and first maximum position.
3. Beside the axial chart, **Outlet radial profile** selects one of the case's
   **400 stored samples** of the separate outlet section. Axial position
   stays at `z/L = 1`. The radial slider moves by sample ordinal, preserving
   each adaptive `y/R` coordinate, including the last stored sample
   `y/R = 0.999999` in these cases. **Wall-log** and **Linear** change only the
   display scale. The radial section comes from the recorded outlet sidecars
   before final axial analysis. The readout gives `y/R`, `r/R = 1 - y/R` and
   `y/δ_CP = (y/R)/(δ_CP/R)`. For each saved run, the dashed guide is
   `y/δ_CP = 1` using its fixed `δ_CP = D₀/(Lₚ ΔP_HP)`. The
   solver uses the input reference TMP for `ΔP_HP`; the displayed mean TMP
   is a separate observable. Imported comparisons have no saved CP scale.
4. The recorded receipts identify identical Level 1 and Level 2 axial tables
   and outlet sidecars. The reader shares those verified replay arrays and
   draws the archive in black, Level 1 in gray and Level 2 in dashed red.
   The extensive-cake cases give mean permeate velocity and cake coverage.
   The available Numba benchmark covers the saved inner-search trials.

The numerical comparison includes the inlet. Physical interpretation of
entrance development requires a separate model. The viewer opens offline
and performs no solver execution.

## Inspect the site hollow-fibre reference

Choose **Site hollow-fibre reference** to inspect the saved Python suspension
case at nominal TMP 5 kPa, particle radius 10 nm and bulk volume fraction 0.001.
Drag the cutaway or use the two sliders to select an axial position and wall
distance. The four plots show both directions through that point for
concentration and the chosen axial or transverse velocity. Cutaway colours
come from all 24 × 80 adaptive concentration cells. The profiles use the
24 × 320 saved display samples, evaluated in physical coordinates.

For this case the fixed site reference length is `δ_CP = 6.398119e-6 m`
with `δ_CP/R = 1.279624e-2`. The pressure reference comes from the saved
Hagen–Poiseuille boundary conditions. In both views `δ_CP` is a characteristic
length scale. The guide shows one such length from the membrane.

## Use a v1 profile

- Select or drag inside the **Hollow-fibre selector** to choose an axial
  station and a radial probe from the wall toward the axis. The sliders
  provide the same selections with the keyboard.
- Choose **Wall-log** or **Linear** for the schematic's radial coordinate.
  Both preserve the exact wall position at `y/R = 0`.
- Use **Quantity** to switch between wall concentration, TMP, permeate
  velocity and other available axial observables. A supplied profile can
  also include an axial velocity field across the radius.
- Show or hide series with their checkboxes. The table displays the
  selected axial value, wall and minimum concentration, flux drift, and
  the nearest stored radial sample with its own `y/R`.
- Select **Near wall** for the first 8% of the radius. With only the fine
  grids visible, the axial particle-flux chart automatically uses their
  smaller variation to set its scale.

### v1 profile format

Set `schema` to `wb.mbla.profile-explorer.v1`. Use the common ordered axial
positions `z_over_L` starting at zero, with the domain length `length_mm`
in millimetres. For each series, provide its wall volume fraction, axial
particle flux in consistent units, radial coordinates as distance from the
wall divided by the radius, and one radial concentration array per axial
node. Each series may have its own radial grid. Peak indices are zero-based
and should come from the diagnostics stated with the data.

Optional `axialQuantities` and `radialQuantities` name extra curves with
units. In each series, `axialValues[j]` is one value per axial node and
`radialValues[j][i]` is one radial row at axial node `i`. The outer
indices match the order of their quantity labels. The
`outletOnly: true` case uses one stored row at `z/L = 1` per radial field.
The original form with wall, particle flux and concentration is also accepted.

```json
{
  "schema": "wb.mbla.profile-explorer.v1",
  "title": "My comparison",
  "subtitle": "Synthetic input example",
  "length_mm": 50,
  "z_over_L": [0, 1],
  "peaks": {"maxFluxDrift": 1, "maxProfileDifference": 1},
  "series": [{
    "label": "Grid A",
    "wall": [0.01, 0.2],
    "particleFlux": [1, 1.02],
    "radialPosition": [0, 1],
    "radialProfiles": [[0.01, 0.01], [0.2, 0.01]]
  }]
}
```

The reader accepts up to four matched-axial-grid series, 400 axial nodes,
4,000 radial nodes per series, and a 32 MiB JSON file. Inlet flux must be
positive. The signed flux deviation is `100 * (Phi(z) / Phi(0) - 1)` for
each series. Negative concentration values remain visible as diagnostics.
The standard Level 1 output table supplies wall concentration and particle
flux. Full radial profiles require a compatible saved field export.

## Renderer source

The v3 reader adapts the original regression renderer, with its research-site
`FiltrationSchematic`, `LinePlot` and radial coordinate transformations. The
v1 profile view uses the same schematic with its existing Canvas charts.
The regression cutaway shades the stored wall and outlet concentration strips.
The site cutaway uses the site's field renderer and sampling routines for its
saved suspension case.

For maintenance from the full source checkout, run
`node tools/build_site_selector.cjs RESEARCH_SITE ESBUILD_MODULE` from
Level 2. Supply the research-site source checkout with its dependencies
and local esbuild. The builder reads the original `CUSTOM_DATA.json`
and `VIEWER_PROVENANCE.json` in `docs/validation-report/`. Readers of the
distributed explorer can open the bundled HTML directly.
The builder checks the regression payload against its recorded SHA-256 before
embedding the original bytes. It packages the renderer, the required fields
from the manifest-verified site preset, and the fixed scales derived from the
matched C++ inputs in `tools/validation-cp-scales.json`. That table records
SI values and the input digests. The HTML includes React dependency notices.

Project Director: Gunwoo Park
