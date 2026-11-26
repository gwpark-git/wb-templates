# Profile explorer

Open `tools/profile-explorer.html` directly in a browser. The template
contains a synthetic example. Select **Open profile data** to view a local
JSON file. The page reads selected values in memory and renders them in Canvas.

- Show or hide each grid with its checkbox.
- Select an axial node with the slider, an axial chart, or one of the jump
  buttons. The table shows wall volume fraction, minimum radial volume
  fraction, and particle-flux deviation at that node.
- Select **Near wall** for the first 8% of the radius. With only the fine
  grids visible, the axial particle-flux chart automatically uses their
  smaller variation to set its scale.

## Data format

Set `schema` to `wb.mbla.profile-explorer.v1`. Use the common ordered axial
positions `z_over_L` starting at zero, with the domain length `length_mm`
in millimetres. For each series, provide its wall volume fraction, axial
particle flux in consistent units, radial coordinates as distance from the
wall divided by the radius, and one radial concentration array per axial
node. Each series may have its own radial grid. Peak indices are zero-based
and should come from the diagnostics stated with the data.

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

The page has no server, CDN, solver or package dependencies. Keep any
research run data separate from the reusable template. A browser plot
inspects the supplied cases under their own scientific scope.
