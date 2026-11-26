## What wb-templates is {#what-wb-templates-is}

wb-templates is the reusable scientific-computing layer of Agentic Workbench
(wb). It brings together code, data readers,
plotting styles and document examples that a researcher can copy into a new
project. The aim is to make it possible to reuse an old calculation or a familiar
figure routine while keeping the assumptions, inputs and checks needed to
understand the result.

Park describes this problem in his [Agentic Workbench project letter](https://github.com/gwpark-git/agentic-workbench-knowledge).
His research codes were tested for particular questions, and he reused his own
Python analysis and visualization routines over many years. When a new question
comes up, the useful parts can outlast the context in which they were written.
The templates show how to recover that context, retain an original comparison
and adapt the parts that the new project needs. The receiving researcher chooses
the scientific model, applicability domain and criteria for using the result.

Within wb, **wb-plugins** supplies reusable procedures for people and agents.
**wb-templates** supplies examples and reusable source. **wb-applications** can
assemble selected parts into an environment for a particular task. The templates
can also be copied and used independently with their declared dependencies.

## Scientific Composability for Multiscale Simulation {#scientific-composability}

A particle simulation and a continuum calculation can share a physical quantity
if both sides agree on its definition, units, applicable conditions and numerical
handling. For example, a continuum filtration solver can query a table of
concentration-dependent transport properties derived from a particle-scale
calculation. The solver needs to know what each ratio means, which concentrations
have physical support and how the table handles queries between or beyond its
recorded points. A composition check observes the two components exchanging and
using the quantities under stated conditions.

The collection documents three levels of modernization. Choose the level that
matches the change you need and the evidence you can collect.

| Level | Work and evidence |
|-------|-------------------|
| 1. Preserve a baseline | Retain a runnable calculation, record its inputs and interfaces, and pin outputs for comparison. |
| 2. Modernize the implementation | Preserve the numerical architecture and check changed operations and selected results against retained references. |
| 3. Compose scientific components | Work toward scientific composability for multiscale simulation. Define the exchanged quantity and run a producer–consumer check under specified physical and numerical conditions. |

The [BD–MC modernization guide](../../wb-bd-mc/docs/MODERNIZATION_LEVELS.org) and
[mBLA modernization overview](../../wb-mbla/docs/MODERNIZATION.md) describe the evidence
available at each level for their respective cases.

### A current Level 3 example

An internal wb-applications project couples the native BD–MC particle model
to the mBLA continuum solver. The live path calls BD–MC to evaluate transport
properties for the filtration calculation. The tested live evaluations are
computationally expensive, so the practical route uses interpolation of an
offline table of calculated properties over its stated range. The application
is internal. Development is extending live evaluation to conditions where
the solver needs a fresh property calculation. The same interface could
support a future ML surrogate after its training data, validity range and
predictions are independently validated.

## Choose a template {#choose-a-template}

| If you want to… | Start here | What you can reuse |
|-----------------|------------|--------------------|
| Make related figures in color and monochrome | [Plot styles](../../plot-styles/README.md) | A Matplotlib style with shared geometry, line weights and typography, plus synthetic layout and legend examples. |
| Reuse the appearance of Park's publication plotting scripts | [wb-style](../../wb-style/README.md) | A standalone `wb_style` preset for ticks, spines, typography, markers and dash patterns, with a neutral CSV example. |
| Read a measurement table with explicit column choices | [Scientific CSV I/O](../../scientific-io/README.md) | Header inspection and numeric validation for selected CSV columns. The caller supplies units and physical interpretation. |
| Share a detailed worked analysis and a short presentation | [Quarto comparison template](../examples/quarto-comparison/README.md) | An independently copyable example with a continuous report, linked Reveal.js slides, a synthetic comparison, a source diff and local build instructions. |
| Recover or adapt a Brownian dynamics or Monte Carlo research code | [BD–MC](../../wb-bd-mc/README.md) | A portable legacy baseline, C++ and Python components, and a standalone equilibrium pressure and structure example. |
| Reproduce or extend a cross-flow filtration calculation | [mBLA](../../wb-mbla/README.md) | Level 1 and Level 2 solver packages, the Level 3 transport interface, a Python/Numba cake search, and an offline viewer for the original regression. |

The [catalog](../../templates.json) lists the available entries. To copy a
scientific template, follow the [copying contract](../../PACKAGING.org).
Copy the whole mBLA package for Level 3 so its provider example can build
against the Level 2 solver. The documentation index has separate editable
presets for [a research report and slides or a Python package guide](../README.md#documentation-families).

The [wb-python package documentation](../../wb-python/README.md) covers the interface preset for future
wb-applications integrations that script Session Profile and Agent Profile
operations. The public catalog entry contains the manual and references.
Running the worked request-readback example requires a separately supplied
wb-python runtime and a configured Sandbox. The optional `wb-jupyter`
application belongs to `wb-applications`; `wb-client` is the separate
Host-side client.

## Worked use cases {#worked-use-cases}

### Inspect the embedded cylinder regression

Open the [mBLA profile explorer](../../wb-mbla/modern-lv2/tools/profile-explorer.html) in a browser.
It starts offline with the embedded cylinder regression at CASE2 (`b = 1`,
nominal transmembrane pressure (TMP) 1.15 kPa). In **Axial profile**, choose
a quantity such as wall concentration and move the **Stored axial station**
slider through the 200 stored stations. Beside it, **Outlet radial profile**
has its own **Stored radial sample** slider for the 400 samples at `z/L = 1`.
The radial readout gives wall distance `y/R`, radius from the axis `r/R`,
and `y/δ_CP` for the case's fixed concentration-polarization length.
The tables give the archive and replay values and signed differences.

Open **Site hollow-fibre reference** to investigate the saved Python
suspension at nominal TMP 5 kPa. Choose a point on the shaded cutaway to
see concentration and the selected velocity component along both directions.
This site preview uses its own fixed CP length and a separate dataset from
the JCP archive comparison and the internal coupling benchmark.

To inspect a compatible file, select **Open profile data**. The
[data guide](../../wb-mbla/modern-lv2/docs/PROFILE_EXPLORER.md) explains the
formats. The explorer displays saved results locally in the browser.

### Update a solver or exchange a quantity

Start from one traceable calculation you trust. Save the original inputs,
recorded outputs and tolerances. In [BD–MC](../../wb-bd-mc/docs/VALIDATION.md), the
portable legacy executable is checked with frozen regressions, and the modern
components have checks for selected topology, transition and clock operations.
The separate [equilibrium EOS example](../../wb-bd-mc/docs/EQUILIBRIUM_EOS.org)
exports observations from a native hard-sphere sampler. A Python process reads
them and estimates pressure. The run covers equilibrium pressure and near-contact
structure under the stated model conditions. For your own exchange, define the
physical quantity and check the supplier and receiving calculation together.

For a small comparison with retained results, open the
[BD–MC equilibrium EOS example](../../examples/bdmc-eos/README.md). It compares
soft-repulsive Brownian dynamics and hard-sphere Monte Carlo with the
Carnahan–Starling equation. The example includes the measured CSV, both
figure formats and the scope of the finite-size checks.

### Read measurements, make figures and share the comparison

Inspect a CSV header and select the numeric columns explicitly with
[Scientific CSV I/O](../../scientific-io/README.md). Choose units, transformations
and axis limits according to the measurements. Apply [plot styles](../../plot-styles/README.md)
to use consistent axes and switch between color and monochrome. Use
[wb-style](../../wb-style/README.md) when you want the separate preset derived
from Park's original plotting scripts. The [collection plotting example](../../README.md#try-plotting)
generates its own synthetic CSV and writes PNG, PDF and rendering metadata.

To communicate the result, copy the [Quarto comparison example](../examples/quarto-comparison/README.md).
Use its report to explain quantities, cases and acceptance criteria in full, and
its slides to show the selected figure and link back to the report. Replace the
synthetic data and source comparison with the inputs from your own study.

## Reproduce and adapt {#reproduce-and-adapt}

1. **Choose the scope.** Decide what you want to reproduce or extend. Read the
   selected template's requirements, test procedure and evidence limits.
2. **Make an independent copy.** Keep the source, small fixtures, tests, license,
   notices and provenance together. The [copying contract](../../PACKAGING.org)
   and [collection commands](../../README.md#copy-verify-and-adapt) describe the
   local exporter and manifest check.
3. **Run the supplied example in the stated environment.** Use a fresh output
   directory. Record the source edition, inputs, command and result so that you
   can rerun the same case later.
4. **Change one part at a time.** Keep your original comparison. For a solver,
   record changes to assumptions, boundary conditions, units, numerical methods or
   estimators. For a figure or document, keep the input data distinct from style
   and layout settings.
5. **Check what you will rely on.** Compare the adapted result with a relevant
   reference or invariant under your own conditions. For an exchanged quantity,
   check that the receiving component interprets the values and applicable range
   as intended.

## Evidence and limits {#evidence-and-limits}

- **BD–MC.** The [modernization record](../../wb-bd-mc/docs/MODERNIZATION_LEVELS.org)
  reports nine frozen current-port regressions, including driven shear cases with
  changes in network topology. The detached Level 3 check covers a native sampler
  and independent pressure analysis. Full modern active-shear correspondence,
  scientific validation of RPY and Stokesian dynamics, and broader post-processing
  validation remain open. The equilibrium MC sweep counts attempted moves.
- **mBLA.** The [modernization overview](../../wb-mbla/docs/MODERNIZATION.md) records
  eight selected cases on the 200 × 400 grid. Level 2 matched Level 1 files
  byte-for-byte and the archived mean pressures and permeate velocities at stored
  precision. A separate non-cake 5 kPa, 10 nm case checked live and pretabulated
  transport providers. The evidence applies to the specified cases, quantities
  and models. Broader operating conditions and mesh convergence need their own
  checks.
- **wb-python.** The [validation reference](../../wb-python/docs/VALIDATION.org)
  gives the tested runtime editions and the scope of the request and Profile checks.
- **Plotting, I/O and publishing.** The supplied examples use synthetic data to
  exercise plotting, numeric selection and publishing. Their evidence covers
  those operations. Original research figures and physical results have separate
  sources and acceptance criteria.

Use the [collection license scope](../../LICENSE_SCOPE.md) and the selected
template's notices when copying code, data or rendered assets.
