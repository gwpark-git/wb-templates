# Quarto comparison documentation template

[AI-generated text](#ai-assistance)

Use this example to publish a long, navigable HTML report from multiple
Markdown chapters and a short Reveal.js deck. Both show the same synthetic
diffusivity, viscosity and osmotic-pressure comparison, with a real source
diff and a linked CSV. Edit the report and slides in a copied folder.

## Render the example

Install Quarto in your chosen environment and make `quarto` available on
`PATH`. From this template folder, with an existing parent for the build:

```sh
sh render.sh /tmp/wb-comparison-build
```

Choose a fresh build path on each run. Open `site/report.html` for the
continuous document with outline and code-copy controls. Open
`site/slides.html` for the presentation. Both render from local assets.
To save a presentation PDF, open `site/slides.html?print-pdf` in Chromium,
print with background graphics, and turn off printer headers and footers.

The Markdown chapter files (prefixed with `_`) are listed in [report.qmd](report.qmd).
[slides.qmd](slides.qmd) names stable slide IDs and the shared figure.
`_quarto.yml` selects render targets and sidecar resources. The inputs
include both PDF and PNG figures and an explicit-license scope.

## Adapt the worked example

Replace the [reference model](models/reference.py) and
[candidate model](models/candidate.py) and set the concentration grid,
quantity labels and tolerance in [tools/make_example.py](tools/make_example.py).
The generator requires Python 3.11+, Matplotlib and the Workbench
`wb_plot_styles` package from the `plot-styles` template. Select the installed
package or add that template checkout to `PYTHONPATH`. The example uses
`personal` in monochrome mode, with red reserved for the revised series.
Physical labels use Mathtext, including `$\eta/\eta_s$`. From the template folder:

```sh
python3 -B tools/make_example.py /tmp/wb-comparison-data
```

The generator stages scientific-notation CSV, data table, source diff and two
figure formats beneath the new destination. Compare them with your intended
model and copy the named outputs into your working copy before rerendering.
The bundled case uses synthetic polynomials throughout.

## Edit the report or slides

Use **Edit** in either rendered view to see which files to open. In the
presentation it opens the editing instructions slide in the same deck. Edit
[slides.qmd](slides.qmd) for slide text, [report.qmd](report.qmd) for report
structure, and the named Markdown chapter for its content. Rebuild from
your edited files. In the existing User-in-Loop Review Conversation, share
the edited source or identify the source file, heading and edition.

For portable PDF links, use the optional Python 3.11+ helper with `pypdf`
available. Pass the exact local deck URL used for printing and a new output.

```sh
python3 tools/portable_pdf.py /tmp/wb-comparison-build/site/slides.html \
  /tmp/slides.browser.pdf /tmp/slides.pdf
```

Keep the resulting PDF beside the rendered report and CSV so their relative
links resolve in PDF readers that support relative file links. See
[rendered asset notices](THIRD_PARTY_NOTICES.md) for browser library terms.

## Author and project information

**Author and Project Director:** Gunwoo Park. [License](LICENSE)

Project overview: <https://github.com/gwpark-git/agentic-workbench-knowledge>  
wb-templates repository (planned): <https://github.com/gwpark-git/wb-templates>

### AI assistance

OpenAI Codex assisted with documentation drafting and editing under Gunwoo
Park's direction. Publication uses wb-plugins Publication Support and
User-in-Loop Review.
