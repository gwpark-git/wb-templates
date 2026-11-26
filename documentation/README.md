# wb-templates documentation

[AI-generated text](collection-guide/report.html#ai-assistance)

[mBLA validation and site profiles](../wb-mbla/modern-lv2/tools/profile-explorer.html) /
[BD–MC pressure comparison](../examples/bdmc-eos/README.md)

## Read on GitHub or in a browser

On GitHub, use this Markdown index and the Markdown and Org source pages
linked from the [collection README](../README.md). For the rendered version,
open the distribution's [HTML home page](../index.html) through a static web
server. Document links in the HTML reading view open the corresponding HTML
page. **Source** opens the maintained Markdown or Org file. **Edit** in
the formatted report, Python guide and slides opens the editing instructions
within that document.

Start with the [collection guide](collection-guide/report.html) for the
project's purpose, practical workflows and the current template catalog.
Its [functional reference](collection-guide/report.html#functional-reference)
links the code, configuration, data contracts and checks for each area.

## Documentation families

| Family | Use it for | Preview | Editable preset |
| --- | --- | --- | --- |
| Basic documentation | Project reports, research workflows and presentations | [HTML report](examples/quarto-comparison/report.html), [Reveal.js slides](examples/quarto-comparison/slides.html), [PDF slides](examples/quarto-comparison/slides.pdf) | [Quarto comparison](examples/quarto-comparison/README.md) |
| Python package guide | Installation, first successful run, API and troubleshooting | [Read a retained Workbench request with wb-python](examples/python-package/preview.html) | [wb-python worked example](examples/python-package/README.md) |

The public previews open without a Workbench runtime. To execute the Python
request-readback example, use a separately supplied runtime and existing
request receipt. The [wb-python manual](../wb-python/README.md) explains setup.

- [mBLA explorer guide](../wb-mbla/modern-lv2/docs/PROFILE_EXPLORER.md) for the saved cases, coordinate definitions and file formats
- [BD–MC equilibrium EOS example](../examples/bdmc-eos/README.md) for retained Carnahan–Starling comparison data
- [Collection catalog and entry points](../README.md#templates)

## Open a local copy in a browser

Keep the full distribution tree together. It contains the PDF, CSV and source
files linked from the previews. If a local file preview restricts navigation,
serve the copied distribution from the machine running the browser. From its root:

```sh
python3 -m http.server --bind 127.0.0.1 8765
```

Then open <http://127.0.0.1:8765/>. The home page, previews and related assets
also work under a site prefix.

[Author and project information](collection-guide/report.html#author-and-project-information)
