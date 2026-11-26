# wb-templates documentation

[AI-generated text](../README.md#ai-assistance)

## Collection guide

- [Purpose, setup and worked scientific examples](collection-guide/_guide.md)
- [Functional reference for code, configuration, data and checks](collection-guide/_reference.md)
- [Template catalog](../README.md#templates) and [copying contract](../PACKAGING.org)
- [Validation and reproduction notes](../provenance/VALIDATION.org)

## Scientific templates

- [mBLA modernization](../wb-mbla/README.md), [inlet velocity validation](../wb-mbla/docs/INLET_VELOCITY.org),
  [profile guide](../wb-mbla/modern-lv2/docs/PROFILE_EXPLORER.md),
  and the [interactive regression and site examples](../wb-mbla/modern-lv2/tools/profile-explorer.html)
- [BD–MC modernization](../wb-bd-mc/docs/MODERNIZATION_LEVELS.org), [build and use](../wb-bd-mc/docs/BUILD_AND_USAGE.md),
  and the [Carnahan–Starling pressure comparison](../examples/bdmc-eos/README.md)
- [Plot styles](../plot-styles/README.md), [wb-style](../wb-style/README.md),
  and [scientific CSV I/O](../scientific-io/README.md)
- [wb-python reader guide and API reference](../wb-python/README.md)

## Documentation families

| Family | Read or copy | Bundled browser example |
| --- | --- | --- |
| Basic documentation | [Quarto report and slide project](examples/quarto-comparison/README.md), [guide chapter](examples/quarto-comparison/_guide.md), [methods](examples/quarto-comparison/_methods.md), [results](examples/quarto-comparison/_results.md), [source comparison](examples/quarto-comparison/_code.md) | [HTML report](examples/quarto-comparison/report.html), [Reveal.js slides](examples/quarto-comparison/slides.html), [PDF slides](examples/quarto-comparison/slides.pdf) |
| Python package guide | [Worked request-readback guide](examples/python-package/guide.md), [render instructions](examples/python-package/README.md) | [Formatted Python guide](examples/python-package/preview.html) |

The public previews open without a Workbench runtime. To execute the Python
request-readback example, use a separately supplied runtime and an existing
request receipt. The [wb-python manual](../wb-python/README.md) explains setup.

## Build and read the HTML edition

Follow the [README build command](../README.md#build-the-html-site). The
[formatted collection handbook](collection-guide/report.html) is also bundled.
Open the generated site from its output directory so relative links to PDF, CSV and
source files resolve. **Source** opens the maintained Markdown or Org file.
**Edit** in a formatted report, Python guide or slide deck opens the
editing instructions within that document.

[Author and project information](../README.md#author-and-project-information)
