# Python package guide example

[AI-generated text](guide.md#ai-assistance)

This example uses the reader layout of the wb-python manual. It shows how a
research script selects an existing Workbench conversation and uses its saved
receipt to check a previous request and read a completed result. The HTML has
a linked outline, copy controls and an Edit route to the Markdown source.
In the public collection, open the preview from the [documentation index](../../README.md).
The preview and local render command need no Workbench runtime or receipt.
Executing the readback example needs an existing request receipt, a configured
Sandbox and a separately supplied wb-python runtime. Use the browser's Print
action to save a PDF.

Edit [guide.md](guide.md), or start a new guide with the
[blank starter](../../templates/python-package-guide/starter.md). From this
example directory, with Quarto on `PATH`, render to a new filename beside the source.

```sh
sh ../../templates/python-package-guide/render.sh guide.md reader.html
```

Open the new local `reader.html`. To rebuild into the same name, move or remove that generated
HTML first. The template uses Quarto's Pandoc without running the documented
Python or shell blocks. Keep PDFs and other resources at the relative locations
used by the Markdown. On an HTTP site in a subdirectory, a relative link such as
`figures/comparison.pdf` keeps the site prefix.

The four front matter fields `title`, `author`, `version` and `date`
provide the title, author metadata, package version and documentation edition
date. `lang` sets the document language. The
[reader layout](../../templates/python-package-guide/reader.html) owns the
responsive outline, print styling and copy controls. Its layout and copy script
are adapted from the MIT-licensed wb-python Reader Guide. See the
[documentation template license](../../../LICENSE).
