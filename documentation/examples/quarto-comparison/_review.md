## Verification and provenance {#verification-and-provenance}

The example data are stored in scientific notation and paired by volume
fraction. The source listing compares the two named input modules. The
build uses Quarto with execution disabled. A separate run of the sample
data tool generates the comparison figures and values.

### Edit the source files {#edit-source}

| To change | Edit |
|---|---|
| A slide | `slides.qmd`. Find the slide heading. Keep its existing section ID when editing a slide that has been cited. |
| The report structure or a section | `report.qmd` lists the chapters. Edit `_guide.md`, `_methods.md`, `_results.md`, `_code.md` or `_review.md` for that section. |
| The comparison values or figure | Edit `models/reference.py`, `models/candidate.py` or the generator `tools/make_example.py`, then regenerate the CSV, figure and generated chapters. |

From the template folder, run `sh render.sh /tmp/wb-comparison-updated` with
a new destination. Share the edited source or identify the file, section and
edition in the existing User-in-Loop Review Conversation.

**Publication method** wb-plugins Publication Support and User-in-Loop Review.
