# Working on this package

Read `README.md`, `docs/AI_UPDATE_SCOPE.md` and `docs/VALIDATION.md` before changing scientific code or claiming validation. Read `NOTICE.md` before redistributing imported material.

Before extending a component, read `docs/CODE_GUIDE.org`, `docs/MODERNIZATION_LEVELS.org` and the relevant checks in `docs/VALIDATION.md`. Org is plain text. Its `*` and `**` headings correspond to Markdown `#` and `##`. Open linked sources directly when they govern the work.

Preserve human and upstream attribution, frozen numerical fixtures and explicit model/units/geometry/estimator contracts. Consult the original source and a separate physical oracle when a result disagrees. Record the tests actually performed; passing finite cases does not establish exhaustive simulation coverage. Keep source changes separate from generated build/run output and do not publish or create a new baseline without the operator's scope.
