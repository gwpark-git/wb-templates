# wb scientific I/O — CSV

Inspect a CSV header before selecting columns. The reader returns finite
numeric values in the order you specify, with descriptive metadata.

[CSV contract](docs/CSV.org) / [License](LICENSE)

With Python 3.10+ and NumPy 1.24+, install from this folder:

```sh
python3 -m pip install .
```

Use your own file and column names in this example:

```python
from wb_scientific_io.csv_table import inspect_csv, read_csv
info = inspect_csv("measurement.csv")
table = read_csv("measurement.csv", columns=["time", "pressure"])
```

Values remain unchanged, column order is explicit and units are not inferred.
Malformed/ragged headers or selected missing/nonfinite/text values fail with
`FormatError`. Unrelated text columns can remain in the same file. This is an
I/O contract. Use the bundled synthetic cases to check it with
`python3 -B -m unittest discover -s tests -v`.
License text and attribution are retained; no original data/library is bundled.
