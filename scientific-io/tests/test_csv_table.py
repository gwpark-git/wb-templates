"""MIT-licensed generic CSV behavior checks; no scientific source fixtures."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

from wb_scientific_io import FormatError
from wb_scientific_io.csv_table import inspect_csv, read_csv


class CSVBehavior(unittest.TestCase):
    """Exercise the CSV reader on synthetic text, malformed records and numeric boundary cases."""
    def setUp(self):
        """Allocate an isolated temporary path for each fixture."""
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "data with spaces.csv"

    def write(self, text):
        """Write the selected CSV string without adding implicit rows."""
        self.path.write_text(text, encoding="utf-8")
        return self.path

    def test_actual_header_shape_type_summary_and_selective_numeric_read(self):
        """Compare actual headers, type counts and finite arrays in requested order."""
        path = self.write('time,pressure,label\n0,10,"condition, A"\n1,12,condition B\n')
        before = path.read_bytes()
        info = inspect_csv(path)
        self.assertEqual(info["header"], ["time", "pressure", "label"])
        self.assertEqual(info["shape"], [2, 3])
        self.assertTrue(info["columns"]["time"]["numeric"])
        self.assertEqual(info["columns"]["label"]["invalid_numeric_count"], 2)
        result = read_csv(path, columns=["pressure", "time"])
        self.assertEqual(list(result["columns"]), ["pressure", "time"])
        np.testing.assert_array_equal(result["columns"]["pressure"], [10, 12])
        self.assertEqual(result["columns"]["time"].dtype, np.float64)
        self.assertFalse(result["metadata"]["unit_inference"])
        self.assertEqual(path.read_bytes(), before)

    def test_missing_unrelated_text_is_reported_but_selected_data_survive(self):
        """Allow unrelated text or missing cells and verify the selected result."""
        path = self.write("x,y,label\n0,1,\n1,2,ok\n")
        result = read_csv(path, columns=["x", "y"])
        self.assertEqual(result["metadata"]["columns"]["label"]["missing_count"], 1)
        self.assertEqual(result["metadata"]["selected_shape"], [2, 2])

    def test_selected_missing_invalid_and_nonfinite_fail_with_column(self):
        """Require column-specific refusal for missing, text and non-finite selected values."""
        for value in ("", "text", "NaN", "Inf", "-Inf"):
            with self.subTest(value=value), self.assertRaisesRegex(FormatError, "selected|missing"):
                read_csv(self.write(f"x,y\n0,{value}\n"), columns=["x", "y"])

    def test_duplicate_empty_headers_ragged_and_malformed_fail(self):
        """Reject invalid names, widths and CSV quoting before numeric selection."""
        for text in ("x,x\n0,1", " x ,x\n0,1", "x,\n0,1", "x,y\n0", 'x,y\n0,"unfinished'):
            with self.subTest(text=text), self.assertRaises(FormatError):
                inspect_csv(self.write(text))

    def test_explicit_delimiter_encoding_blank_records_and_no_data(self):
        """Check explicit parsing choices and reject a numeric selection without rows."""
        path = self.write("x;y\n\n0;1\n")
        self.assertEqual(inspect_csv(path, delimiter=";")["blank_record_count"], 1)
        self.assertEqual(read_csv(path, columns=["x"], delimiter=";")["columns"]["x"].shape, (1,))
        self.write("x,y\n")
        self.assertEqual(inspect_csv(path)["row_count"], 0)
        with self.assertRaises(FormatError):
            read_csv(path, columns=["x"])

    def test_missing_file_column_selection_and_unsupported_policy(self):
        """Report missing input, unknown columns and unsupported missing policies."""
        with self.assertRaises(FileNotFoundError):
            inspect_csv(self.path)
        path = self.write("x,y\n0,1\n")
        for names in ([], "x", ["x", "x"], ["absent"]):
            with self.subTest(names=names), self.assertRaises(FormatError):
                read_csv(path, columns=names)
        with self.assertRaises(FormatError):
            read_csv(path, columns=["x"], missing="drop")
        with self.assertRaises(FormatError):
            inspect_csv(path, delimiter="::")


if __name__ == "__main__":
    unittest.main()
