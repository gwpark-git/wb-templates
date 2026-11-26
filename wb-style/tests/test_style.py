"""Behavioral checks with neutral synthetic fixtures only.

"""
import importlib.util
from pathlib import Path
from shutil import which
import tempfile
import unittest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from wb_style import configuration, dash_pattern, font_selection, style_context

path = Path(__file__).resolve().parents[1] / "examples/plot_example.py"
spec = importlib.util.spec_from_file_location("neutral_style_example", path)
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


class StyleBehavior(unittest.TestCase):
    """Check style restoration, actual artists, CSV data integrity and portable rendering."""
    def test_context_restores_normal_and_error_exit(self):
        """Compare complete caller rcParams after normal and exceptional exits."""
        before = matplotlib.rcParams.copy()
        with style_context():
            self.assertEqual(matplotlib.rcParams["xtick.direction"], "in")
            self.assertEqual(matplotlib.rcParams["axes.linewidth"], 1.5)
        self.assertEqual(dict(matplotlib.rcParams), dict(before))
        with self.assertRaisesRegex(RuntimeError, "chosen failure"):
            with style_context():
                raise RuntimeError("chosen failure")
        self.assertEqual(dict(matplotlib.rcParams), dict(before))

    def test_actual_artists_and_unchanged_raw_values(self):
        """Read figure size, strokes and markers from the rendered Matplotlib artists."""
        data = example.synthetic_series()
        d = configuration()["plot_defaults"]
        with style_context():
            fig, ax = plt.subplots(figsize=d["figsize_inches"])
            try:
                line, = ax.plot(data[:, 0], data[:, 1], color=d["line_color"], lw=d["line_width_pt"])
                marker, = ax.plot(data[::8, 0], data[::8, 1], ls="none", marker="o", mfc=d["marker_face_color"], ms=d["marker_size_pt"], mew=d["marker_edge_width_pt"])
                np.testing.assert_array_equal(fig.get_size_inches(), [8, 6])
                np.testing.assert_array_equal(line.get_ydata(), data[:, 1])
                self.assertEqual(line.get_linewidth(), 3)
                self.assertEqual(marker.get_markerfacecolor(), "white")
                self.assertEqual(marker.get_markersize(), 8)
                self.assertEqual(ax.spines["left"].get_linewidth(), 1.5)
            finally:
                plt.close(fig)

    def test_configuration_is_fresh_and_dash_patterns_are_explicit(self):
        """Mutate a configuration copy and compare named dashes and font outcomes."""
        d = configuration()
        d["rc_params"]["axes.linewidth"] = 99
        self.assertEqual(configuration()["rc_params"]["axes.linewidth"], 1.5)
        self.assertEqual(dash_pattern("densely dashed"), (0, (5, 1)))
        with self.assertRaises(KeyError):
            dash_pattern("unknown")
        self.assertIn(font_selection()["selected"], ["Arial", "DejaVu Sans"])

    def test_neutral_csv_and_actual_nonuniform_coordinates(self):
        """Require exact CSV roundtrip and marker selection on original nonuniform x values."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "series.csv"
            x = np.linspace(0, 1, 31) ** 2
            data = np.column_stack((x, 2 * x + 1, 3 - x))
            np.savetxt(path, data, delimiter=",", header="x,series_a,series_b", comments="", fmt="%.17g")
            actual = example.load_series(path)
            np.testing.assert_array_equal(actual, data)
            indices = example.marker_indices(len(data))
            np.testing.assert_array_equal(indices, [0, 8, 16, 24, 30])
            np.testing.assert_array_equal(actual[indices, 0], x[[0, 8, 16, 24, 30]])

    def test_invalid_csv_refuses_before_output(self):
        """Exercise malformed and non-finite input and assert no output directory exists."""
        for text in ["x,other\n0,1\n1,2\n", "x,series_a,series_b\n0,nan,1\n1,2,3\n", "x,series_a,series_b\n0,1\n1,2,3\n", "x,series_a,series_b\n0,1,2\n"]:
            with self.subTest(text=text), tempfile.TemporaryDirectory() as directory:
                source = Path(directory) / "input.csv"
                source.write_text(text)
                output = Path(directory) / "output"
                with self.assertRaises(ValueError):
                    example.render(output, source)
                self.assertFalse(output.exists())

    def test_real_render_roundtrip_and_no_overwrite(self):
        """Check PNG/PDF artifacts, exact numeric readback and existing-file preservation."""
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "demo"
            summary = example.render(output)
            self.assertEqual(summary["data_origin"], "synthetic")
            self.assertEqual(summary["rows"], 81)
            self.assertGreater((output / "series.png").stat().st_size, 1000)
            self.assertTrue((output / "series.pdf").read_bytes().startswith(b"%PDF"))
            np.testing.assert_array_equal(example.load_series(output / "series.csv"), example.synthetic_series())
            before = (output / "series.csv").read_bytes()
            with self.assertRaises(FileExistsError):
                example.render(output)
            self.assertEqual((output / "series.csv").read_bytes(), before)

    def test_default_needs_no_tex(self):
        """Verify Mathtext works and a missing TeX executable leads to explicit refusal."""
        with style_context() as rendering:
            self.assertEqual(rendering["text_engine"], "Mathtext")
        if which("latex") is None:
            with self.assertRaisesRegex(RuntimeError, "requires latex"):
                with style_context(tex=True):
                    pass


if __name__ == "__main__":
    unittest.main()
