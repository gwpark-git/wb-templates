import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb_plot_styles import (COLOR_MODES, PRESETS, batlow_colors, dash_pattern, plot_defaults, style_context)
from wb_plot_styles.helpers import HandlerXoffset, draw_matrix_legend


class PresetBehavior(unittest.TestCase):
    """Check actual rc state, artists, color resources and the finite appearance API."""
    def test_one_common_style_and_removed_year_selectors_fail_before_mutation(self):
        """Require one public preset and preserve rcParams for obsolete selectors."""
        self.assertEqual(PRESETS, ("personal",))
        for name in ("park-2020", "park-2026", "2020", "2026", "jcp2020", "jcp2026"):
            with self.subTest(name=name):
                before = dict(matplotlib.rcParams)
                with self.assertRaises(ValueError):
                    with style_context(name):
                        pass
                self.assertEqual(dict(matplotlib.rcParams), before)

    def test_source_personal_values_and_cycle_without_global_registration(self):
        """Compare fixed personal settings and batlow samples without registry mutation."""
        before = set(matplotlib.colormaps)
        with style_context() as metadata:
            self.assertEqual(matplotlib.rcParams["xtick.major.width"], 1.2)
            self.assertEqual(matplotlib.rcParams["lines.linewidth"], 2.5)
            self.assertEqual(matplotlib.rcParams["axes.labelsize"], 15)
            np.testing.assert_allclose(matplotlib.rcParams["axes.prop_cycle"].by_key()["color"], batlow_colors([.2, .7]))
            self.assertEqual(metadata["defaults"]["dpi"], 300)
        self.assertEqual(set(matplotlib.colormaps), before)

    def test_color_mode_changes_color_without_changing_geometry_or_data(self):
        """Compare real Line2D, error-bar and scatter outputs for the same input geometry."""
        observed = {}
        for mode in COLOR_MODES:
            with style_context(color_mode=mode) as metadata:
                fig, ax = plt.subplots()
                lines = [ax.plot([0, 1], [i, i + 1])[0] for i in range(4)]
                error = ax.errorbar([0, 1], [0, 1], yerr=[.1, .2], fmt="o-")
                scatter = ax.scatter([0, 1], [1, 2])
                observed[mode] = (metadata["defaults"], [line.get_color() for line in lines])
                self.assertEqual(metadata["color_mode"], mode)
                np.testing.assert_array_equal(lines[0].get_ydata(), [0, 1])
                np.testing.assert_array_equal(error.lines[0].get_ydata(), [0, 1])
                if mode == "monochrome":
                    self.assertTrue(all(line.get_color() == "black" for line in lines))
                    self.assertEqual([line.get_linestyle() for line in lines], ["-", "--", "-.", ":"])
                    np.testing.assert_allclose(scatter.get_facecolors()[0], [0, 0, 0, 1])
                else:
                    np.testing.assert_allclose(lines[0].get_color(), batlow_colors([.2])[0])
                plt.close(fig)
        self.assertEqual(observed["color"][0], observed["monochrome"][0])

    def test_invalid_modes_fail_without_mutation_and_custom_dashes_remain(self):
        """Check mode rejection before mutation and preserve a specific custom dash."""
        for mode in ("2020", "bw", "unknown", None, []):
            before = dict(matplotlib.rcParams)
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                with style_context(color_mode=mode):
                    pass
            self.assertEqual(dict(matplotlib.rcParams), before)
        self.assertEqual(dash_pattern("-.."), (0, (3, 1, 1, 1, 1, 1)))

    def test_modes_restore_nested_state_and_allow_explicit_color_override(self):
        """Keep the outer style across nested use and allow explicit artist styling."""
        before = dict(matplotlib.rcParams)
        with style_context(color_mode="monochrome"):
            outer = dict(matplotlib.rcParams)
            with style_context(color_mode="color"):
                pass
            self.assertEqual(dict(matplotlib.rcParams), outer)
            fig, ax = plt.subplots()
            line, = ax.plot([0, 1], [1, 2], color="red", ls=":")
            self.assertEqual(line.get_color(), "red")
            self.assertEqual(line.get_linestyle(), ":")
            plt.close(fig)
        self.assertEqual(dict(matplotlib.rcParams), before)

    def test_friendly_and_rc_overrides_restore_after_error(self):
        """Compare friendly defaults and rc values and require restoration on exception."""
        before = dict(matplotlib.rcParams)
        with self.assertRaisesRegex(RuntimeError, "example"):
            with style_context(overrides={"markersize": 12, "linewidth": 4, "figsize": [5, 3], "axes.linewidth": 2}) as metadata:
                self.assertEqual(matplotlib.rcParams["lines.markersize"], 12)
                self.assertEqual(matplotlib.rcParams["figure.figsize"], [5, 3])
                self.assertEqual(metadata["defaults"]["markersize"], 12)
                raise RuntimeError("example")
        self.assertEqual(dict(matplotlib.rcParams), before)

    def test_every_preset_restores_arbitrary_caller_style_and_is_nested(self):
        """Test restoration with nondefault caller font and axes settings."""
        for name in PRESETS:
            with self.subTest(name=name), matplotlib.rc_context({"axes.linewidth": 7, "font.size": 31}):
                before = dict(matplotlib.rcParams)
                with style_context(name):
                    outer = dict(matplotlib.rcParams)
                    with style_context("personal", {"lines.markersize": 17}) as metadata:
                        self.assertEqual(metadata["defaults"]["markersize"], 17)
                    self.assertEqual(dict(matplotlib.rcParams), outer)
                self.assertEqual(dict(matplotlib.rcParams), before)

    def test_rejects_nonappearance_unknown_and_invalid_overrides_without_mutation(self):
        """Require rejection before mutation for invalid, unsafe or nonappearance keys."""
        invalid = [{"xlim": [0, 1]}, {"backend": "Agg"}, {"text.usetex": True},
                   {"markersize": -1}, {"dpi": float("nan")}, {"figsize": [2]},
                   {"lines.color": "not a color"}, {"bbox_inches": "whatever"}]
        for override in invalid:
            with self.subTest(override=override):
                before = dict(matplotlib.rcParams)
                with self.assertRaises(ValueError):
                    with style_context(overrides=override):
                        pass
                self.assertEqual(dict(matplotlib.rcParams), before)
        with self.assertRaises(ValueError):
            plot_defaults("unknown")

    def test_defaults_are_fresh_and_have_no_scientific_axis_limits(self):
        """Mutate one returned mapping and check the complete keys in a new result."""
        defaults = plot_defaults()
        defaults["figsize"][0] = 999
        self.assertEqual(plot_defaults()["figsize"], [6, 6])
        self.assertEqual(set(defaults), {"figsize", "dpi", "linewidth", "markersize", "labelsize", "tick_labelsize", "bbox_inches"})

    def test_batlow_rejects_nonfinite_or_outside_positions(self):
        """Reject NaN and coordinates outside the normalized color interval."""
        for positions in ([float("nan")], [-.1], [1.1]):
            with self.assertRaises(ValueError):
                batlow_colors(positions)


class HelperBehavior(unittest.TestCase):
    """Check marker and label geometry, and refusal without partial artists."""
    def test_matrix_artist_count_geometry_headers_and_optional_empty_cell(self):
        """Compare artists with actual cell centers and check empty-cell handling."""
        fig, ax = plt.subplots()
        artists = draw_matrix_legend(ax, ["one", "two"], ["A", "B"], ["o", "s"],
                                      [["red", "blue"], [None, "black"]], row_header_text="row", border=True)
        self.assertEqual(len(artists["markers"]), 3)
        self.assertEqual(len(artists["headers"]), 3)
        cell = artists["table"][(0, 0)]
        first = artists["markers"][0]
        np.testing.assert_allclose([first.get_xdata()[0], first.get_ydata()[0]],
                                    [cell.get_x() + cell.get_width() / 2, cell.get_y() + cell.get_height() / 2])
        self.assertIsNotNone(artists["border"])
        plt.close(fig)

    def test_matrix_invalid_shape_fails_before_creating_artists(self):
        """Require an empty axes table list after a dimension error."""
        fig, ax = plt.subplots()
        with self.assertRaises(ValueError):
            draw_matrix_legend(ax, ["one"], ["A", "B"], ["o"], [["red"]])
        self.assertEqual(len(ax.tables), 0)
        plt.close(fig)

    def test_offset_handler_retains_opposite_marker_positions(self):
        """Check the 12-unit difference between equal and opposite offsets."""
        fig, ax = plt.subplots()
        line, = ax.plot([0, 1], [0, 1])
        legend = ax.legend([line], ["curve"])
        left = HandlerXoffset(x_offset=6).get_xdata(legend, 0, 0, 10, 10, 15)
        right = HandlerXoffset(x_offset=-6).get_xdata(legend, 0, 0, 10, 10, 15)
        np.testing.assert_allclose(right[0] - left[0], 12)
        np.testing.assert_allclose(right[1] - left[1], 12)
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
