"""Small modernized legend helpers; source lineage/rights are in NOTICE.md."""
from __future__ import annotations

import numpy as np
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.markers import MarkerStyle
from matplotlib.colors import is_color_like


class HandlerXoffset(HandlerLine2D):
    """Retain the source's horizontal shift for grouped legend symbols."""

    def __init__(self, marker_pad=.3, numpoints=1, x_offset=0, **kwargs):
        """Store a finite horizontal displacement in Matplotlib legend coordinates."""
        super().__init__(marker_pad=marker_pad, numpoints=numpoints, **kwargs)
        if not np.isfinite(x_offset):
            raise ValueError("x_offset must be finite")
        self.x_offset = x_offset

    def get_xdata(self, legend, xdescent, ydescent, width, height, fontsize):
        """Translate line and marker coordinates by the same signed offset."""
        xdata, markers = super().get_xdata(legend, xdescent, ydescent, width, height, fontsize)
        return np.asarray(xdata) - self.x_offset, np.asarray(markers) - self.x_offset


def draw_matrix_legend(ax, row_labels, col_labels, shapes, colors, *, fontsize=12,
                       xscale=.5, yscale=1.3, row_header_text=None, border=False,
                       gaps=(.1, .25, .02), ms=12, facecolor="none", mec="k", alpha=1):
    """Draw the source's table-aligned symbol legend; return created artists.

    Input dimensions/colors/markers are checked before artists are created.
    None colors deliberately leave cells empty. No plotted data are transformed.
    """
    nrows, ncols = len(row_labels), len(col_labels)
    if nrows < 1 or ncols < 1 or len(shapes) != ncols or len(colors) != nrows or any(len(row) != ncols for row in colors):
        raise ValueError("legend needs nonempty labels, one shape/column and a matching colors matrix")
    for shape in shapes:
        MarkerStyle(shape)
    if any(color is not None and not is_color_like(color) for row in colors for color in row):
        raise ValueError("invalid matrix legend color")
    if len(gaps) != 3 or any(not np.isfinite(v) or v <= 0 for v in (*gaps, fontsize, xscale, yscale, ms)):
        raise ValueError("gaps and legend dimensions must be positive finite values")
    if not 0 <= alpha <= 1:
        raise ValueError("alpha must be in [0,1]")
    table = ax.table(cellText=[[""] * ncols for _ in range(nrows)], rowLabels=row_labels,
                     cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)
    table.scale(xscale, yscale)
    for cell in table.get_celld().values():
        cell.set_linewidth(0)
    for i in range(nrows):
        table[(i, -1)].get_text().set_ha("right")
        table[(i, -1)].set_width(gaps[0])
        for j in range(ncols):
            table[(i, j)].set_width(gaps[1] if j == 0 else gaps[2])
    # Finalize table cell geometry before positioning markers and headers.
    ax.figure.canvas.draw()
    markers, headers = [], []
    for i in range(nrows):
        for j in range(ncols):
            if colors[i][j] is None:
                continue
            cell = table[(i, j)]
            x, y = cell.get_xy()
            marker, = ax.plot(x + cell.get_width() / 2, y + cell.get_height() / 2,
                              marker=shapes[j], linestyle="none", ms=ms, mfc=colors[i][j],
                              mec=mec, mew=1.2, alpha=alpha, transform=ax.transAxes,
                              clip_on=False, zorder=3)
            markers.append(marker)
    for j, label in enumerate(col_labels):
        cell = table[(0, j)]
        x, y = cell.get_xy()
        headers.append(ax.text(x + cell.get_width() / 2, y + 1.2 * cell.get_height(), label,
                               ha="center", va="bottom", transform=ax.transAxes, zorder=4, fontsize=fontsize))
    if row_header_text:
        cell = table[(0, -1)]
        x, y = cell.get_xy()
        headers.append(ax.text(x + cell.get_width() / 2, y + 1.2 * cell.get_height(), row_header_text,
                               ha="center", va="bottom", transform=ax.transAxes, zorder=4, fontsize=fontsize))
    patch = None
    if border:
        # Measure table/text in display pixels before converting the padded frame to axes units.
        from matplotlib.patches import Rectangle
        from matplotlib.transforms import Bbox
        renderer = ax.figure.canvas.get_renderer()
        bbox = Bbox.union([table.get_window_extent(renderer)] + [text.get_window_extent(renderer) for text in headers])
        bbox = bbox.expanded(1.06, 1.18)
        bbox = Bbox.from_extents(bbox.x0 - 10, bbox.y0 - 8, bbox.x1 + 6, bbox.y1 + 8)
        bbox_ax = ax.transAxes.inverted().transform_bbox(bbox)
        patch = Rectangle((bbox_ax.x0, bbox_ax.y0), bbox_ax.width, bbox_ax.height,
                           facecolor=facecolor, edgecolor="k", lw=1, transform=ax.transAxes, zorder=2)
        ax.add_patch(patch)
    return {"table": table, "markers": markers, "headers": headers, "border": patch}
