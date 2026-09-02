from __future__ import annotations

from typing import Union
from xml.etree import ElementTree as ET

import numpy as np

from .._colormap import colormap_lookup
from .._facecolor import (
    ColorSpec,
    faces_flat,
    faces_from_projection,
    resolve_face_colors,
)
from .._projection import ISO_AZIMUTH, ISO_ELEVATION, project_surface
from .._skin import _extract_skin_py, _has_skinnable_cells

# The extra viewBox width a colorbar occupies, as a multiple of the drawing's
# width. Only the viewBox grows: the scaling factor, the stroke width and every
# mesh coordinate come from the unmodified width. Twins of the kSvgColorbar*
# constants in src/cpp/src/formats/svg.cpp.
_COLORBAR_GAP = 0.08
_COLORBAR_WIDTH = 0.06
_COLORBAR_LABELS = 0.30
_COLORBAR_SEGMENTS = 32


def _color_hex(rgb) -> str:
    """SVG's colour vocabulary. Lowercase hex, as C's "%02x" also produces."""
    r, g, b = rgb
    return "#%02x%02x%02x" % (r, g, b)


def _view_width(width: float, colorbar: bool) -> float:
    if not colorbar:
        return width
    return width * (1.0 + _COLORBAR_GAP + _COLORBAR_WIDTH + _COLORBAR_LABELS)


def _append_colorbar(svg, colors, min_x, min_y, width, height, float_fmt):
    """A vertical gradient to the right of the drawing plus min/max labels.

    ``<rect>``, not ``<path>``, because the document-level ``path {...}`` rule
    would otherwise stroke every swatch. Attribute insertion order is fixed and
    matches the C++ emission order (ElementTree writes attributes in insertion
    order). Twin of ``svg_append_colorbar``.
    """
    num = f"{{:{float_fmt}}}"
    x0 = min_x + width + width * _COLORBAR_GAP
    bar_w = width * _COLORBAR_WIDTH
    seg_h = height / float(_COLORBAR_SEGMENTS)
    for i in range(_COLORBAR_SEGMENTS):
        y0 = min_y + height * (float(i) / float(_COLORBAR_SEGMENTS))
        # y grows downward in SVG, so the FIRST segment is the top one and must
        # carry the HIGH end of the range.
        c = colormap_lookup(
            colors.table,
            float(_COLORBAR_SEGMENTS - 1 - i) / float(_COLORBAR_SEGMENTS - 1),
        )
        ET.SubElement(
            svg,
            "rect",
            x=num.format(x0),
            y=num.format(y0),
            width=num.format(bar_w),
            height=num.format(seg_h),
            fill=_color_hex(c),
        )
    label_x = x0 + bar_w + width * 0.02
    font_size = width * 0.04
    hi = ET.SubElement(
        svg,
        "text",
        x=num.format(label_x),
        y=num.format(min_y + font_size),
        **{"font-size": num.format(font_size)},
    )
    hi.text = num.format(colors.vmax)
    lo = ET.SubElement(
        svg,
        "text",
        x=num.format(label_x),
        y=num.format(min_y + height),
        **{"font-size": num.format(font_size)},
    )
    lo.text = num.format(colors.vmin)


def write(
    filename,
    mesh,
    float_fmt: str = ".3f",
    stroke_width: Union[str, None] = None,
    # Use a default image_width (not None). If set to None, images will come out at the
    # width of the mesh (which is okay). Some viewers (e.g., eog) have problems
    # displaying SVGs of width around 1 since they interpret it as the width in pixels.
    image_width: Union[int, float, None] = 100,
    # ParaView's default colors
    fill: str = "#c8c5bd",
    stroke: str = "#000080",
    # Orthographic camera for genuinely 3D input (ignored for flat meshes);
    # the default is the classic CAD isometric view.
    azimuth: float = ISO_AZIMUTH,
    elevation: float = ISO_ELEVATION,
    roll: float = 0.0,
    # Data-driven colouring. `color_by` unset keeps the flat `fill` and leaves
    # the output byte-identical to previous releases.
    color_by: Union[str, None] = None,
    component: Union[int, None] = None,
    cmap: str = "viridis",
    vmin: Union[float, None] = None,
    vmax: Union[float, None] = None,
    nan_color: str = "#808080",
    colorbar: bool = False,
):
    spec = ColorSpec(color_by, component, cmap, vmin, vmax, colorbar)

    if mesh.points.shape[1] == 3 and not np.allclose(
        mesh.points[:, 2], 0.0, rtol=0.0, atol=1.0e-14
    ):
        # Genuinely non-flat 3D mesh: render the projected skin instead.
        return _write_projected(
            filename,
            mesh,
            float_fmt,
            stroke_width,
            image_width,
            fill,
            stroke,
            azimuth,
            elevation,
            roll,
            spec,
            nan_color,
        )

    pts = mesh.points[:, :2].copy()
    min_x = np.min(pts[:, 0]) if len(pts) > 0 else 0.0
    max_x = np.max(pts[:, 0]) if len(pts) > 0 else 0.0
    min_y = np.min(pts[:, 1]) if len(pts) > 0 else 0.0
    max_y = np.max(pts[:, 1]) if len(pts) > 0 else 0.0

    pts[:, 1] = max_y + min_y - pts[:, 1]

    width = max_x - min_x
    height = max_y - min_y

    if image_width is not None and width != 0:
        scaling_factor = image_width / width
        min_x *= scaling_factor
        min_y *= scaling_factor
        width *= scaling_factor
        height *= scaling_factor
        pts *= scaling_factor

    if stroke_width is None:
        # "%g", not str(): the C++ writer formats this with snprintf("%g", ...)
        # and the two must agree byte for byte (tests/python/test_svg.py pins it).
        stroke_width = "%g" % (width / 100)

    # On the flat path the drawn mesh IS the source mesh, so a face's cell
    # index needs no provenance indirection.
    colors = resolve_face_colors(spec, mesh, mesh, faces_flat(mesh))
    show_bar = colors.active and colorbar

    fmt = " ".join(4 * [f"{{:{float_fmt}}}"])
    svg = ET.Element(
        "svg",
        xmlns="http://www.w3.org/2000/svg",
        version="1.1",
        viewBox=fmt.format(min_x, min_y, _view_width(width, show_bar), height),
    )

    style = ET.SubElement(svg, "style")
    opts = [
        f"fill: {fill}",
        f"stroke: {stroke}",
        f"stroke-width: {stroke_width}",
        "stroke-linejoin:bevel",
    ]
    # Use path, not polygon, because svgo converts polygons to paths and doesn't convert
    # the style alongside. No problem if it's paths all the way.
    style.text = "path {" + "; ".join(opts) + "}"

    face_index = 0
    for cell_block in mesh.cells:
        if cell_block.type not in ["line", "triangle", "quad"]:
            continue

        if cell_block.type == "line":
            fmt = (
                f"M {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
            )
        elif cell_block.type == "triangle":
            fmt = (
                f"M {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
                + "Z"
            )
        elif cell_block.type == "quad":
            fmt = (
                f"M {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
                + "Z"
            )
        for cell in cell_block.data:
            f = face_index
            face_index += 1
            # A per-path fill attribute alone does NOT override the
            # document-level rule: per the CSS cascade, a presentation
            # attribute has zero specificity and loses to ANY stylesheet
            # rule, even the plain "path {fill: ...}" type selector above -
            # every renderer that honours the cascade (browsers, cairosvg)
            # paints every face in the flat fallback colour regardless of
            # this attribute. The redundant inline `style` wins the cascade
            # and actually colours the face; `fill` stays for inspection
            # and because tests key off it. Without colouring, neither
            # attribute is emitted and output is unchanged.
            if colors.active and cell_block.type != "line":
                c = colors.color(f)
                hex_color = _color_hex(c) if c is not None else nan_color
                ET.SubElement(
                    svg,
                    "path",
                    d=fmt.format(*pts[cell].flatten()),
                    fill=hex_color,
                    style=f"fill:{hex_color}",
                )
            else:
                ET.SubElement(
                    svg,
                    "path",
                    d=fmt.format(*pts[cell].flatten()),
                )

    if show_bar:
        _append_colorbar(svg, colors, min_x, min_y, width, height, float_fmt)

    tree = ET.ElementTree(svg)
    tree.write(filename)


def _write_projected(
    filename,
    mesh,
    float_fmt,
    stroke_width,
    image_width,
    fill,
    stroke,
    azimuth,
    elevation,
    roll,
    spec,
    nan_color,
):
    # 3D rendering path: extract the skin of any volume cells, project with
    # the orthographic camera, and paint the faces back-to-front. Mirrors the
    # C++ core's svg_proj_write (bbox -> y-flip -> scaling -> emission).
    if _has_skinnable_cells(mesh):
        # Colouring by cell data needs to know which input cell each skin facet
        # came from; only ask for it when it will actually be read.
        need_parents = spec.active and spec.color_by in mesh.cell_data
        draw_mesh = _extract_skin_py(
            mesh, linearize=True, record_parent_ids=need_parents
        )
    else:
        draw_mesh = mesh
    x, y, faces = project_surface(draw_mesh, azimuth, elevation, roll)

    colors = resolve_face_colors(spec, mesh, draw_mesh, faces_from_projection(faces))
    show_bar = colors.active and spec.colorbar

    min_x = np.min(x) if len(x) > 0 else 0.0
    max_x = np.max(x) if len(x) > 0 else 0.0
    min_y = np.min(y) if len(y) > 0 else 0.0
    max_y = np.max(y) if len(y) > 0 else 0.0

    # Flip y (projected math convention y-up -> SVG screen convention y-down).
    y = max_y + min_y - y

    width = max_x - min_x
    height = max_y - min_y

    if image_width is not None and width != 0:
        scaling_factor = image_width / width
        min_x *= scaling_factor
        min_y *= scaling_factor
        width *= scaling_factor
        height *= scaling_factor
        x = x * scaling_factor
        y = y * scaling_factor

    if stroke_width is None:
        # "%g", not str(): the C++ writer formats this with snprintf("%g", ...)
        # and the two must agree byte for byte (tests/python/test_svg.py pins it).
        stroke_width = "%g" % (width / 100)

    fmt = " ".join(4 * [f"{{:{float_fmt}}}"])
    svg = ET.Element(
        "svg",
        xmlns="http://www.w3.org/2000/svg",
        version="1.1",
        viewBox=fmt.format(min_x, min_y, _view_width(width, show_bar), height),
    )

    style = ET.SubElement(svg, "style")
    opts = [
        f"fill: {fill}",
        f"stroke: {stroke}",
        f"stroke-width: {stroke_width}",
        "stroke-linejoin:bevel",
    ]
    style.text = "path {" + "; ".join(opts) + "}"

    point_fmt = f"{{:{float_fmt}}} {{:{float_fmt}}}"
    for f, (nodes, is_line, _) in enumerate(faces):
        parts = []
        for k, p in enumerate(nodes):
            parts.append(("M " if k == 0 else "L ") + point_fmt.format(x[p], y[p]))
        d = "".join(parts)
        if not is_line:
            d += "Z"
        # See the twin comment in write(): the inline `style` is what
        # actually wins the CSS cascade against the document-level rule;
        # `fill` stays for inspection/test compatibility.
        if colors.active and not is_line:
            c = colors.color(f)
            hex_color = _color_hex(c) if c is not None else nan_color
            ET.SubElement(svg, "path", d=d, fill=hex_color, style=f"fill:{hex_color}")
        else:
            ET.SubElement(svg, "path", d=d)

    if show_bar:
        _append_colorbar(svg, colors, min_x, min_y, width, height, float_fmt)

    tree = ET.ElementTree(svg)
    tree.write(filename)
