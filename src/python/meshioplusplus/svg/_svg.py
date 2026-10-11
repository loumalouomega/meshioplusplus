from __future__ import annotations

from typing import Union
from xml.etree import ElementTree as ET

import numpy as np

from .. import _drawing
from .._colormap import colormap_lookup
from .._common import warn
from .._exceptions import ReadError
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


def _write_tree(tree, filename):
    """Writes the tree with LF line ends on every platform. ElementTree opens a
    filename in text mode, which turns the newlines of the marker comment into
    CRLF on Windows and breaks byte identity with the C++ writer, which writes
    binary. A binary file object keeps the newlines as they are."""
    with open(filename, "wb") as f:
        tree.write(f)


def _append_marker(svg, projected, azimuth, elevation, roll, scale, offset_x, offset_y):
    """The marker comment right after the opening tag: the provenance tag, the
    ``Drawing: ...`` line, then whatever else the provenance block carries.
    Always the multi-line form, with "--" broken apart as the XML comment rules
    need (twin of ``svg_marker_comment``)."""
    line = _drawing.marker_line(
        projected, azimuth, elevation, roll, (scale, offset_x, offset_y)
    )
    svg.append(ET.Comment(_drawing.xml_comment_text(_drawing.comment_lines(line))))


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

    y_flip = max_y + min_y
    pts[:, 1] = max_y + min_y - pts[:, 1]

    width = max_x - min_x
    height = max_y - min_y

    applied_scale = 1.0
    if image_width is not None and width != 0:
        scaling_factor = image_width / width
        applied_scale = scaling_factor
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

    _append_marker(svg, False, 0.0, 0.0, 0.0, applied_scale, 0.0, y_flip)

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
    skipped_types = []
    for cell_block in mesh.cells:
        if cell_block.type not in ["line", "triangle", "quad"]:
            if cell_block.type not in skipped_types:
                skipped_types.append(cell_block.type)
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

    if skipped_types:
        warn(
            "SVG: cell type(s) "
            + ", ".join(skipped_types)
            + " are not representable (only line/triangle/quad); skipping."
        )

    tree = ET.ElementTree(svg)
    _write_tree(tree, filename)


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
    y_flip = max_y + min_y
    y = max_y + min_y - y

    width = max_x - min_x
    height = max_y - min_y

    applied_scale = 1.0
    if image_width is not None and width != 0:
        scaling_factor = image_width / width
        applied_scale = scaling_factor
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

    _append_marker(svg, True, azimuth, elevation, roll, applied_scale, 0.0, y_flip)

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
    _write_tree(tree, filename)


# ---------------------------------------------------------------------------
# Reader. Inverts exactly what the writer above emits, and nothing more; the
# twin of ``read_svg`` (src/cpp/src/formats/svg.cpp) -- see ``_drawing.py``.
# ---------------------------------------------------------------------------

_SVG_NS = "http://www.w3.org/2000/svg"


class _Node:
    """One element of the parsed document: its tag, attributes and children
    (elements, and comments as ``(text,)`` tuples)."""

    __slots__ = ("tag", "attrib", "kids")

    def __init__(self, tag, attrib):
        self.tag = tag
        self.attrib = attrib
        self.kids = []


class _Target:
    """A parser target that keeps the root's comments and element children, and
    notes whether the root declares a namespace prefix."""

    def __init__(self):
        self.stack = []
        self.root = None
        self.prefixed = False

    def start(self, tag, attrib):
        node = _Node(tag, attrib)
        if self.stack:
            self.stack[-1].kids.append(node)
        else:
            self.root = node
        self.stack.append(node)
        return node

    def end(self, tag):
        self.stack.pop()

    def data(self, data):
        pass

    def comment(self, text):
        if self.stack:
            self.stack[-1].kids.append((text,))

    def start_ns(self, prefix, uri):
        if prefix and not self.stack:
            self.prefixed = True

    def end_ns(self, prefix):
        pass

    def close(self):
        return self.root


def _local(node):
    """``(namespace, local name)`` of an element tag."""
    tag = node.tag
    if tag.startswith("{"):
        uri, _, name = tag[1:].partition("}")
        return uri, name
    return "", tag


def _check_element(node):
    uri, name = _local(node)
    if "transform" in node.attrib:
        raise ReadError(
            f"SVG: a transform on <{_drawing.safe(name)}> is not supported "
            "(meshio++ draws absolute coordinates)"
        )
    if uri not in ("", _SVG_NS):
        raise ReadError(
            f"SVG: <{_drawing.safe(name)}> is in a namespace other than SVG's"
        )


def _view_box(text):
    """The four numbers of a viewBox."""
    box = []
    pos = 0
    for _ in range(4):
        while pos < len(text) and (text[pos] in _drawing._BLANKS or text[pos] == ","):
            pos += 1
        m = _drawing.NUMBER.match(text, pos)
        if m is None:
            raise ReadError("SVG: malformed viewBox")
        box.append(_drawing.number_value(m.group(), "SVG"))
        pos = m.end()
    while pos < len(text) and (text[pos] in _drawing._BLANKS or text[pos] == ","):
        pos += 1
    if pos != len(text):
        raise ReadError("SVG: malformed viewBox")
    return box


def _read_path(d, builder, fill):
    """One path's ``d``: ``M x y (L x y)* Z?``, absolute, one subpath."""
    pos = 0
    n = len(d)

    def skip():
        nonlocal pos
        while pos < n and (d[pos] in _drawing._BLANKS or d[pos] == ","):
            pos += 1

    def pair():
        nonlocal pos
        skip()
        mx = _drawing.NUMBER.match(d, pos)
        if mx is None:
            raise ReadError("SVG: a path command is missing its coordinates")
        pos = mx.end()
        skip()
        my = _drawing.NUMBER.match(d, pos)
        if my is None:
            raise ReadError("SVG: a path command is missing its coordinates")
        pos = my.end()
        x, y = mx.group(), my.group()
        return builder.weld(
            x,
            y,
            lambda: (_drawing.number_value(x, "SVG"), _drawing.number_value(y, "SVG")),
        )

    ids = []
    started = False
    closed = False
    while True:
        skip()
        if pos >= n:
            break
        c = d[pos]
        if closed:
            raise ReadError(
                "SVG: a path continues after Z (several subpaths are not supported)"
            )
        if c == "M":
            if started:
                raise ReadError("SVG: a path with several subpaths is not supported")
            pos += 1
            ids.append(pair())
            started = True
        elif c == "L":
            if not started:
                raise ReadError("SVG: a path must start with M")
            pos += 1
            ids.append(pair())
        elif c == "Z":
            if not started:
                raise ReadError("SVG: a path must start with M")
            pos += 1
            closed = True
        elif c in "mlhvcsqtaz":
            raise ReadError(
                f"SVG: the relative path command '{c}' is not supported "
                "(meshio++ writes absolute M, L and Z)"
            )
        elif c in "CSQTA":
            raise ReadError(f"SVG: the curve command '{c}' is not supported")
        elif c in "HV":
            raise ReadError(f"SVG: the path command '{c}' is not supported")
        elif _drawing.NUMBER.match(d, pos) is not None:
            raise ReadError(
                "SVG: repeated coordinates after one path command are not supported"
            )
        else:
            raise ReadError(f"SVG: unexpected character '{_drawing.safe(c)}' in a path")
    if not started:
        raise ReadError("SVG: an empty path")
    builder.add_cell(ids, closed, fill)


def read(filename, points_only=False, arrays=None, strict=True):
    """Read an SVG drawing written by :func:`write` back as a planar mesh.

    The twin of ``read_svg``; the format page has the details.
    """
    if hasattr(filename, "read"):
        data = filename.read()
    else:
        with open(filename, "rb") as f:
            data = f.read()
    if isinstance(data, str):
        data = data.encode("utf-8")
    target = _Target()
    parser = ET.XMLParser(target=target)
    try:
        parser.feed(data)
        root = parser.close()
    except ET.ParseError as exc:
        raise ReadError(f"SVG: could not parse {filename} ({exc})") from exc
    if root is None or _local(root)[1] != "svg":
        raise ReadError("SVG: the root element is not <svg>")
    _check_element(root)
    if target.prefixed:
        raise ReadError("SVG: namespace prefixes are not supported")

    # The marker: the comments that lead the root's children.
    head = []
    for kid in root.kids:
        if isinstance(kid, _Node):
            break
        head.extend(kid[0].split("\n"))
    marker = _drawing.parse_marker(head, "SVG", True)
    if not marker.found and strict:
        raise _drawing.unmarked_error("SVG")

    view_box_attr = root.attrib.get("viewBox")
    if view_box_attr is None:
        raise ReadError("SVG: the <svg> has no viewBox")
    view_box = _view_box(view_box_attr)

    builder = _drawing.Builder("SVG", "svg:fill")
    builder.want_fill = not points_only and (arrays is None or "svg:fill" in arrays)
    for node in root.kids:
        if not isinstance(node, _Node):
            continue
        name = _local(node)[1]
        _check_element(node)
        if name == "path":
            d = node.attrib.get("d")
            if d is None:
                raise ReadError("SVG: a <path> has no d attribute")
            _read_path(d, builder, _drawing.hex_color(node.attrib.get("fill")))
        elif name in ("style", "rect", "text"):
            # The stylesheet is the writer's own; <rect> and <text> are the colour bar.
            pass
        elif name == "polygon":
            raise ReadError(
                "SVG: <polygon> is not supported (meshio++ draws faces as <path>; "
                "tools such as svgo rewrite polygons to paths)"
            )
        else:
            raise ReadError(
                f"SVG: the element <{_drawing.safe(name)}> is not supported"
            )

    # model = (x / S + OX, OY - y / S). Without a marker S = 1 and the y-flip is
    # undone from the viewBox: it was about the middle of the drawing's y extent.
    scale, offset_x, offset_y = 1.0, 0.0, 2.0 * view_box[1] + view_box[3]
    if marker.found:
        scale, offset_x, offset_y = marker.scale, marker.offset_x, marker.offset_y
    mesh = builder.finish(True, scale, offset_x, offset_y)
    _drawing.read_warning("SVG", marker)
    return mesh
