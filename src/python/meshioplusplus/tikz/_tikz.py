from __future__ import annotations

from typing import Union

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
from .._files import open_file
from .._projection import ISO_AZIMUTH, ISO_ELEVATION, project_surface
from .._skin import _extract_skin_py, _has_skinnable_cells

# Number of gradient swatches the colorbar is built from (twin of
# kTikzColorbarSegments).
_COLORBAR_SEGMENTS = 32


def _marker_lines(projected, azimuth, elevation, roll):
    """The comment that opens the file: the provenance tag, the ``Drawing: ...``
    line, then whatever else the provenance block carries, each as a ``%`` line
    (twin of ``tikz_marker_lines``)."""
    return [
        "% " + line
        for line in _drawing.comment_lines(
            _drawing.marker_line(projected, azimuth, elevation, roll)
        )
    ]


def _color_rgb(rgb) -> str:
    """TikZ's inline colour vocabulary.

    The braces are required: without them the commas inside would split the
    surrounding key=value option list.
    """
    r, g, b = rgb
    return f"{{rgb,255:red,{r};green,{g};blue,{b}}}"


def _fill_style(fill_color: str, draw: str, line_width) -> str:
    """The per-path option list for a filled face.

    Factored out so the constant-fill and per-face-colour paths build it
    identically; with colouring off this reproduces the previous string
    exactly (twin of ``tikz_fill_style``).
    """
    opts = [f"fill={fill_color}", f"draw={draw}"]
    if line_width is not None:
        opts.append(f"line width={line_width}")
    return ", ".join(opts)


def _append_colorbar(lines, colors, min_x, max_x, min_y, max_y, float_fmt):
    """A vertical gradient bar to the right of the drawing, plus min/max labels.

    Twin of ``tikz_append_colorbar``; emitted only when asked for, so the
    default output is untouched.
    """
    num = f"{{:{float_fmt}}}"
    w = max_x - min_x
    h = max_y - min_y
    x0 = max_x + w * 0.08
    x1 = x0 + w * 0.06
    for i in range(_COLORBAR_SEGMENTS):
        y0 = min_y + h * (float(i) / float(_COLORBAR_SEGMENTS))
        y1 = min_y + h * (float(i + 1) / float(_COLORBAR_SEGMENTS))
        c = colormap_lookup(colors.table, float(i) / float(_COLORBAR_SEGMENTS - 1))
        lines.append(
            f"  \\fill[{_color_rgb(c)}] ({num.format(x0)},{num.format(y0)}) "
            f"rectangle ({num.format(x1)},{num.format(y1)});"
        )
    lines.append(
        f"  \\node[anchor=west, font=\\tiny] at ({num.format(x1)},"
        f"{num.format(min_y)}) {{{num.format(colors.vmin)}}};"
    )
    lines.append(
        f"  \\node[anchor=west, font=\\tiny] at ({num.format(x1)},"
        f"{num.format(max_y)}) {{{num.format(colors.vmax)}}};"
    )


def write(
    filename,
    mesh,
    float_fmt: str = ".6f",
    # Emit a full, directly-compilable LaTeX document by default. Set to False to
    # emit only the bare `tikzpicture` environment for \input into a larger document.
    standalone: bool = True,
    # TikZ line width for the edges, e.g. "0.4pt". None -> use TikZ's default.
    line_width: Union[str, None] = None,
    # xcolor spec for the filled faces (triangles/quads) and the edge stroke.
    fill: str = "gray!30",
    draw: str = "black",
    # Optional \begin{tikzpicture}[scale=...]. None -> omit (coordinates verbatim).
    scale: Union[int, float, None] = None,
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
    nan_color: str = "gray",
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
            standalone,
            line_width,
            fill,
            draw,
            scale,
            azimuth,
            elevation,
            roll,
            spec,
            nan_color,
        )

    # TikZ/PGF uses the math convention (y grows upward), so unlike the SVG
    # writer there is no y-flip: the first two columns map straight to (x, y).
    pts = mesh.points[:, :2]

    coord = f"({{:{float_fmt}}},{{:{float_fmt}}})"

    # Per-path style options.
    fill_style = _fill_style(fill, draw, line_width)
    line_opts = [f"draw={draw}"]
    if line_width is not None:
        line_opts.append(f"line width={line_width}")
    line_style = ", ".join(line_opts)

    # On the flat path the drawn mesh IS the source mesh, so a face's cell
    # index needs no provenance indirection.
    colors = resolve_face_colors(spec, mesh, mesh, faces_flat(mesh))

    lines: list[str] = []
    face_index = 0
    skipped_types = []
    for cell_block in mesh.cells:
        if cell_block.type not in ["line", "triangle", "quad"]:
            if cell_block.type not in skipped_types:
                skipped_types.append(cell_block.type)
            continue

        for cell in cell_block.data:
            path = " -- ".join(coord.format(x, y) for x, y in pts[cell])
            f = face_index
            face_index += 1
            if cell_block.type == "line":
                lines.append(f"  \\draw[{line_style}] {path};")
            elif colors.active:
                c = colors.color(f)
                fill_color = _color_rgb(c) if c is not None else nan_color
                style = _fill_style(fill_color, draw, line_width)
                lines.append(f"  \\draw[{style}] {path} -- cycle;")
            else:
                # triangle / quad: closed, filled face
                lines.append(f"  \\draw[{fill_style}] {path} -- cycle;")

    if skipped_types:
        warn(
            "TikZ: cell type(s) "
            + ", ".join(skipped_types)
            + " are not representable (only line/triangle/quad); skipping."
        )

    if colors.active and colorbar and len(mesh.points) > 0:
        _append_colorbar(
            lines,
            colors,
            float(np.min(pts[:, 0])),
            float(np.max(pts[:, 0])),
            float(np.min(pts[:, 1])),
            float(np.max(pts[:, 1])),
            float_fmt,
        )

    pic_opts = []
    if scale is not None:
        pic_opts.append(f"scale={scale}")
    if line_width is not None:
        pic_opts.append(f"line width={line_width}")
    pic_opt_str = f"[{', '.join(pic_opts)}]" if pic_opts else ""

    body = [f"\\begin{{tikzpicture}}{pic_opt_str}", *lines, "\\end{tikzpicture}"]

    if standalone:
        out = [
            "\\documentclass{standalone}",
            "\\usepackage{tikz}",
            "\\begin{document}",
            *body,
            "\\end{document}",
        ]
    else:
        out = body
    out = [*_marker_lines(False, 0.0, 0.0, 0.0), *out]

    with open_file(filename, "w", newline="\n") as f:
        f.write("\n".join(out) + "\n")


def _write_projected(
    filename,
    mesh,
    float_fmt,
    standalone,
    line_width,
    fill,
    draw,
    scale,
    azimuth,
    elevation,
    roll,
    spec,
    nan_color,
):
    # 3D rendering path: extract the skin of any volume cells, project with
    # the orthographic camera, and draw the faces back-to-front. Mirrors the
    # C++ core's tikz_proj_write byte for byte.
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

    coord = f"({{:{float_fmt}}},{{:{float_fmt}}})"

    fill_style = _fill_style(fill, draw, line_width)
    line_opts = [f"draw={draw}"]
    if line_width is not None:
        line_opts.append(f"line width={line_width}")
    line_style = ", ".join(line_opts)

    lines: list[str] = []
    for f, (nodes, is_line, _) in enumerate(faces):
        path = " -- ".join(coord.format(x[p], y[p]) for p in nodes)
        if is_line:
            lines.append(f"  \\draw[{line_style}] {path};")
        elif colors.active:
            c = colors.color(f)
            fill_color = _color_rgb(c) if c is not None else nan_color
            lines.append(
                f"  \\draw[{_fill_style(fill_color, draw, line_width)}] {path} -- cycle;"
            )
        else:
            lines.append(f"  \\draw[{fill_style}] {path} -- cycle;")

    if colors.active and spec.colorbar and len(x) > 0:
        _append_colorbar(
            lines,
            colors,
            float(np.min(x)),
            float(np.max(x)),
            float(np.min(y)),
            float(np.max(y)),
            float_fmt,
        )

    pic_opts = []
    if scale is not None:
        pic_opts.append(f"scale={scale}")
    if line_width is not None:
        pic_opts.append(f"line width={line_width}")
    pic_opt_str = f"[{', '.join(pic_opts)}]" if pic_opts else ""

    body = [f"\\begin{{tikzpicture}}{pic_opt_str}", *lines, "\\end{tikzpicture}"]

    if standalone:
        out = [
            "\\documentclass{standalone}",
            "\\usepackage{tikz}",
            "\\begin{document}",
            *body,
            "\\end{document}",
        ]
    else:
        out = body
    out = [*_marker_lines(True, azimuth, elevation, roll), *out]

    with open_file(filename, "w", newline="\n") as f:
        f.write("\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# Reader. Inverts exactly what the writer above emits, and nothing more; the
# twin of ``read_tikz`` (src/cpp/src/formats/tikz.cpp) -- see ``_drawing.py``.
# ---------------------------------------------------------------------------

_LETTERS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _construct_at(text, pos):
    """A name for the construct at the cursor, for a refusal: the next run of
    non-blank characters, cut short."""
    end = pos
    while end < len(text) and text[end] not in _drawing._BLANKS and end - pos < 16:
        end += 1
    return _drawing.safe(text[pos:end])


def _split_options(text):
    """Splits an option list on the commas outside braces."""
    items = []
    depth = 0
    start = 0
    for i, c in enumerate(text):
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
        elif c == "," and depth == 0:
            items.append(text[start:i])
            start = i + 1
    items.append(text[start:])
    return items


def _strip(text):
    return text.strip(_drawing._BLANKS)


def _rgb(value):
    """``{rgb,255:red,R;green,G;blue,B}`` -> ``(R, G, B)``; anything else
    (``gray!30``, ...) -> ``None``."""
    value = _strip(value)
    if len(value) < 2 or value[0] != "{" or value[-1] != "}":
        return None
    value = value[1:-1]
    head = "rgb,255:red,"
    if not value.startswith(head):
        return None
    value = value[len(head) :]
    rgb = []
    for sep in (";green,", ";blue,", ""):
        end = 0
        while end < len(value) and "0" <= value[end] <= "9":
            end += 1
        if end == 0 or end > 3:
            return None
        v = int(value[:end])
        if v > 255:
            return None
        rgb.append(v)
        value = value[end:]
        if not value.startswith(sep):
            return None
        value = value[len(sep) :]
    if value:
        return None
    return tuple(rgb)


def _picture_options(options):
    """``scale=`` and ``line width=`` only: anything else changes the drawing."""
    for item in _split_options(options):
        item = _strip(item)
        if not item:
            continue
        key = _strip(item.partition("=")[0])
        if key not in ("scale", "line width"):
            raise ReadError(
                f"TikZ: the tikzpicture option '{_drawing.safe(key)}' is not supported"
            )


def _draw_options(options):
    """``fill``, ``draw`` and ``line width``; the fill is returned when it is an
    RGB one."""
    fill = None
    for item in _split_options(options):
        item = _strip(item)
        if not item:
            continue
        key, eq, value = item.partition("=")
        key = _strip(key)
        if key not in ("fill", "draw", "line width"):
            raise ReadError(
                f"TikZ: the \\draw option '{_drawing.safe(key)}' is not supported"
            )
        if key == "fill" and eq:
            fill = _rgb(value)
    return fill


class _Reader:
    """One reading of a figure: the text with its comments cut, and a cursor."""

    def __init__(self, text, builder):
        self.text = text
        self.pos = 0
        self.builder = builder

    def run(self):
        text = self.text
        in_document = False
        in_picture = False
        seen_picture = False
        while True:
            self._skip_blanks()
            if self.pos >= len(text):
                break
            if text[self.pos] != "\\":
                raise ReadError(
                    f"TikZ: unexpected text '{_construct_at(text, self.pos)}'"
                )
            self.pos += 1
            end = self.pos
            while end < len(text) and text[end] in _LETTERS:
                end += 1
            name = text[self.pos : end]
            if not name:
                raise ReadError(
                    "TikZ: unexpected control sequence "
                    f"'\\{_construct_at(text, self.pos)}'"
                )
            self.pos = end
            if name in ("documentclass", "usepackage"):
                self._skip_blanks()
                if self.pos < len(text) and text[self.pos] == "[":
                    self._balanced("[", "]")
                self._skip_blanks()
                self._balanced("{", "}")
            elif name in ("begin", "end"):
                self._skip_blanks()
                env = self._balanced("{", "}")
                if env == "document":
                    if in_document if name == "begin" else not in_document:
                        raise ReadError("TikZ: unbalanced document environment")
                    in_document = name == "begin"
                elif env == "tikzpicture":
                    if name == "begin":
                        if in_picture or seen_picture:
                            raise ReadError("TikZ: only one tikzpicture is read")
                        self._skip_blanks()
                        if self.pos < len(text) and text[self.pos] == "[":
                            _picture_options(self._balanced("[", "]"))
                        in_picture = seen_picture = True
                    else:
                        if not in_picture:
                            raise ReadError(
                                "TikZ: \\end{tikzpicture} without its begin"
                            )
                        in_picture = False
                else:
                    raise ReadError(
                        f"TikZ: the environment '{_drawing.safe(env)}' is not supported"
                    )
            elif name == "draw":
                self._require_picture(in_picture, name)
                self._draw()
            elif name == "fill":
                self._require_picture(in_picture, name)
                self._colorbar_swatch()
            elif name == "node":
                self._require_picture(in_picture, name)
                self._skip_statement()
            else:
                raise ReadError(f"TikZ: the command \\{name} is not supported")
        if in_picture:
            raise ReadError("TikZ: the tikzpicture is not closed")
        if not seen_picture:
            raise ReadError("TikZ: no tikzpicture")

    def _skip_blanks(self):
        text = self.text
        while self.pos < len(text) and text[self.pos] in _drawing._BLANKS:
            self.pos += 1

    @staticmethod
    def _require_picture(in_picture, name):
        if not in_picture:
            raise ReadError(f"TikZ: \\{name} outside a tikzpicture")

    def _balanced(self, open_, close):
        """The text between the bracket at the cursor and its match (braces nest
        and hide brackets); the cursor ends just past the closing one."""
        text = self.text
        if self.pos >= len(text) or text[self.pos] != open_:
            raise ReadError(f"TikZ: expected '{open_}'")
        start = self.pos + 1
        braces = 0
        for i in range(start, len(text)):
            c = text[i]
            if c == "{":
                braces += 1
            elif c == "}" and braces > 0:
                braces -= 1
            elif c == close and braces == 0:
                self.pos = i + 1
                return text[start:i]
        raise ReadError(f"TikZ: '{open_}' is never closed")

    def _coordinate(self, keep=True):
        """``(x,y)`` -> the welded vertex; x and y are plain numbers, nothing
        else. ``keep`` false (a colour-bar corner): checked, not welded."""
        text = self.text
        self._skip_blanks()
        if self.pos >= len(text):
            raise ReadError("TikZ: a path ends where a coordinate is expected")
        if text[self.pos] == "+":
            raise ReadError("TikZ: relative coordinates (+, ++) are not supported")
        if text[self.pos] != "(":
            raise ReadError(
                f"TikZ: unsupported path construct '{_construct_at(text, self.pos)}'"
            )
        close = text.find(")", self.pos)
        if close < 0:
            raise ReadError("TikZ: a coordinate is never closed")
        inside = text[self.pos + 1 : close]
        self.pos = close + 1
        first, comma, second = inside.partition(",")
        x = _strip(first)
        y = _strip(second) if comma else ""
        if (
            not comma
            or not x
            or not y
            or _drawing.NUMBER.fullmatch(x) is None
            or _drawing.NUMBER.fullmatch(y) is None
        ):
            raise ReadError(
                f"TikZ: the coordinate ({_drawing.safe(inside)}) is not a plain (x,y) pair of "
                "numbers (named, polar, calc, 3-D and unit-carrying coordinates "
                "are not supported)"
            )
        if not keep:
            _drawing.number_value(x, "TikZ")
            _drawing.number_value(y, "TikZ")
            return -1
        return self.builder.weld(
            x,
            y,
            lambda: (
                _drawing.number_value(x, "TikZ"),
                _drawing.number_value(y, "TikZ"),
            ),
        )

    def _draw(self):
        """``\\draw[opts] (x,y) -- (x,y) ... [-- cycle];``"""
        text = self.text
        self._skip_blanks()
        fill = None
        if self.pos < len(text) and text[self.pos] == "[":
            fill = _draw_options(self._balanced("[", "]"))
        ids = []
        closed = False
        while True:
            ids.append(self._coordinate())
            self._skip_blanks()
            if self.pos < len(text) and text[self.pos] == ";":
                self.pos += 1
                break
            if text[self.pos : self.pos + 2] != "--":
                raise ReadError(
                    "TikZ: unsupported path construct "
                    f"'{_construct_at(text, self.pos)}'"
                )
            self.pos += 2
            self._skip_blanks()
            if text[self.pos : self.pos + 5] == "cycle" and (
                self.pos + 5 >= len(text) or text[self.pos + 5] not in _LETTERS
            ):
                self.pos += 5
                closed = True
                self._skip_blanks()
                if self.pos >= len(text) or text[self.pos] != ";":
                    raise ReadError("TikZ: expected ';' after cycle")
                self.pos += 1
                break
        self.builder.add_cell(ids, closed, fill)

    def _colorbar_swatch(self):
        """``\\fill[colour] (x,y) rectangle (x,y);`` -- the colour bar's swatches."""
        text = self.text
        self._skip_blanks()
        if self.pos < len(text) and text[self.pos] == "[":
            self._balanced("[", "]")
        self._coordinate(False)
        self._skip_blanks()
        if text[self.pos : self.pos + 9] != "rectangle":
            raise ReadError(
                "TikZ: \\fill is only read as a colour-bar rectangle, found "
                f"'{_construct_at(text, self.pos)}'"
            )
        self.pos += 9
        self._coordinate(False)
        self._skip_blanks()
        if self.pos >= len(text) or text[self.pos] != ";":
            raise ReadError("TikZ: expected ';' after the rectangle")
        self.pos += 1

    def _skip_statement(self):
        """``\\node ... ;`` -- the colour bar's labels: up to the next ';'
        outside braces."""
        text = self.text
        braces = 0
        while self.pos < len(text):
            c = text[self.pos]
            if c == "{":
                braces += 1
            elif c == "}" and braces > 0:
                braces -= 1
            elif c == ";" and braces == 0:
                self.pos += 1
                return
            self.pos += 1
        raise ReadError("TikZ: a \\node is never closed")


def read(filename, points_only=False, arrays=None, strict=True):
    """Read a TikZ figure written by :func:`write` back as a planar mesh.

    The twin of ``read_tikz``; the format page has the details.
    """
    if hasattr(filename, "read"):
        data = filename.read()
    else:
        with open(filename, "rb") as f:
            data = f.read()
    # One character per byte, as the C++ reader sees the file.
    text = data.decode("latin-1") if isinstance(data, (bytes, bytearray)) else data

    # Comments are cut at the first unescaped '%'. The marker is the comment
    # lines that lead the file, before any code.
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()  # no empty last line after a final newline
    body = []
    head = []
    leading = True
    for line in lines:
        if line.endswith("\r"):
            line = line[:-1]
        cut = -1
        i = 0
        while i < len(line):
            if line[i] == "\\":
                i += 1
            elif line[i] == "%":
                cut = i
                break
            i += 1
        code = line if cut < 0 else line[:cut]
        if not _strip(code):
            if cut >= 0 and leading:
                head.append(_strip(line[cut:].lstrip("%")))
        else:
            leading = False
        body.append(code + "\n")

    marker = _drawing.parse_marker(head, "TikZ", False)
    if not marker.found and strict:
        raise _drawing.unmarked_error("TikZ")

    builder = _drawing.Builder("TikZ", "tikz:fill")
    builder.want_fill = not points_only and (arrays is None or "tikz:fill" in arrays)
    _Reader("".join(body), builder).run()

    # Coordinates are verbatim: no flip, no scale.
    mesh = builder.finish(False, 1.0, 0.0, 0.0)
    _drawing.read_warning("TikZ", marker)
    return mesh
