"""What the ``svg`` and ``tikz`` drawings share: the marker comment their writers
emit and their readers require, and the builder that turns a drawing's paths
back into a planar mesh.

The Python twin of ``src/cpp/src/detail/drawing.hpp``; its doc comment has the
design. In short: both drawings carry a comment right after the provenance tag
(``Written by meshio++ v...``) holding one line, ``Drawing: planar`` or
``Drawing: projected azimuth=A elevation=E roll=R``, and for ``svg`` also
``scale=S offset=OX,OY`` (numbers are ``%.17g``, so the reader restores them
exactly). A reader refuses a file without it unless told ``strict=False``,
because a drawing is not a mesh: only the structure meshio++'s own writers emit
is read back, and the marker is how a reader knows it is looking at one.

The model, as the ``svg`` reader restores it::

    model_x = svg_x / S + OX        model_y = OY - svg_y / S

with ``S`` the scale ``image_width`` applied (1 when none was) and ``OY`` the
``max_y + min_y`` the writer flipped y about. Without a marker ``S = 1``,
``OX = 0`` and ``OY = 2 * viewBox_y + viewBox_height``.
"""

from __future__ import annotations

import math
import re

import numpy as np

from . import _provenance
from ._common import warn
from ._exceptions import ReadError
from ._mesh import CellBlock, Mesh

#: The prefix of the marker line, shared by both drawings.
MARKER_PREFIX = "Drawing: "

#: An SVG/TikZ number: ``[+-]? (digits [. digits?] | . digits) ([eE] [+-]? digits)?``
#: -- the twin of ``drawing_scan_number``. ASCII digits only (``\d`` would take others).
NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")

_BLANKS = " \t\n\r\v\f"
_BLANK_RUN = re.compile(r"[ \t\n\r\v\f]+")


def split_blanks(text: str) -> list:
    """The tokens of ``text`` split on the C++ reader's blanks (ASCII only)."""
    return [t for t in _BLANK_RUN.split(text) if t]


def safe(text: str) -> str:
    """File text for an error message: anything but printable ASCII becomes ``?``
    (the twin of ``drawing_safe``)."""
    return "".join(c if " " <= c <= "~" else "?" for c in text)


def num17(value) -> str:
    """``%.17g`` -- what the marker prints numbers with, so the reader restores them."""
    return "%.17g" % value


def numg(value) -> str:
    """``%g`` -- what the read warning prints numbers with."""
    return "%g" % value


def marker_line(projected, azimuth, elevation, roll, scale_offset=None) -> str:
    """The marker's ``Drawing: ...`` line. ``scale_offset`` is ``(S, OX, OY)``
    or ``None`` (tikz)."""
    line = MARKER_PREFIX
    if projected:
        line += (
            f"projected azimuth={num17(azimuth)} elevation={num17(elevation)}"
            f" roll={num17(roll)}"
        )
    else:
        line += "planar"
    if scale_offset is not None:
        scale, offset_x, offset_y = scale_offset
        line += f" scale={num17(scale)} offset={num17(offset_x)},{num17(offset_y)}"
    return line


def comment_lines(line: str) -> list:
    """The provenance lines of a drawing's comment, the marker line right after
    the tag (so it is there under every provenance mode)."""
    lines = list(_provenance.lines(_provenance.SlotTier.BLOCK))
    lines.insert(1, line)
    return lines


def xml_comment_text(lines) -> str:
    """The text of an XML comment holding ``lines``, in the multi-line form with
    every ``--`` broken apart (the comment rules forbid it)."""
    out = "\n"
    for line in lines:
        while "--" in line:
            line = line.replace("--", "- -")
        out += line + "\n"
    return out


class Marker:
    """What a drawing's marker said (the twin of ``DrawingMarker``)."""

    def __init__(self):
        self.found = False
        self.projected = False
        self.azimuth = 0.0
        self.elevation = 0.0
        self.roll = 0.0
        self.scale = 1.0
        self.offset_x = 0.0
        self.offset_y = 0.0


def _marker_number(token: str, fmt: str, key: str) -> float:
    m = NUMBER.fullmatch(token)
    value = float(token) if m else None
    if value is None or not math.isfinite(value):
        raise ReadError(f"{fmt}: malformed drawing marker, {key} is not a number")
    return value


def parse_marker(lines, fmt: str, need_scale: bool) -> Marker:
    """The marker among ``lines`` (the comment text of a drawing's head).

    Found when both the provenance tag and a ``Drawing: `` line are present. A
    ``Drawing: `` line that is malformed is a ``ReadError`` (the file claims to
    be ours and is not); keys this version does not know are skipped.
    """
    tag = False
    drawing = None
    for raw in lines:
        line = raw.strip(_BLANKS)
        if line.startswith("Written by meshio++ v"):
            tag = True
        elif drawing is None and line.startswith(MARKER_PREFIX):
            drawing = line
    marker = Marker()
    if not tag or drawing is None:
        return marker
    tokens = split_blanks(drawing[len(MARKER_PREFIX) :])
    if not tokens:
        raise ReadError(f"{fmt}: malformed drawing marker, empty")
    if tokens[0] == "planar":
        marker.projected = False
    elif tokens[0] == "projected":
        marker.projected = True
    else:
        raise ReadError(
            f"{fmt}: malformed drawing marker, unknown kind '{safe(tokens[0])}'"
        )
    seen = set()
    for token in tokens[1:]:
        key, eq, value = token.partition("=")
        if not eq:
            continue
        if key == "azimuth":
            marker.azimuth = _marker_number(value, fmt, "azimuth")
        elif key == "elevation":
            marker.elevation = _marker_number(value, fmt, "elevation")
        elif key == "roll":
            marker.roll = _marker_number(value, fmt, "roll")
        elif key == "scale":
            marker.scale = _marker_number(value, fmt, "scale")
        elif key == "offset":
            first, comma, second = value.partition(",")
            if not comma:
                raise ReadError(
                    f"{fmt}: malformed drawing marker, offset needs two numbers"
                )
            marker.offset_x = _marker_number(first, fmt, "offset")
            marker.offset_y = _marker_number(second, fmt, "offset")
        else:
            continue
        seen.add(key)
    if marker.projected and not {"azimuth", "elevation", "roll"} <= seen:
        raise ReadError(
            f"{fmt}: malformed drawing marker, a projection needs azimuth, "
            "elevation, roll"
        )
    if need_scale and not {"scale", "offset"} <= seen:
        raise ReadError(
            f"{fmt}: malformed drawing marker, scale and offset are missing"
        )
    if "scale" in seen and not marker.scale > 0.0:
        raise ReadError(f"{fmt}: malformed drawing marker, scale must be > 0")
    marker.found = True
    return marker


def read_warning(fmt: str, marker: Marker) -> None:
    """The one warning a drawing read raises: what the file was, and what is lost."""
    if not marker.found:
        warn(
            f"{fmt}: no meshio++ drawing marker; read as a plain drawing in its own "
            "units (depth, scale and offset are unknown), z = 0 assumed"
        )
    elif marker.projected:
        warn(
            f"{fmt}: projected from 3-D (azimuth {numg(marker.azimuth)}, "
            f"elevation {numg(marker.elevation)}, roll {numg(marker.roll)}), "
            "depth lost; z = 0 assumed"
        )
    else:
        warn(f"{fmt}: planar drawing, z = 0 assumed")


def unmarked_error(fmt: str) -> ReadError:
    """The refusal for a drawing without a marker."""
    return ReadError(
        f"{fmt}: no meshio++ drawing marker. Only drawings written by meshio++ are "
        "read back (a drawing is not a mesh); pass strict=False (ReadOptions::"
        "mLenient) to read this file as a plain drawing anyway"
    )


def hex_color(text):
    """An ``#rgb`` or ``#rrggbb`` colour as ``(r, g, b)``, or ``None``."""
    text = text.strip(_BLANKS) if text else ""
    if len(text) not in (4, 7) or text[0] != "#":
        return None
    digits = "0123456789abcdefABCDEF"
    rgb = []
    for k in range(3):
        if len(text) == 4:
            pair = text[1 + k] * 2
        else:
            pair = text[1 + 2 * k : 3 + 2 * k]
        if pair[0] not in digits or pair[1] not in digits:
            return None
        rgb.append(int(pair, 16))
    return tuple(rgb)


def number_value(token: str, fmt: str) -> float:
    """A scanned number's value, or ``ReadError``."""
    value = float(token)
    if not math.isfinite(value):
        raise ReadError(f"{fmt}: bad number '{safe(token)}'")
    return value


class Builder:
    """Collects a drawing's paths and builds the planar mesh (the twin of
    ``DrawingBuilder``).

    Points are welded on the exact text of their coordinates, in order of first
    appearance; consecutive cells of one type and size form one block. A fill is
    kept as an ``int32`` ``(n, 3)`` cell array named by the caller, ``-1 -1 -1``
    for a cell without an RGB fill, and only when some cell has one.
    """

    def __init__(self, fmt: str, fill_name: str):
        self._fmt = fmt
        self._fill_name = fill_name
        self._index = {}
        self._x = []
        self._y = []
        self._blocks = []  # [type, width, conn, fill]
        self._any_fill = False
        self.want_fill = True

    def weld(self, x_text: str, y_text: str, parse) -> int:
        """The welded index of the vertex whose coordinate texts are given.
        ``parse()`` returns ``(x, y)`` and is called for a new vertex only."""
        key = x_text + " " + y_text
        found = self._index.get(key)
        if found is not None:
            return found
        x, y = parse()
        idx = len(self._x)
        self._x.append(x)
        self._y.append(y)
        self._index[key] = idx
        return idx

    def add_cell(self, ids, closed: bool, fill) -> None:
        """One cell. ``closed`` is a closed path (a face), otherwise a line."""
        n = len(ids)
        if not closed:
            if n != 2:
                raise ReadError(
                    f"{self._fmt}: an open path of {n} vertices is not a line "
                    "(only 2-vertex open paths are read)"
                )
            ctype = "line"
        elif n == 3:
            ctype = "triangle"
        elif n == 4:
            ctype = "quad"
        elif n > 4:
            ctype = "polygon"
        else:
            raise ReadError(f"{self._fmt}: a closed path of {n} vertices is not a face")
        if not self._blocks or self._blocks[-1][0] != ctype or self._blocks[-1][1] != n:
            self._blocks.append([ctype, n, [], []])
        block = self._blocks[-1]
        block[2].extend(ids)
        if fill is not None and self.want_fill:
            self._any_fill = True
            block[3].extend(fill)
        else:
            block[3].extend((-1, -1, -1))

    def finish(self, flip: bool, scale: float, offset_x: float, offset_y: float):
        """The mesh, points mapped back to model units.

        With ``flip``: ``x / scale + offset_x``, ``offset_y - y / scale``.
        Without it and with the identity scale/offset the coordinates are left
        bit-exact.
        """
        x = np.array(self._x, dtype=np.float64)
        y = np.array(self._y, dtype=np.float64)
        if not (not flip and scale == 1.0 and offset_x == 0.0 and offset_y == 0.0):
            x = x / scale + offset_x
            y = offset_y - y / scale if flip else y / scale + offset_y
        if not (np.isfinite(x).all() and np.isfinite(y).all()):
            raise ReadError(f"{self._fmt}: a coordinate is not finite")
        points = np.zeros((len(x), 3), dtype=np.float64)
        points[:, 0] = x
        points[:, 1] = y
        cells = []
        fills = []
        for ctype, width, conn, fill in self._blocks:
            n = len(conn) // width
            cells.append(
                CellBlock(ctype, np.array(conn, dtype=np.int64).reshape(n, width))
            )
            if self._any_fill:
                fills.append(np.array(fill, dtype=np.int32).reshape(n, 3))
        cell_data = {self._fill_name: fills} if self._any_fill else {}
        return Mesh(points, cells, cell_data=cell_data)
