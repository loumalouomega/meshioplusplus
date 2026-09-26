"""I/O for OpenRadioss / Radioss starter decks (``*_0000.rad``).

The pure-Python twin of ``src/cpp/src/formats/radioss.cpp``: both engines read
the same meshes.

``#RADIOSS STARTER``, ``/BEGIN`` (run name, input version, units), then
``/KEYWORD/option/id`` blocks up to ``/END``; ``#``/``$`` lines are comments and
``#include file`` inlines another file. Fields are 10 (integer) and 20 (real)
columns from format 51 on, 8 and 16 before; a line with a comma is split on
commas. Nodes and elements make the mesh; parts, subsets, groups and ``/SURF/SEG``
surfaces become regions.
"""

import math
import os

import numpy as np

from .. import _provenance
from .._common import num_nodes_per_cell, warn
from .._exceptions import ReadError, WriteError
from .._facets import FacetIndex, facet_nodes
from .._mesh import Mesh, topological_dimension
from .._node_order import node_order
from .._regions import Region
from ..lsdyna._cards import format_real_fit, format_real_short
from ..lsdyna._lsdyna import _collapse_solid, _expand_solid

__all__ = ["read", "write"]

_MAX_INCLUDE_DEPTH = 8

# keyword -> (group family, nodes, meshio++ type)
_ELEMENTS = {
    "BRICK": ("BRIC", 8, "hexahedron"),
    "PENTA6": ("BRIC", 6, "wedge"),
    "TETRA4": ("BRIC", 4, "tetra"),
    "TETRA10": ("BRIC", 10, "tetra10"),
    "BRIC20": ("BRIC", 20, "hexahedron20"),
    "SHELL": ("SHEL", 4, "quad"),
    "SH3N": ("SH3N", 3, "triangle"),
    "QUAD": ("QUAD", 4, "quad"),
    "TRIA": ("TRIA", 3, "triangle"),
    "BEAM": ("BEAM", 2, "line"),
    "TRUSS": ("TRUS", 2, "line"),
    "SPRING": ("SPRI", 2, "line"),
}

_GROUPS = {
    "GRNOD": "NODE",
    "GRBRIC": "BRIC",
    "GRSHEL": "SHEL",
    "GRSH3N": "SH3N",
    "GRQUAD": "QUAD",
    "GRTRIA": "TRIA",
    "GRBEAM": "BEAM",
    "GRTRUS": "TRUS",
    "GRSPRI": "SPRI",
}

_UNREAD = ("TSHELL", "TSH3N", "SHEL16", "SPHCEL", "RIVET", "XELEM")

_SURF_KINDS = (
    "BOX",
    "BOX2",
    "PART",
    "SUBSET",
    "MAT",
    "PROP",
    "GRBRIC",
    "GRSHEL",
    "GRSH3N",
    "GRTRIA",
    "SURF",
)


def _collect(path, depth, out):
    """The deck's lines with includes inlined and comments dropped; True at /END."""
    if depth > _MAX_INCLUDE_DEPTH:
        raise ReadError(f"Radioss: #include nested deeper than {_MAX_INCLUDE_DEPTH}")
    try:
        with open(path, "rb") as f:
            text = f.read().decode("latin-1")
    except OSError:
        raise ReadError(f"Radioss: cannot open {path}")
    label = os.path.basename(path)
    for number, line in enumerate(text.split("\n"), start=1):
        if line.endswith("\r"):
            line = line[:-1]
        if line[:9].upper().startswith("#INCLUDE"):
            name = line[8:].strip(" \t\r")
            if len(name) >= 2 and name[0] in "\"'" and name[-1] == name[0]:
                name = name[1:-1]
            name = name.replace("\\", "/")
            inc = (
                name
                if os.path.isabs(name)
                else os.path.join(os.path.dirname(path), name)
            )
            if not os.path.isfile(inc):
                warn(f"Radioss: include file '{name}' not found; skipped")
                continue
            _collect(inc, depth + 1, out)
            continue
        if depth > 0 and line.strip(" \t\r").upper().startswith("#ENDDATA"):
            return False
        # Blank lines stay: a blank /PART or group title is still its title line.
        if line and line[0] in "#$":
            continue
        if line.strip(" \t\r").upper() == "/END":
            if depth == 0:
                return True
            continue
        out.append((line, f"{label}:{number}"))
    return False


def _fields(line, width, count):
    if "," in line:
        return [f.strip(" \t\r") for f in line.split(",")]
    return [line[k * width : (k + 1) * width].strip(" \t\r") for k in range(count)]


def _int(text, where):
    text = text.strip()
    if not text:
        return 0
    try:
        return int(text)
    except ValueError:
        try:
            v = float(text.replace("D", "E").replace("d", "e"))
        except ValueError:
            raise ReadError(f"Radioss: bad integer field '{text}' ({where})")
        if v != int(v):
            raise ReadError(f"Radioss: bad integer field '{text}' ({where})")
        return int(v)


def _real(text, where):
    text = text.strip()
    if not text:
        return 0.0
    t = text.replace("D", "E").replace("d", "e")
    try:
        return float(t)
    except ValueError:
        # Fortran's 1.5-3 (no exponent letter).
        k = max(t.rfind("+"), t.rfind("-"))
        if k > 0 and t[k - 1] not in "eE":
            try:
                return float(t[:k] + "E" + t[k:])
            except ValueError:
                pass
        raise ReadError(f"Radioss: bad real field '{text}' ({where})")


_IMPERIAL = {"in": 0.0254, "ft": 0.3048, "yd": 0.9144, "mi": 1609.344}
_SI_PREFIXES = {
    "": 1.0,
    "y": 1e-24,
    "z": 1e-21,
    "a": 1e-18,
    "f": 1e-15,
    "p": 1e-12,
    "n": 1e-9,
    "mu": 1e-6,
    "m": 1e-3,
    "c": 1e-2,
    "d": 1e-1,
    "da": 1e1,
    "h": 1e2,
    "k": 1e3,
    "M": 1e6,
    "G": 1e9,
    "T": 1e12,
    "P": 1e15,
    "E": 1e18,
    "Z": 1e21,
    "Y": 1e24,
}


def _length_unit(unit):
    """A /BEGIN length unit in metres (an SI prefix and "m", a few imperial
    names, or a number); 0 when unknown."""
    if not unit:
        return 0.0
    if unit in _IMPERIAL:
        return _IMPERIAL[unit]
    if unit.endswith("m") and unit[:-1] in _SI_PREFIXES:
        return _SI_PREFIXES[unit[:-1]]
    try:
        return _real(unit, "")
    except ReadError:
        return 0.0


def _slice(line, at, width, comma_field):
    """A fixed-width slice of a line (the comma field when the line has commas)."""
    if "," in line:
        f = [x.strip(" \t\r") for x in line.split(",")]
        return f[comma_field] if comma_field < len(f) else ""
    return line[at : at + width].strip(" \t\r")


def _in_box(boxes, ident, x, node_point, skewed, missing, depth=0):
    """Whether point `x` is inside box `ident` (/BOX/BOX: its positive boxes,
    minus its negative ones; a skewed box contains nothing)."""
    b = boxes.get(ident)
    if b is None or depth > 16:
        missing.add(ident)
        return False
    kind = b["kind"]
    if kind == "BOX":
        inside = any(
            _in_box(boxes, c, x, node_point, skewed, missing, depth + 1)
            for c in b["children"]
            if c > 0
        )
        if inside and any(
            _in_box(boxes, -c, x, node_point, skewed, missing, depth + 1)
            for c in b["children"]
            if c < 0
        ):
            inside = False
        return inside
    p1 = node_point(b["n1"], b["p1"])
    if b["skew"] and kind == "RECTA":
        # Edges along the skew's axes: compare in the skew frame.
        axes = b.get("axes")
        if axes is None:
            skewed.add(ident)
            return False
        p2 = node_point(b["n2"], b["p2"])
        for a in axes:
            v = a[0] * x[0] + a[1] * x[1] + a[2] * x[2]
            u1 = a[0] * p1[0] + a[1] * p1[1] + a[2] * p1[2]
            u2 = a[0] * p2[0] + a[1] * p2[1] + a[2] * p2[2]
            if v < min(u1, u2) or v > max(u1, u2):
                return False
        return True
    if kind == "SPHER":
        return sum((x[d] - p1[d]) ** 2 for d in range(3)) <= 0.25 * b["diameter"] ** 2
    p2 = node_point(b["n2"], b["p2"])
    if kind == "RECTA":
        return all(min(p1[d], p2[d]) <= x[d] <= max(p1[d], p2[d]) for d in range(3))
    if kind == "CYLIN":
        a = [p2[d] - p1[d] for d in range(3)]
        q = [x[d] - p1[d] for d in range(3)]
        aa = sum(v * v for v in a)
        qa = sum(q[d] * a[d] for d in range(3))
        if aa == 0.0 or qa < 0.0 or qa > aa:
            return False
        r2 = sum((q[d] - qa / aa * a[d]) ** 2 for d in range(3))
        return r2 <= 0.25 * b["diameter"] ** 2
    return False


def _header_version(head):
    """The input version on a ``#RADIOSS STARTER`` line (41 in
    ``#RADIOSS STARTER      41BAR2V41B``), 0 if none."""
    at = head.find("#RADIOSS STARTER")
    if at < 0:
        return 0
    rest = head[at + 16 :].split("\n", 1)[0].lstrip(" \t")
    digits = ""
    for ch in rest[:4]:
        if not ch.isdigit():
            break
        digits += ch
    return int(digits) if digits else 0


def _engine_fields(path):
    """An engine deck's keywords and the numbers of their lines, as
    ``radioss:engine:<keyword>`` field data (an output request is empty)."""
    lines = []
    _collect(path, 0, lines)
    out = []
    for text, _ in lines:
        t = text.strip(" \t\r")
        if not t:
            continue
        if t[0] == "/":
            out.append(("radioss:engine:" + t[1:], []))
            continue
        if not out:
            continue
        for tok in t.split():
            try:
                out[-1][1].append(float(tok.replace("D", "E").replace("d", "e")))
            except ValueError:
                pass
    return out


def _add_engine_fields(field_data, path):
    for name, values in _engine_fields(path):
        if name in field_data:  # a repeated keyword: its numbers appended
            values = list(np.asarray(field_data[name]).ravel()) + values
        field_data[name] = np.array(values, dtype=np.float64)


def _cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def _unit(a):
    n = math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    return ([x / n for x in a] if n > 0.0 else list(a)), n > 0.0


def _skew_axes(skew):
    """The unit axes of a fixed skew (X, Z = X x Y, Y = Z x X), as the rows
    of a rotation into the skew frame; None when X and Y are parallel."""
    x, okx = _unit(skew["x"])
    z, okz = _unit(_cross(x, skew["y"]))
    if not okx or not okz:
        return None
    return [x, _cross(z, x), z]


def read(filename):
    path = os.fspath(filename)
    try:
        with open(path, "rb") as f:
            head = f.read(512).decode("latin-1").upper()
    except OSError:
        raise ReadError(f"Radioss: cannot open {path}")
    if "#RADIOSS ENGINE" in head:
        mesh = Mesh(np.zeros((0, 3)), [])
        _add_engine_fields(mesh.field_data, path)
        return mesh
    header_version = _header_version(head)
    lines = []
    _collect(path, 0, lines)

    version = 2019
    iw, rw = 10, 20
    # Before input version 5.1 (4.1, 4.4) there is no /BEGIN: the version is on
    # the #RADIOSS STARTER line, fields are 8 and 16 columns wide, and a title
    # is the keyword's last part rather than a line of its own.
    titles_in_path = 0 < header_version < 51
    if titles_in_path:
        version = header_version
        iw, rw = 8, 16
    node_ids, coords = [], []
    elements = []  # (family, id, part, type, nodes, where)
    parts = {}
    part_order = []
    groups = []  # (keyword, subtype, title, id, ids)
    # {id, title, segs, kind (SEG, PART, SUBSET, MAT, PROP, GRBRIC, GRSHEL,
    # GRSH3N, GRTRIA, SURF), mode (EXT, ALL, FREE or ""), ids}
    surfaces = []
    subsets = {}
    boxes = {}
    skews = {}
    analytic = []  # (name, values) of /SURF/PLANE and /SURF/ELLIPS
    ellipsoid_skew = {}
    length_scale = 1.0
    work_length = 0.0  # the work length unit in metres, 0 if /BEGIN names none
    # /UNIT/<id>: the local length units in metres, read first (a keyword can
    # name a unit defined further down).
    unit_length = {}
    for k in range(len(lines) - 2):
        t = lines[k][0].strip(" \t\r").upper()
        if not t.startswith("/UNIT/"):
            continue
        uid = t[6:].strip()
        if not uid.isdigit():
            continue
        f = _fields(lines[k + 2][0], 20, 3)
        length = _length_unit(f[1]) if len(f) > 1 else 0.0
        if length > 0.0:
            unit_length[int(uid)] = length
        else:
            warn(f"Radioss: /UNIT/{uid} has no length unit meshio++ knows; ignored")
    warned_units = set()
    skipped_keywords = set()
    zero_springs = linear_bric20 = 0

    n = len(lines)
    i = 0
    while i < n:
        text, where = lines[i]
        if not text or text[0] != "/":
            i += 1
            continue
        raw_path = [p.strip(" \t\r") for p in text.strip(" \t\r")[1:].split("/")]
        path_ = [p.upper() for p in raw_path]
        key = path_[0]
        body = i + 1
        end = body
        while end < n and (not lines[end][0] or lines[end][0][0] != "/"):
            end += 1

        # The keyword's id: its first integer field (`/BOX/RECTA/3/1` is box 3 in
        # unit system 1; `/SURF/PART/EXT/12` is surface 12).
        def last_id():
            for part in path_[1:]:
                if part and all(ch in "+-0123456789" for ch in part):
                    return _int(part, where)
            return _int(path_[-1], where) if len(path_) >= 2 else 0

        # The unit system a keyword names: the integer after its option id
        # (`/BOX/RECTA/3/1`), or its only integer when it has none (`/NODE/1`).
        def unit_of(option_id):
            ints = [
                _int(part, where)
                for part in path_[1:]
                if part and all(ch in "+-0123456789" for ch in part)
            ]
            at = 1 if option_id else 0
            return ints[at] if len(ints) > at else 0

        # Lengths of this keyword into the work units: its /UNIT's, else /BEGIN's.
        def scale_of(unit):
            if unit == 0:
                return length_scale
            if unit not in unit_length or work_length <= 0.0:
                if unit not in warned_units:
                    warned_units.add(unit)
                    why = (
                        "not defined"
                        if unit not in unit_length
                        else "used without /BEGIN units"
                    )
                    warn(
                        f"Radioss: unit system {unit} is {why}; its lengths are read in "
                        "the input units"
                    )
                return length_scale
            return unit_length[unit] / work_length

        # A keyword's title: its own line from input version 5.1 on; before,
        # the keyword's last part (`/PART/1/CUIVRE`), the data then starting on
        # the next line. Returns (title, first data line).
        def title_of(k):
            if titles_in_path:
                last = raw_path[-1]
                named = len(raw_path) >= 3 and not all(
                    ch in "+-0123456789" for ch in last
                )
                return (last if named else ""), k
            return (lines[k][0].strip(" \t\r") if k < end else ""), k + 1

        if key == "BEGIN":
            # A deck with /BEGIN has its titles on lines of their own, whatever
            # its input version.
            titles_in_path = False
            if body + 1 < end:
                f = _fields(lines[body + 1][0], 10, 2)
                if f and f[0]:
                    version = _int(f[0], lines[body + 1][1])
            iw, rw = (10, 20) if version >= 51 else (8, 16)
            # Input and work units (mass, length, time; 20 columns each): the
            # solver works in the work units, so lengths are converted.
            if body + 2 < end:
                unit_in = _fields(lines[body + 2][0], 20, 3)
                work = _fields(lines[body + 3][0], 20, 3) if body + 3 < end else []
                li = unit_in[1] if len(unit_in) > 1 else ""
                lw = work[1] if len(work) > 1 and work[1] else li
                fi, fw = _length_unit(li), _length_unit(lw)
                work_length = fw
                if li and (fi <= 0.0 or fw <= 0.0):
                    warn(
                        f"Radioss: unknown length unit '{li if fi <= 0.0 else lw}' in "
                        "/BEGIN; lengths are read as written"
                    )
                elif li:
                    length_scale = fi / fw
        elif key == "BOX" and len(path_) >= 3:
            kind = path_[1]
            b = {
                "kind": kind,
                "skew": 0,
                "n1": 0,
                "n2": 0,
                "p1": [0.0] * 3,
                "p2": [0.0] * 3,
                "diameter": 0.0,
                "children": [],
            }
            _, k = title_of(body)
            box_scale = scale_of(unit_of(True))

            def real3(line_no):
                if line_no >= end:
                    return [0.0] * 3
                f = _fields(lines[line_no][0], rw, 3)
                return [
                    (
                        _real(f[d], lines[line_no][1]) * box_scale
                        if d < len(f) and f[d]
                        else 0.0
                    )
                    for d in range(3)
                ]

            def int_at(line_no, field):
                if line_no >= end:
                    return 0
                t = _slice(lines[line_no][0], field * iw, iw, field)
                return _int(t, lines[line_no][1]) if t else 0

            def diameter(line_no):
                if line_no >= end:
                    return 0.0
                t = _slice(lines[line_no][0], 3 * iw, rw, 2 if kind == "SPHER" else 3)
                return _real(t, lines[line_no][1]) * box_scale if t else 0.0

            if kind == "RECTA":
                b["n1"], b["n2"], b["skew"] = int_at(k, 0), int_at(k, 1), int_at(k, 2)
                b["p1"], b["p2"] = real3(k + 1), real3(k + 2)
            elif kind == "CYLIN":
                b["n1"], b["n2"], b["diameter"] = (
                    int_at(k, 0),
                    int_at(k, 1),
                    diameter(k),
                )
                b["p1"], b["p2"] = real3(k + 1), real3(k + 2)
            elif kind == "SPHER":
                b["n1"], b["diameter"] = int_at(k, 0), diameter(k)
                b["p1"] = real3(k + 1)
            elif kind == "BOX":
                for kk in range(k, end):
                    for f in _fields(lines[kk][0], iw, 10):
                        if f:
                            b["children"].append(_int(f, lines[kk][1]))
            else:
                skipped_keywords.add(f"/BOX/{kind}")
            boxes[last_id()] = b
        elif key == "SKEW" and len(path_) >= 3 and path_[1] == "FIX":
            # origin (from 5.1), X axis, Y axis
            _, k = title_of(body)
            sk_scale = scale_of(unit_of(True))

            def vec(line_no, scale):
                if line_no >= end:
                    return [0.0] * 3
                f = _fields(lines[line_no][0], rw, 3)
                w = lines[line_no][1]
                return [
                    _real(f[d], w) * scale if d < len(f) and f[d] else 0.0
                    for d in range(3)
                ]

            origin = [0.0] * 3
            if not titles_in_path:  # 4.x skews have no origin line
                origin = vec(k, sk_scale)
                k += 1
            skews[last_id()] = {
                "origin": origin,
                "x": vec(k, 1.0),
                "y": vec(k + 1, 1.0),
            }
        elif key == "SURF" and len(path_) >= 3 and path_[1] in ("PLANE", "ELLIPS"):
            # Analytical surfaces: no segments, their definition as field data.
            sid = last_id()
            _, k = title_of(body)
            sc = scale_of(unit_of(True))

            def reals(line_no, count):
                if line_no >= end:
                    return [0.0] * count
                f = _fields(lines[line_no][0], rw, count)
                w = lines[line_no][1]
                return [
                    _real(f[d], w) * sc if d < len(f) and f[d] else 0.0
                    for d in range(count)
                ]

            if path_[1] == "PLANE":
                analytic.append(
                    (f"radioss:surf_plane:{sid}", reals(k, 3) + reals(k + 1, 3))
                )
            else:
                f = _fields(lines[k][0], iw, 2) if k < end else []
                skew = _int(f[0], lines[k][1]) if f and f[0] else 0
                degree = _int(f[1], lines[k][1]) if len(f) > 1 and f[1] else 2
                ellipsoid_skew[sid] = skew
                analytic.append(
                    (
                        f"radioss:surf_ellips:{sid}",
                        [float(max(degree, 2))] + reals(k + 1, 3) + reals(k + 2, 3),
                    )
                )
        elif key == "NODE":
            node_scale = scale_of(unit_of(False))
            for k in range(body, end):
                ln, w = lines[k]
                if "," in ln:
                    f = _fields(ln, iw, 4)
                else:
                    f = [ln[:iw].strip()] + [
                        ln[iw + rw * d : iw + rw * (d + 1)].strip() for d in range(3)
                    ]
                if not f or not f[0]:
                    continue
                node_ids.append(_int(f[0], w))
                coords.append(
                    [
                        _real(f[d], w) * node_scale if d < len(f) else 0.0
                        for d in (1, 2, 3)
                    ]
                )
        elif key in _ELEMENTS:
            family, count, ctype = _ELEMENTS[key]
            part = last_id()
            per_record = 11 if key == "TETRA10" else count + 1
            values = []
            rec_where = where
            for k in range(body, end):
                ln, w = lines[k]
                if not ln.strip(" \t\r"):
                    continue
                if key == "TETRA10":
                    take = 1 if not values else 10
                elif key == "BRIC20":
                    take = 9 if not values else (8 if len(values) == 9 else 4)
                else:
                    take = per_record
                f = _fields(ln, iw, take)
                if not values:
                    rec_where = w
                for j in range(take):
                    if len(values) < per_record:
                        values.append(_int(f[j], w) if j < len(f) and f[j] else 0)
                if len(values) < per_record and key in ("TETRA10", "BRIC20"):
                    continue
                nodes = values[1 : 1 + count]
                t = ctype
                keep = True
                if key == "BRICK":
                    t, nodes = _collapse_solid(nodes)
                elif key == "SHELL" and (nodes[3] == nodes[2] or nodes[3] == 0):
                    t, nodes = "triangle", nodes[:3]
                elif key == "BRIC20":
                    if 0 in nodes[8:]:
                        linear_bric20 += 1
                        t, nodes = "hexahedron", nodes[:8]
                    else:
                        order = node_order("radioss", "hexahedron20")
                        nodes = [nodes[j] for j in order.to_meshio]
                elif key == "SPRING" and nodes[1] == 0:
                    zero_springs += 1
                    keep = False
                if keep:
                    elements.append(
                        (
                            family,
                            values[0],
                            part,
                            t,
                            list(nodes),
                            rec_where,
                            key in ("TETRA4", "TETRA10"),
                        )
                    )
                values = []
            if values:
                raise ReadError(
                    f"Radioss: /{key} element {values[0]} is cut short ({rec_where})"
                )
        elif key == "PART":
            pid = last_id()
            prop, mat, subset = 0, 0, 0
            title, k = title_of(body)
            if k < end:
                f = _fields(lines[k][0], iw, 3)
                w = lines[k][1]
                prop = _int(f[0], w) if f and f[0] else 0
                mat = _int(f[1], w) if len(f) > 1 and f[1] else 0
                subset = _int(f[2], w) if len(f) > 2 and f[2] else 0
            if pid not in parts:
                parts[pid] = (title, prop, mat, subset)
                part_order.append(pid)
        elif key in _GROUPS and len(path_) >= 3:
            title, k = title_of(body)
            ids = []
            for kk in range(k, end):
                for f in _fields(lines[kk][0], iw, 10):
                    if f:
                        ids.append(_int(f, lines[kk][1]))
            groups.append((key, path_[1], title, last_id(), ids))
        elif key == "SURF" and len(path_) >= 3 and path_[1] == "SEG":
            title, k = title_of(body)
            segs = []
            for kk in range(k, end):
                f = _fields(lines[kk][0], iw, 5)
                w = lines[kk][1]
                seg = [
                    _int(f[j + 1], w) if j + 1 < len(f) and f[j + 1] else 0
                    for j in range(4)
                ]
                if seg[0] or seg[1] or seg[2]:
                    segs.append(seg)
            surfaces.append(
                {
                    "id": last_id(),
                    "title": title,
                    "segs": segs,
                    "kind": "SEG",
                    "mode": "",
                    "ids": [],
                }
            )
        elif key == "SURF" and len(path_) >= 3 and path_[1] in _SURF_KINDS:
            title, k = title_of(body)
            ids = []
            for kk in range(k, end):
                for f in _fields(lines[kk][0], iw, 10):
                    if f:
                        ids.append(_int(f, lines[kk][1]))
            surfaces.append(
                {
                    "id": last_id(),
                    "title": title,
                    "segs": [],
                    "kind": path_[1],
                    "mode": (
                        path_[2]
                        if len(path_) >= 4 and path_[2] in ("EXT", "ALL", "FREE")
                        else ""
                    ),
                    "ids": ids,
                }
            )
        elif key == "SURF":
            skipped_keywords.add("/SURF/" + (path_[1] if len(path_) > 1 else ""))
        elif key == "SUBSET":
            title, k = title_of(body)
            children = []
            for kk in range(k, end):
                for f in _fields(lines[kk][0], iw, 10):
                    if f:
                        children.append(_int(f, lines[kk][1]))
            subsets[last_id()] = (title, children)
        elif key.startswith("GR") or key in _UNREAD:
            skipped_keywords.add("/" + key)
        i = end
    if skipped_keywords:
        warn(f"Radioss: mesh keywords not read: {', '.join(sorted(skipped_keywords))}")
    if zero_springs:
        warn(f"Radioss: {zero_springs} /SPRING element(s) with a single node skipped")
    if linear_bric20:
        warn(
            f"Radioss: {linear_bric20} /BRIC20 element(s) with missing mid-edge nodes "
            "read as hexahedra"
        )

    for b in boxes.values():
        if b["skew"]:
            b["axes"] = _skew_axes(skews[b["skew"]]) if b["skew"] in skews else None

    node_index = {}
    for p, nid in enumerate(node_ids):
        if nid in node_index:
            raise ReadError(f"Radioss: node {nid} is defined twice")
        node_index[nid] = p
    points = np.array(coords, dtype=np.float64).reshape(-1, 3)

    # Tetrahedra come in either winding (every /TETRA4 of OpenRadioss's INT_25 QA
    # deck is inverted, gmsh writes them positive): an inverted one is mirrored.
    reoriented = 0
    for el in elements:
        if not el[6]:
            continue
        nodes = el[4]
        idx = [node_index.get(v) for v in nodes[:4]]
        if None in idx:
            continue
        p = points[idx]
        if np.linalg.det(np.array([p[1] - p[0], p[2] - p[0], p[3] - p[0]])) >= 0.0:
            continue
        nodes[1], nodes[2] = nodes[2], nodes[1]
        if len(nodes) == 10:  # edges 0-1 <-> 0-2 and 1-3 <-> 2-3
            nodes[4], nodes[6] = nodes[6], nodes[4]
            nodes[8], nodes[9] = nodes[9], nodes[8]
        reoriented += 1
    if reoriented:
        warn(f"Radioss: {reoriented} inverted tetrahedra reoriented")

    order_of_types = []
    by_type = {}
    for e, el in enumerate(elements):
        if el[3] not in by_type:
            by_type[el[3]] = []
            order_of_types.append(el[3])
        by_type[el[3]].append(e)
    owner = {}
    part_cells = {}
    cell_dim = []
    cell_family = []
    cell_part = []
    cell_conn = []
    cells = []
    part_blocks, prop_blocks, mat_blocks = [], [], []
    for t in order_of_types:
        members = by_type[t]
        k = num_nodes_per_cell[t]
        conn = np.empty((len(members), k), dtype=np.int64)
        pa = np.empty(len(members), dtype=np.int64)
        pr = np.empty(len(members), dtype=np.int64)
        ma = np.empty(len(members), dtype=np.int64)
        for r, e in enumerate(members):
            family, eid, part, _, nodes, w, _ = elements[e]
            for j in range(k):
                idx = node_index.get(nodes[j])
                if idx is None:
                    raise ReadError(
                        f"Radioss: element {eid} names undefined node {nodes[j]} ({w})"
                    )
                conn[r, j] = idx
            cell = len(cell_dim)
            fam = owner.setdefault(family, {})
            if eid in fam:
                raise ReadError(f"Radioss: element {eid} is defined twice ({w})")
            fam[eid] = cell
            p = parts.get(part)
            pa[r] = part
            pr[r] = p[1] if p else 0
            ma[r] = p[2] if p else 0
            part_cells.setdefault(part, []).append(cell)
            cell_dim.append(topological_dimension[t])
            cell_family.append(family)
            cell_part.append(part)
            cell_conn.append(conn[r].tolist())
        cells.append((t, conn))
        part_blocks.append(pa)
        prop_blocks.append(pr)
        mat_blocks.append(ma)

    mesh = Mesh(points, cells)
    mesh.field_data["radioss:version"] = np.array(version, dtype=np.int64)
    if length_scale != 1.0:
        mesh.field_data["radioss:length_scale"] = np.array(length_scale)
    if cells:
        mesh.cell_data["radioss:part"] = part_blocks
        mesh.cell_data["radioss:property"] = prop_blocks
        mesh.cell_data["radioss:material"] = mat_blocks

    regions = []
    seen = set()

    def add_region(name, kind, tag, keyword, entries):
        dim = -1
        if kind == "cell":
            dim = max((cell_dim[c] for c in entries), default=-1)
        elif kind == "side":
            dim = max((cell_dim[c] - 1 for c, _ in entries), default=-1)
        if (kind, name, tag) in seen:
            name += f" [{keyword}]"
        seen.add((kind, name, tag))
        arr = np.array(entries, dtype=np.int64)
        if kind == "side":
            arr = arr.reshape(-1, 2)
        regions.append(Region(name, kind, arr, dim, tag))

    for pid in part_cells:
        if pid not in parts:
            parts[pid] = ("", 0, 0, 0)
            part_order.append(pid)
    for pid in part_order:
        title = parts[pid][0]
        add_region(
            title or f"Part {pid}", "cell", pid, "PART", list(part_cells.get(pid, []))
        )

    def subset_closure(root):
        ids = set()
        stack = [root]
        while stack:
            c = stack.pop()
            if c in ids:
                continue
            ids.add(c)
            if c in subsets:
                stack.extend(subsets[c][1])
        return ids

    for sid in sorted(subsets):
        title, _ = subsets[sid]
        ids = subset_closure(sid)
        members = []
        for pid in part_order:
            if parts[pid][3] in ids:
                members.extend(part_cells.get(pid, []))
        add_region(title or f"Subset {sid}", "cell", sid, "SUBSET", members)

    group_of = {(g[0], g[3]): k for k, g in enumerate(groups)}
    skipped_subtypes = set()
    dropped = 0
    skewed_boxes, missing_boxes = set(), set()

    def node_point(node, fallback):
        p = node_index.get(node) if node else None
        return fallback if p is None else points[p].tolist()

    def resolve(g, visiting):
        nonlocal dropped
        out = set()
        if g in visiting:
            return out
        visiting.add(g)
        keyword, subtype, _, _, ids = groups[g]
        family = _GROUPS[keyword]
        nodes = family == "NODE"
        # GENE: `first last` id ranges; GEN_INCR: `first last step`.
        if subtype in ("GENE", "GEN_INCR"):
            w = 2 if subtype == "GENE" else 3
            spans = [ids[k : k + w] for k in range(0, len(ids) - w + 1, w)]

            def take(ident):
                for span in spans:
                    first, last = span[0], span[1]
                    step = span[2] if w == 3 else 1
                    if (
                        first <= ident <= last
                        and step > 0
                        and (ident - first) % step == 0
                    ):
                        return True
                return False

            source = node_index if nodes else owner.get(family, {})
            out = {entity for ident, entity in source.items() if take(ident)}
            visiting.discard(g)
            return out
        removed = set()
        for raw in ids:
            ident = abs(raw)
            hits = set()
            if subtype == family:
                idx = (
                    node_index.get(ident) if nodes else owner.get(family, {}).get(ident)
                )
                if idx is None:
                    dropped += 1
                else:
                    hits.add(idx)
            elif subtype in ("PART", "SUBSET"):
                # The parts' nodes for GRNOD, else their cells of the group's family.
                pids = [ident]
                if subtype == "SUBSET":
                    closure = subset_closure(ident)
                    pids = [pid for pid in part_order if parts[pid][3] in closure]
                for pid in pids:
                    for c in part_cells.get(pid, []):
                        if nodes:
                            hits.update(cell_conn[c])
                        elif cell_family[c] == family:
                            hits.add(c)
            elif subtype in ("BOX", "BOX2"):
                # Nodes inside; elements with all (BOX) or any (BOX2) node inside.
                def inside(p):
                    return _in_box(
                        boxes, ident, points[p], node_point, skewed_boxes, missing_boxes
                    )

                if nodes:
                    hits.update(p for p in range(len(points)) if inside(p))
                else:
                    for cell in owner.get(family, {}).values():
                        flags = [inside(p) for p in cell_conn[cell]]
                        if any(flags) if subtype == "BOX2" else all(flags):
                            hits.add(cell)
            elif nodes and subtype == "SURF":
                for surf in surfaces:
                    if surf["id"] == ident:
                        for seg in surf["segs"]:
                            for v in seg:
                                if v and v in node_index:
                                    hits.add(node_index[v])
            elif nodes and _GROUPS.get(subtype, "NODE") != "NODE":
                # A GRNOD of element groups: those elements' nodes.
                sub = group_of.get((subtype, ident))
                if sub is not None:
                    for c in resolve(sub, visiting):
                        hits.update(cell_conn[c])
            elif subtype == keyword:
                sub = group_of.get((keyword, ident))
                if sub is not None:
                    hits.update(resolve(sub, visiting))
            else:
                skipped_subtypes.add(f"/{keyword}/{subtype}")
                break
            (removed if raw < 0 else out).update(hits)
        out -= removed
        visiting.discard(g)
        return out

    for g, (keyword, _, title, gid, _) in enumerate(groups):
        entries = sorted(resolve(g, set()))
        add_region(
            title or f"{keyword}_{gid}",
            "point" if _GROUPS[keyword] == "NODE" else "cell",
            gid,
            keyword,
            entries,
        )
    for t in sorted(skipped_subtypes):
        warn(
            f"Radioss: {t} groups are not resolved (only entity ids, parts and other "
            "groups are); their regions are left empty or partial"
        )

    if surfaces:
        facets = FacetIndex(mesh, surface_self=True)
        surface_of = {}
        for k, surf in enumerate(surfaces):
            surface_of.setdefault(surf["id"], k)
        done = {}
        model_faces = None

        def solid_faces(cell):
            out = []
            f = 0
            while True:
                hit = facet_nodes(mesh, cell, f)
                if hit is None:
                    return out
                ftype, fnodes = hit
                corners = 3 if ftype.startswith("triangle") else 4
                out.append((f, tuple(sorted(fnodes[:corners]))))
                f += 1

        def side_set(index, depth):
            nonlocal dropped, model_faces
            if index in done:
                return done[index]
            surf = surfaces[index]
            out = set()
            kind = surf["kind"]
            if kind == "SURF":
                minus = set()
                for raw in surf["ids"]:
                    sub = surface_of.get(abs(raw))
                    if sub is None or depth > 16:
                        dropped += 1
                        continue
                    (minus if raw < 0 else out).update(side_set(sub, depth + 1))
                out -= minus
            elif kind in ("BOX", "BOX2"):
                # Shell faces with all (BOX) or any (BOX2) node in the box; with
                # EXT the model's external solid faces, with ALL every solid
                # face, likewise.
                box = surf["ids"][0] if surf["ids"] else 0

                def inside(pts):
                    n = sum(
                        1
                        for p in pts
                        if _in_box(
                            boxes,
                            box,
                            points[p],
                            node_point,
                            skewed_boxes,
                            missing_boxes,
                        )
                    )
                    return n > 0 if kind == "BOX2" else n == len(pts)

                mode = surf["mode"]
                if mode == "EXT" and model_faces is None:
                    model_faces = {}
                    for c in range(len(cell_dim)):
                        if cell_dim[c] == 3:
                            for _, key in solid_faces(c):
                                model_faces[key] = model_faces.get(key, 0) + 1
                for c in range(len(cell_conn)):
                    if cell_family[c] in ("SHEL", "SH3N", "TRIA"):
                        if inside(cell_conn[c]):
                            out.add((c, 0))
                    elif cell_dim[c] == 3 and mode:
                        for f, key in solid_faces(c):
                            if (mode == "ALL" or model_faces[key] == 1) and inside(key):
                                out.add((c, f))
            elif kind == "SEG":
                for seg in surf["segs"]:
                    # 2 nodes in a 2-D analysis, 3 for a triangle (n4 blank or n3).
                    if seg[2] == 0:
                        count = 2
                    else:
                        count = 3 if seg[3] == 0 or seg[3] == seg[2] else 4
                    idx = [node_index.get(v) for v in seg[:count]]
                    hit = facets.find(idx) if None not in idx else None
                    if hit is None:
                        dropped += 1
                        continue
                    out.add(tuple(hit.first))
            else:
                # The cells the surface draws from.
                cells = set()
                if kind.startswith("GR"):
                    for raw in surf["ids"]:
                        sub = group_of.get((kind, abs(raw)))
                        if sub is None:
                            dropped += 1
                            continue
                        cells.update(resolve(sub, set()))
                else:
                    for raw in surf["ids"]:
                        ident = abs(raw)
                        closure = subset_closure(ident) if kind == "SUBSET" else set()
                        for c, pid in enumerate(cell_part):
                            part = parts[pid]
                            if kind == "PART":
                                pick = pid == ident
                            elif kind == "SUBSET":
                                pick = part[3] in closure
                            elif kind == "MAT":
                                pick = part[2] == ident
                            else:
                                pick = part[1] == ident
                            if pick:
                                cells.add(c)
                # Shells: their own face. Solids: with EXT the faces no other
                # chosen solid shares, with FREE those no solid of the model
                # shares, with ALL every face (GRBRIC without a mode: EXT).
                mode = "EXT" if kind == "GRBRIC" and not surf["mode"] else surf["mode"]
                chosen = {}
                solids = []
                for c in sorted(cells):
                    if cell_family[c] in ("SHEL", "SH3N", "TRIA"):
                        out.add((c, 0))
                    elif cell_dim[c] == 3 and mode:
                        faces = solid_faces(c)
                        solids.append((c, faces))
                        for _, key in faces:
                            chosen[key] = chosen.get(key, 0) + 1
                if mode == "FREE" and model_faces is None:
                    model_faces = {}
                    for c in range(len(cell_dim)):
                        if cell_dim[c] == 3:
                            for _, key in solid_faces(c):
                                model_faces[key] = model_faces.get(key, 0) + 1
                for c, faces in solids:
                    for f, key in faces:
                        if mode == "FREE":
                            shared = model_faces[key]
                        elif mode == "EXT":
                            shared = chosen[key]
                        else:
                            shared = 1
                        if shared == 1:
                            out.add((c, f))
            done[index] = out
            return out

        for index, surf in enumerate(surfaces):
            entries = sorted(side_set(index, 0))
            add_region(
                surf["title"] or f"SURF_{surf['id']}",
                "side",
                surf["id"],
                "SURF",
                entries,
            )
    for name, values in analytic:
        # An ellipsoid's orientation: its skew's axes (the identity without one).
        sid = int(name.rsplit(":", 1)[1])
        if name.startswith("radioss:surf_ellips:"):
            axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            skew = ellipsoid_skew.get(sid, 0)
            if skew:
                got = _skew_axes(skews[skew]) if skew in skews else None
                if got is None:
                    warn(
                        f"Radioss: /SURF/ELLIPS {sid} names a skew that is not a /SKEW/FIX"
                    )
                else:
                    axes = got
            values = values + [x for a in axes for x in a]
        mesh.field_data[name] = np.array(values, dtype=np.float64)
    for t in sorted(skewed_boxes):
        warn(
            f"Radioss: /BOX {t} names a skew that is not a /SKEW/FIX; it contains nothing"
        )
    for t in sorted(missing_boxes):
        warn(f"Radioss: /BOX {t} is not defined; it contains nothing")
    if dropped:
        warn(
            f"Radioss: {dropped} group or surface entries name undefined ids or no cell "
            "facet and were dropped"
        )
    mesh.regions = regions
    # A starter deck `<run>_0000.rad`: its engine deck `<run>_0001.rad`'s controls.
    stem, ext = os.path.splitext(path)
    if stem.endswith("_0000"):
        engine = stem[:-5] + "_0001" + ext
        if os.path.isfile(engine):
            _add_engine_fields(mesh.field_data, engine)
    return mesh


# ----------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------

# meshio++ type -> the starter card that writes it. A pyramid and a wedge are
# degenerate /BRICKs (`1 2 3 4 5 5 5 5`, `1 3 6 4 2 2 5 5`), which the reader
# collapses back and which need no formulation of their own, unlike /PENTA6.
_WRITE_CARD = {
    "hexahedron": "BRICK",
    "pyramid": "BRICK",
    "wedge": "BRICK",
    "tetra": "TETRA4",
    "tetra10": "TETRA10",
    "hexahedron20": "BRIC20",
    "quad": "SHELL",
    "triangle": "SH3N",
    "line": "TRUSS",
}
# A group family's /GR keyword; /GRNOD is the point one.
_GROUP_KEYWORD = {
    "BRIC": "GRBRIC",
    "SHEL": "GRSHEL",
    "SH3N": "GRSH3N",
    "TRUS": "GRTRUS",
}
# The stub property type a card's part needs.
_STUB_KIND = {
    "BRICK": "SOLID",
    "TETRA4": "SOLID",
    "TETRA10": "SOLID",
    "BRIC20": "SOLID",
    "SHELL": "SHELL",
    "SH3N": "SHELL",
    "TRUSS": "TRUSS",
}


def _stub_need(cards):
    """The stub property cards ``cards`` can share, ``(type, Isolid)``, or None
    when they need different types. /BRIC20 needs Isolid 16, the other solids
    the default 0; a part mixing /BRIC20 and /BRICK gets 0, which the starter
    turns to 16 for the /BRIC20 with a warning."""
    kinds = {_STUB_KIND[c] for c in cards}
    if len(kinds) != 1:
        return None
    bricks = {c for c in cards if c in ("BRICK", "BRIC20")}
    return kinds.pop(), 16 if bricks == {"BRIC20"} else 0


_WRITE_VERSION = 2022
_WRITE_UNITS = ("kg", "m", "s")
_TITLE_WIDTH = 100


def _i10(*values):
    return "".join(f"{int(v):10d}" for v in values)


def _real20(value):
    """A real in an F20 field: exact when its shortest spelling fits, else
    the 20 columns' worth of digits."""
    text = format_real_short(value)
    return text if len(text) <= 20 else format_real_fit(value, 20)


def _f20(*values):
    return "".join(f"{_real20(v):>20}" for v in values)


def _title(name):
    """A title line: one line, at most 100 characters, never read as a comment
    or a keyword (a leading `#`, `$` or `/` gets a space in front)."""
    t = str(name).replace("\r", " ").replace("\n", " ").strip(" \t\x0b\x0c")
    if t[:1] in ("#", "$", "/"):
        t = " " + t
    return t[:_TITLE_WIDTH]


def _ids_lines(ids):
    """Ids, ten I10 fields per line."""
    return [_i10(*ids[k : k + 10]) for k in range(0, len(ids), 10)]


class _Ids:
    """One id namespace: a tag when positive and free, else the next free id."""

    def __init__(self):
        self.used = set()
        self.next = 1

    def take(self, tag=0):
        tag = int(tag)
        if tag > 0 and tag not in self.used:
            self.used.add(tag)
            return tag
        while self.next in self.used:
            self.next += 1
        self.used.add(self.next)
        return self.next


def write(filename, mesh, stubs=False):
    """Write ``mesh`` as an OpenRadioss starter deck (see ``radioss/__init__.py``)."""
    path = os.fspath(filename)
    points = np.asarray(mesh.points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] > 3:
        raise WriteError("Radioss writer: points must have 1 to 3 coordinates")
    npts = len(points)
    if npts == 0:
        raise WriteError(
            "Radioss writer: a starter deck needs nodes; the mesh has none"
        )

    # Cells: which are written, as what, with which element id.
    cell_info = []  # per global cell: (block, row, card) or None when dropped
    dropped_types = set()
    for b, block in enumerate(mesh.cells):
        card = _WRITE_CARD.get(block.type)
        ragged = isinstance(block.data, list)
        if card is None or ragged:
            dropped_types.add(block.type)
        for r in range(len(block.data)):
            cell_info.append(None if card is None or ragged else (b, r, card))
    for t in sorted(dropped_types):
        warn(f"Radioss writer: '{t}' cells have no starter element card; dropped")
        _provenance.note(
            "cells-dropped", f"Radioss has no element card for '{t}' cells"
        )
    ncells = len(cell_info)

    def family(c):
        return _ELEMENTS[cell_info[c][2]][0]

    # The C++ core keeps regions in (kind, name, dim, tag) order; walk them in
    # that order so both engines write the same bytes.
    kind_order = {"point": 0, "cell": 1, "side": 2}
    regions = sorted(
        getattr(mesh, "regions", None) or [],
        key=lambda r: (kind_order.get(r.kind, 3), r.name.encode("utf-8"), r.dim, r.tag),
    )
    used_regions = set()

    # Parts: radioss:part when the mesh carries it, else cell regions of one
    # family, else one part per block.
    part_ids = _Ids()
    cell_part = [0] * ncells
    parts = []  # [pid, title, prop, mat]

    def cell_array(name):
        if name not in mesh.cell_data or len(mesh.cell_data[name]) != len(mesh.cells):
            return None
        return [np.asarray(a).reshape(-1) for a in mesh.cell_data[name]]

    radioss_part = cell_array("radioss:part")
    if radioss_part is not None:
        # Part ids come from the data; an id of 0 or less gets a new one.
        raw = [0] * ncells
        g = 0
        for b, block in enumerate(mesh.cells):
            for r in range(len(block.data)):
                raw[g] = int(radioss_part[b][r])
                if cell_info[g] is not None and raw[g] > 0:
                    part_ids.used.add(raw[g])
                g += 1
        renumbered = {}
        for c in range(ncells):
            if cell_info[c] is None:
                continue
            pid = raw[c]
            if pid <= 0:
                if pid not in renumbered:
                    renumbered[pid] = part_ids.take()
                pid = renumbered[pid]
            cell_part[c] = pid
        if renumbered:
            warn("Radioss writer: cells with a part id of 0 or less put in new parts")
        prop = cell_array("radioss:property")
        mat = cell_array("radioss:material")
        members = {}
        g = 0
        first = {}
        for b, block in enumerate(mesh.cells):
            for r in range(len(block.data)):
                if cell_info[g] is not None:
                    pid = cell_part[g]
                    members.setdefault(pid, []).append(g)
                    first.setdefault(pid, (b, r))
                g += 1
        for pid, cells in members.items():
            title = ""
            for k, reg in enumerate(regions):
                if (
                    k not in used_regions
                    and reg.kind == "cell"
                    and int(reg.tag) == pid
                    and sorted(np.asarray(reg.entries).reshape(-1).tolist()) == cells
                ):
                    title = reg.name
                    used_regions.add(k)
                    break
            b, r = first[pid]
            pr = int(prop[b][r]) if prop is not None else pid
            ma = int(mat[b][r]) if mat is not None else 1
            parts.append([pid, title, pr, ma])
    else:
        claimed = [False] * ncells
        for k, reg in enumerate(regions):
            if reg.kind != "cell":
                continue
            cells = sorted(set(int(c) for c in np.asarray(reg.entries).reshape(-1)))
            if not cells or any(
                not 0 <= c < ncells or cell_info[c] is None or claimed[c] for c in cells
            ):
                continue
            # One stub property must be able to serve a part.
            if _stub_need({cell_info[c][2] for c in cells}) is None:
                continue
            pid = part_ids.take(reg.tag)
            for c in cells:
                claimed[c] = True
                cell_part[c] = pid
            used_regions.add(k)
            parts.append([pid, reg.name, pid, 1])
        g = 0
        for block in mesh.cells:
            pid = None
            for _ in range(len(block.data)):
                if cell_info[g] is not None and not claimed[g]:
                    if pid is None:
                        pid = part_ids.take()
                        parts.append([pid, block.type, pid, 1])
                    cell_part[g] = pid
                g += 1

    # Subsets: a remaining cell region that is exactly a union of whole parts.
    part_cells = {}
    for c in range(ncells):
        if cell_info[c] is not None:
            part_cells.setdefault(cell_part[c], set()).add(c)
    part_subset = {}
    subset_ids = _Ids()
    subsets = []  # (sid, title, children)
    subset_parts = {}
    # Smallest first, so a subset nested in another is its child.
    candidates = []
    for k, reg in enumerate(regions):
        if k in used_regions or reg.kind != "cell":
            continue
        cells = set(int(c) for c in np.asarray(reg.entries).reshape(-1))
        if cells and all(0 <= c < ncells and cell_info[c] is not None for c in cells):
            candidates.append((len(cells), k, cells))
    for _, k, cells in sorted(candidates):
        reg = regions[k]
        inside = [
            p[0] for p in parts if part_cells.get(p[0]) and part_cells[p[0]] <= cells
        ]
        if set().union(*(part_cells[p] for p in inside)) != cells:
            continue
        children = sorted({part_subset[p] for p in inside if p in part_subset})
        if any(not subset_parts[s] <= set(inside) for s in children):
            continue
        sid = subset_ids.take(reg.tag)
        for p in inside:
            part_subset.setdefault(p, sid)
        subset_parts[sid] = set(inside)
        subsets.append((sid, reg.name, children))
        used_regions.add(k)

    # Element ids, in the order the cards write them.
    elem_id = [0] * ncells
    next_elem = 1
    cards = []  # (card, pid, [global cells])
    g = 0
    for b, block in enumerate(mesh.cells):
        run = None
        for _ in range(len(block.data)):
            if cell_info[g] is not None:
                card = cell_info[g][2]
                if run is None or run[1] != cell_part[g]:
                    run = (card, cell_part[g], [])
                    cards.append(run)
                run[2].append(g)
                elem_id[g] = next_elem
                next_elem += 1
            g += 1

    # Groups, then surfaces.
    group_ids = {}
    groups = []  # (keyword, subtype, gid, title, ids)
    dropped_sides = 0
    other_regions = 0
    surf_ids = _Ids()
    surfaces = []  # (sid, title, segs)
    for k, reg in enumerate(regions):
        if k in used_regions:
            continue
        entries = np.asarray(reg.entries)
        if reg.kind == "point":
            ids = sorted({int(p) + 1 for p in entries.reshape(-1) if 0 <= p < npts})
            gid = group_ids.setdefault("GRNOD", _Ids()).take(reg.tag)
            groups.append(("GRNOD", "NODE", gid, reg.name, ids))
        elif reg.kind == "cell":
            by_family = {}
            for c in sorted(set(int(c) for c in entries.reshape(-1))):
                if 0 <= c < ncells and cell_info[c] is not None:
                    by_family.setdefault(family(c), []).append(elem_id[c])
            if not by_family:
                by_family = {"BRIC": []}
            for fam, ids in by_family.items():
                keyword = _GROUP_KEYWORD[fam]
                gid = group_ids.setdefault(keyword, _Ids()).take(reg.tag)
                groups.append((keyword, fam, gid, reg.name, ids))
        elif reg.kind == "side":
            segs = []
            # Sorted and without repeats, as the C++ core keeps a region.
            for c, f in sorted(set(map(tuple, entries.reshape(-1, 2).tolist()))):
                if not (0 <= c < ncells) or cell_info[c] is None:
                    dropped_sides += 1
                    continue
                b, r, card = cell_info[c]
                fam = family(c)
                if fam in ("SHEL", "SH3N"):
                    nodes = [int(v) for v in np.asarray(mesh.cells[b].data[r])]
                elif fam == "BRIC":
                    hit = facet_nodes(mesh, c, int(f))
                    if hit is None:
                        dropped_sides += 1
                        continue
                    ftype, fnodes = hit
                    nodes = fnodes[: 3 if ftype.startswith("triangle") else 4]
                else:
                    dropped_sides += 1
                    continue
                seg = [v + 1 for v in nodes]
                if len(seg) == 3:
                    seg.append(seg[2])
                segs.append(seg)
            sid = surf_ids.take(reg.tag)
            surfaces.append((sid, reg.name, segs))
        else:
            other_regions += 1
    if dropped_sides:
        warn(
            f"Radioss writer: {dropped_sides} side region entries on cells with no "
            "face segment dropped"
        )
        _provenance.note(
            "regions-dropped", "a /SURF/SEG segment is a solid face or a shell"
        )
    if other_regions:
        warn(f"Radioss writer: {other_regions} region(s) of unknown kind dropped")

    # Data: only the part, property, material and deck field data are written.
    planes = []
    dropped_data = len(mesh.point_data) + sum(
        1
        for name in mesh.cell_data
        if name not in ("radioss:part", "radioss:property", "radioss:material")
    )
    for name, value in mesh.field_data.items():
        if name in ("radioss:version", "radioss:length_scale"):
            continue
        if name.startswith("radioss:surf_plane:"):
            v = np.asarray(value, dtype=np.float64).reshape(-1)
            ident = name.rsplit(":", 1)[1]
            if len(v) == 6 and ident.isdigit():
                planes.append((int(ident), v.tolist()))
                continue
        dropped_data += 1
    if dropped_data:
        warn(
            "Radioss writer: a starter deck holds no data arrays; point, cell and "
            "field data other than the radioss: parts and analytical surfaces dropped"
        )
        _provenance.note("data-dropped", "a Radioss starter deck holds no data arrays")

    version = _WRITE_VERSION
    if "radioss:version" in mesh.field_data:
        v = int(np.asarray(mesh.field_data["radioss:version"]).reshape(-1)[0])
        if 2017 <= v <= 2026:
            version = v

    run = os.path.splitext(os.path.basename(path))[0]
    if run.endswith("_0000"):
        run = run[:-5]
    out = ["#RADIOSS STARTER"]
    out += ["# " + line for line in _provenance.lines(_provenance.SlotTier.BLOCK)]
    # The same input and work units, so the starter converts nothing and the
    # reader reads the coordinates as written.
    units = "".join(f"{u:>20}" for u in _WRITE_UNITS)
    out += ["/BEGIN", _title(run or "meshio")[:80], _i10(version, 0), units, units]

    if stubs:
        # One stub property per property id, as its parts' cards need it
        # (_STUB_NEED). A part whose property already serves a part that needs
        # another stub gets a property of its own; material 0 becomes 1.
        prop_ids = _Ids()
        prop_ids.used.update(p[2] for p in parts if p[2] > 0)
        props = {}  # prop id -> (kind, isolid)
        mats = []
        for p in parts:
            need = _stub_need({cell_info[c][2] for c in part_cells.get(p[0], ())})
            if need is None:
                raise WriteError(
                    f"Radioss writer: part {p[0]} holds solid, shell or truss elements "
                    "together, which no one stub property serves"
                )
            if p[2] <= 0:
                p[2] = prop_ids.take(p[0])
            if props.setdefault(p[2], need) != need:
                p[2] = prop_ids.take()
                props[p[2]] = need
            if p[3] <= 0:
                p[3] = 1
            if p[3] not in mats:
                mats.append(p[3])
        for mid in sorted(mats):
            out += [
                f"/MAT/LAW1/{mid}",
                f"meshio++ stub material {mid}",
                _f20(1.0, 0.0),
                _f20(1.0, 0.3),
            ]
        for pr in sorted(props):
            kind, isolid = props[pr]
            out += [f"/PROP/{kind}/{pr}", f"meshio++ stub property {pr}"]
            if kind == "SOLID":
                out += [
                    _i10(isolid, 0)
                    + " " * 10
                    + _i10(0)
                    + " " * 10
                    + _i10(0, 0, 0)
                    + _f20(0.0),
                    _f20(0.0, 0.0, 0.0, 0.0, 0.0),
                    _f20(0.0) + _i10(0, 0),
                ]
            elif kind == "SHELL":
                out += [
                    _i10(0, 0, 0),
                    _f20(0.0, 0.0, 0.0, 0.0, 0.0),
                    _i10(1, 0) + _f20(1.0),
                ]
            else:
                out += [_f20(1.0, 0.0)]

    for pid, title, pr, ma in parts:
        out += [f"/PART/{pid}", _title(title), _i10(pr, ma, part_subset.get(pid, 0))]
    for sid, title, children in subsets:
        out += [f"/SUBSET/{sid}", _title(title)] + _ids_lines(children)

    out.append("/NODE")
    for p in range(npts):
        xyz = [float(points[p, d]) if d < points.shape[1] else 0.0 for d in range(3)]
        out.append(_i10(p + 1) + _f20(*xyz))

    bric20 = node_order("radioss", "hexahedron20")
    for card, pid, cells in cards:
        out.append(f"/{card}/{pid}")
        for c in cells:
            b, r, _ = cell_info[c]
            nodes = [int(v) + 1 for v in np.asarray(mesh.cells[b].data[r])]
            eid = elem_id[c]
            if card == "BRICK":
                nodes = _expand_solid(mesh.cells[b].type, nodes)
            if card == "TETRA10":
                out += [_i10(eid), _i10(*nodes)]
            elif card == "BRIC20":
                file_nodes = [nodes[j] for j in bric20.from_meshio]
                out += [
                    _i10(eid, *file_nodes[:8]),
                    _i10(*file_nodes[8:16]),
                    _i10(*file_nodes[16:]),
                ]
            else:
                out.append(_i10(eid, *nodes))

    for keyword, subtype, gid, title, ids in groups:
        out += [f"/{keyword}/{subtype}/{gid}", _title(title)] + _ids_lines(ids)
    for sid, title, segs in surfaces:
        out += [f"/SURF/SEG/{sid}", _title(title)]
        out += [_i10(k + 1, *seg) for k, seg in enumerate(segs)]
    for sid, v in sorted(planes):
        out += [
            f"/SURF/PLANE/{surf_ids.take(sid)}",
            f"PLANE_{sid}",
            _f20(*v[:3]),
            _f20(*v[3:]),
        ]
    out.append("/END")
    # UTF-8, as the C++ writer's strings are (a Radioss deck is ASCII).
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(out) + "\n")
