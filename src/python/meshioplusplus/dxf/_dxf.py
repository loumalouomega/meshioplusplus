"""DXF R2000 (AC1015) drawings of line meshes and sketches (roadmap 7.6.4).

The writer emits an ASCII DXF R2000 file: the HEADER (``$ACADVER``,
``$HANDSEED``, ``$INSUNITS``), the TABLES and BLOCKS a reader expects (the
model and paper space block records), the ENTITIES, and an OBJECTS root
dictionary. A mesh carrying a sketch (``fit_curves``' ``sketch:*`` field data)
is written as its pieces -- ``LINE``, ``ARC`` and ``CIRCLE`` in world
coordinates, the arcs and circles in the object coordinate system of their
sketch plane's normal -- and any other line mesh as one ``LINE`` per segment.
The reader takes ``LINE``, ``ARC``, ``CIRCLE``, ``LWPOLYLINE`` and ``SPLINE``
(arcs, circles, bulges and splines sampled) and refuses every other entity by
name. The C++ core writes the same bytes.
"""

from __future__ import annotations

import math

import numpy as np

from .. import _provenance
from .._exceptions import ReadError, WriteError
from .._mesh import CellBlock, Mesh
from .._units import LENGTH_UNIT_KEY

_UNITS = (
    (1, 0.0254),
    (2, 0.3048),
    (3, 1609.344),
    (4, 0.001),
    (5, 0.01),
    (6, 1.0),
    (7, 1000.0),
    (10, 0.9144),
    (12, 1e-9),
    (13, 1e-6),
    (14, 0.1),
)
_READ = ("LINE", "ARC", "CIRCLE", "LWPOLYLINE", "SPLINE")
_CIRCLE_SEGMENTS = 64


def _num(v) -> str:
    v = float(v)
    if v == 0.0:
        v = 0.0  # no "-0"
    return format(v, ".17g")


def _insunits(mesh) -> int:
    fd = getattr(mesh, "field_data", None) or {}
    value = fd.get(LENGTH_UNIT_KEY)
    if value is None:
        return 0
    flat = np.asarray(value, dtype=float).ravel()
    if flat.size != 1:
        return 0
    for code, metres in _UNITS:
        if abs(flat[0] - metres) <= 1e-9 * metres:
            return code
    return 0


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _unit(a):
    n = math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    return (a[0] / n, a[1] / n, a[2] / n)


def _ocs(n):
    """The DXF arbitrary-axis algorithm: the object coordinate system's x and
    y axes for the extrusion direction ``n`` (scalar arithmetic, as the core)."""
    n = _unit(tuple(float(v) for v in n))
    if abs(n[0]) < 1.0 / 64.0 and abs(n[1]) < 1.0 / 64.0:
        ax = _cross((0.0, 1.0, 0.0), n)
    else:
        ax = _cross((0.0, 0.0, 1.0), n)
    ax = _unit(ax)
    ay = _cross(n, ax)
    return ax, ay, n


class _Out:
    def __init__(self):
        self.lines = []
        self.handle = 0x10

    def pair(self, code, value):
        self.lines.append(f"{code:>3}")
        self.lines.append(str(value))

    def next(self) -> str:
        h = format(self.handle, "X")
        self.handle += 1
        return h


def _pieces(mesh):
    """The sketch pieces in world coordinates, or ``None`` without a sketch."""
    fd = mesh.field_data
    need = ("sketch:type", "sketch:params", "sketch:frame", "sketch:frames")
    if not all(k in fd for k in need):
        return None
    types = np.asarray(fd["sketch:type"], dtype=np.int64).ravel()
    params = np.asarray(fd["sketch:params"], dtype=np.float64).reshape(-1, 10)
    frame_of = np.asarray(fd["sketch:frame"], dtype=np.int64).ravel()
    frames = np.asarray(fd["sketch:frames"], dtype=np.float64).reshape(-1, 4, 4)
    out = []
    for t, p, f in zip(types, params, frame_of):
        m = frames[int(f)]
        e1 = tuple(float(m[r, 0]) for r in range(3))
        e2 = tuple(float(m[r, 1]) for r in range(3))
        n = tuple(float(m[r, 2]) for r in range(3))
        o = tuple(float(m[r, 3]) for r in range(3))

        def world(u, v, e1=e1, e2=e2, o=o):
            u, v = float(u), float(v)
            return tuple(o[r] + u * e1[r] + v * e2[r] for r in range(3))

        out.append(
            (int(t), world(p[0], p[1]), world(p[2], p[3]), world(p[4], p[5]), p, n)
        )
    return out


def _entities(out, mesh):
    pieces = _pieces(mesh)
    if pieces is not None:
        for t, start, end, centre, p, n in pieces:
            if t == 1:
                out.pair(0, "LINE")
                out.pair(5, out.next())
                out.pair(330, "1F")
                out.pair(100, "AcDbEntity")
                out.pair(8, "0")
                out.pair(100, "AcDbLine")
                for k, v in enumerate(start):
                    out.pair(10 + 10 * k, _num(v))
                for k, v in enumerate(end):
                    out.pair(11 + 10 * k, _num(v))
                continue
            ax, ay, nn = _ocs(n)
            c = (_dot(centre, ax), _dot(centre, ay), _dot(centre, nn))
            out.pair(0, "CIRCLE" if t == 3 else "ARC")
            out.pair(5, out.next())
            out.pair(330, "1F")
            out.pair(100, "AcDbEntity")
            out.pair(8, "0")
            out.pair(100, "AcDbCircle")
            for k, v in enumerate(c):
                out.pair(10 + 10 * k, _num(v))
            out.pair(40, _num(p[6]))
            for k, v in enumerate(nn):
                out.pair(210 + 10 * k, _num(v))
            if t == 2:
                ds = tuple(start[r] - centre[r] for r in range(3))
                de = tuple(end[r] - centre[r] for r in range(3))
                a0 = math.atan2(_dot(ds, ay), _dot(ds, ax)) * 180.0 / math.pi
                a1 = math.atan2(_dot(de, ay), _dot(de, ax)) * 180.0 / math.pi
                if p[7] < 0.0:  # clockwise: DXF arcs run counter-clockwise
                    a0, a1 = a1, a0
                out.pair(100, "AcDbArc")
                out.pair(50, _num(a0))
                out.pair(51, _num(a1))
        return
    pts = np.asarray(mesh.points, dtype=np.float64)
    if pts.ndim == 2 and pts.shape[1] < 3:
        pts = np.hstack([pts, np.zeros((len(pts), 3 - pts.shape[1]))])
    for block in mesh.cells:
        if block.type == "vertex":
            continue
        if block.type != "line":
            raise WriteError(
                f"meshio++: dxf: only line cells or a sketch are written, not {block.type} "
                "(slice or sections give lines; fit_curves a sketch)"
            )
        for a, b in np.asarray(block.data, dtype=np.int64):
            out.pair(0, "LINE")
            out.pair(5, out.next())
            out.pair(330, "1F")
            out.pair(100, "AcDbEntity")
            out.pair(8, "0")
            out.pair(100, "AcDbLine")
            for k in range(3):
                out.pair(10 + 10 * k, _num(pts[a, k]))
            for k in range(3):
                out.pair(11 + 10 * k, _num(pts[b, k]))


def _table(out, name, handle, entries):
    out.pair(0, "TABLE")
    out.pair(2, name)
    out.pair(5, handle)
    out.pair(330, "0")
    out.pair(100, "AcDbSymbolTable")
    out.pair(70, len(entries))
    for entry in entries:
        entry(out, handle)
    out.pair(0, "ENDTAB")


def write(filename, mesh):
    """Write ``mesh`` as an ASCII DXF R2000 drawing."""
    out = _Out()
    body = _Out()
    body.handle = 0x30
    _entities(body, mesh)
    out.pair(0, "SECTION")
    out.pair(2, "HEADER")
    out.pair(9, "$ACADVER")
    out.pair(1, "AC1015")
    out.pair(9, "$HANDSEED")
    out.pair(5, format(body.handle, "X"))
    out.pair(9, "$INSUNITS")
    out.pair(70, _insunits(mesh))
    out.pair(0, "ENDSEC")
    out.pair(0, "SECTION")
    out.pair(2, "TABLES")

    def ltype(name, desc):
        def put(o, owner):
            o.pair(0, "LTYPE")
            o.pair(5, {"ByBlock": "14", "ByLayer": "15", "Continuous": "16"}[name])
            o.pair(330, owner)
            o.pair(100, "AcDbSymbolTableRecord")
            o.pair(100, "AcDbLinetypeTableRecord")
            o.pair(2, name)
            o.pair(70, 0)
            o.pair(3, desc)
            o.pair(72, 65)
            o.pair(73, 0)
            o.pair(40, "0.0")

        return put

    def layer(o, owner):
        o.pair(0, "LAYER")
        o.pair(5, "10")
        o.pair(330, owner)
        o.pair(100, "AcDbSymbolTableRecord")
        o.pair(100, "AcDbLayerTableRecord")
        o.pair(2, "0")
        o.pair(70, 0)
        o.pair(62, 7)
        o.pair(6, "Continuous")

    def style(o, owner):
        o.pair(0, "STYLE")
        o.pair(5, "11")
        o.pair(330, owner)
        o.pair(100, "AcDbSymbolTableRecord")
        o.pair(100, "AcDbTextStyleTableRecord")
        o.pair(2, "Standard")
        o.pair(70, 0)
        o.pair(40, "0.0")
        o.pair(41, "1.0")
        o.pair(50, "0.0")
        o.pair(71, 0)
        o.pair(42, "2.5")
        o.pair(3, "txt")
        o.pair(4, "")

    def appid(o, owner):
        o.pair(0, "APPID")
        o.pair(5, "12")
        o.pair(330, owner)
        o.pair(100, "AcDbSymbolTableRecord")
        o.pair(100, "AcDbRegAppTableRecord")
        o.pair(2, "ACAD")
        o.pair(70, 0)

    def dimstyle(o, owner):
        o.pair(0, "DIMSTYLE")
        o.pair(105, "13")
        o.pair(330, owner)
        o.pair(100, "AcDbSymbolTableRecord")
        o.pair(100, "AcDbDimStyleTableRecord")
        o.pair(2, "Standard")
        o.pair(70, 0)

    def block_record(name, handle, layout):
        def put(o, owner):
            o.pair(0, "BLOCK_RECORD")
            o.pair(5, handle)
            o.pair(330, owner)
            o.pair(100, "AcDbSymbolTableRecord")
            o.pair(100, "AcDbBlockTableRecord")
            o.pair(2, name)
            o.pair(340, layout)

        return put

    _table(out, "VPORT", "8", [])
    _table(
        out,
        "LTYPE",
        "5",
        [
            ltype("ByBlock", ""),
            ltype("ByLayer", ""),
            ltype("Continuous", "Solid line"),
        ],
    )
    _table(out, "LAYER", "2", [layer])
    _table(out, "STYLE", "3", [style])
    _table(out, "VIEW", "6", [])
    _table(out, "UCS", "7", [])
    _table(out, "APPID", "9", [appid])
    out.pair(0, "TABLE")
    out.pair(2, "DIMSTYLE")
    out.pair(5, "A")
    out.pair(330, "0")
    out.pair(100, "AcDbSymbolTable")
    out.pair(70, 1)
    out.pair(100, "AcDbDimStyleTable")
    dimstyle(out, "A")
    out.pair(0, "ENDTAB")
    _table(
        out,
        "BLOCK_RECORD",
        "1",
        [
            block_record("*Model_Space", "1F", "22"),
            block_record("*Paper_Space", "1B", "1E"),
        ],
    )
    out.pair(0, "ENDSEC")
    out.pair(0, "SECTION")
    out.pair(2, "BLOCKS")
    for name, rec, begin, end in (
        ("*Model_Space", "1F", "20", "21"),
        ("*Paper_Space", "1B", "1C", "1D"),
    ):
        out.pair(0, "BLOCK")
        out.pair(5, begin)
        out.pair(330, rec)
        out.pair(100, "AcDbEntity")
        out.pair(8, "0")
        out.pair(100, "AcDbBlockBegin")
        out.pair(2, name)
        out.pair(70, 0)
        out.pair(10, "0.0")
        out.pair(20, "0.0")
        out.pair(30, "0.0")
        out.pair(3, name)
        out.pair(1, "")
        out.pair(0, "ENDBLK")
        out.pair(5, end)
        out.pair(330, rec)
        out.pair(100, "AcDbEntity")
        out.pair(8, "0")
        out.pair(100, "AcDbBlockEnd")
    out.pair(0, "ENDSEC")
    out.pair(0, "SECTION")
    out.pair(2, "ENTITIES")
    out.lines += body.lines
    out.pair(0, "ENDSEC")
    out.pair(0, "SECTION")
    out.pair(2, "OBJECTS")
    for handle, owner, entries in (
        ("C", "0", (("ACAD_GROUP", "D"), ("ACAD_LAYOUT", "1A"))),
        ("D", "C", ()),
        ("1A", "C", (("Layout1", "1E"), ("Model", "22"))),
    ):
        out.pair(0, "DICTIONARY")
        out.pair(5, handle)
        out.pair(330, owner)
        out.pair(100, "AcDbDictionary")
        out.pair(281, 1)
        for name, target in entries:
            out.pair(3, name)
            out.pair(350, target)
    for handle, name, record, order in (
        ("22", "Model", "1F", 0),
        ("1E", "Layout1", "1B", 1),
    ):
        out.pair(0, "LAYOUT")
        out.pair(5, handle)
        out.pair(330, "1A")
        out.pair(100, "AcDbPlotSettings")
        out.pair(1, "")
        out.pair(4, "")
        out.pair(6, "")
        out.pair(40, "0.0")
        out.pair(41, "0.0")
        out.pair(42, "0.0")
        out.pair(43, "0.0")
        out.pair(44, "0.0")
        out.pair(45, "0.0")
        out.pair(46, "0.0")
        out.pair(47, "0.0")
        out.pair(48, "0.0")
        out.pair(49, "0.0")
        out.pair(140, "0.0")
        out.pair(141, "0.0")
        out.pair(142, "1.0")
        out.pair(143, "1.0")
        out.pair(70, 688)
        out.pair(72, 0)
        out.pair(73, 0)
        out.pair(74, 5)
        out.pair(7, "")
        out.pair(75, 16)
        out.pair(147, "1.0")
        out.pair(148, "0.0")
        out.pair(149, "0.0")
        out.pair(100, "AcDbLayout")
        out.pair(1, name)
        out.pair(70, 1)
        out.pair(71, order)
        out.pair(10, "0.0")
        out.pair(20, "0.0")
        out.pair(11, "420.0")
        out.pair(21, "297.0")
        out.pair(12, "0.0")
        out.pair(22, "0.0")
        out.pair(32, "0.0")
        out.pair(14, "0.0")
        out.pair(24, "0.0")
        out.pair(34, "0.0")
        out.pair(15, "0.0")
        out.pair(25, "0.0")
        out.pair(35, "0.0")
        out.pair(146, "0.0")
        out.pair(13, "0.0")
        out.pair(23, "0.0")
        out.pair(33, "0.0")
        out.pair(16, "1.0")
        out.pair(26, "0.0")
        out.pair(36, "0.0")
        out.pair(17, "0.0")
        out.pair(27, "1.0")
        out.pair(37, "0.0")
        out.pair(76, 0)
        out.pair(330, record)
    out.pair(0, "ENDSEC")
    out.pair(0, "EOF")
    # The provenance block leads as 999 comment pairs, which every DXF reader skips.
    text = _provenance.render_lines(_provenance.SlotTier.BLOCK, "999\n")
    text += "\n".join(out.lines) + "\n"
    if hasattr(filename, "write"):
        filename.write(text)
        return
    with open(filename, "w", encoding="ascii", newline="\n") as f:
        f.write(text)


def _pairs(text):
    lines = text.splitlines()
    if len(lines) % 2:
        lines = lines[:-1] if not lines[-1].strip() else lines
    if len(lines) % 2:
        raise ReadError("meshio++: dxf: the file ends inside a group code pair")
    for k in range(0, len(lines), 2):
        code = lines[k].strip()
        try:
            yield int(code), lines[k + 1].strip()
        except ValueError:
            raise ReadError(
                f"meshio++: dxf: line {k + 1}: {code!r} is not a group code"
            )


def _float(value, code):
    try:
        v = float(value)
    except ValueError:
        raise ReadError(f"meshio++: dxf: group {code}: {value!r} is not a number")
    if not math.isfinite(v):
        raise ReadError(f"meshio++: dxf: group {code}: not finite")
    return v


def _bspline(degree, knots, ctrl, weights, n):
    """``n + 1`` points of the (rational) B-spline by de Boor's algorithm."""
    p = degree
    ctrl = np.asarray(ctrl, dtype=np.float64)
    w = np.ones(len(ctrl)) if weights is None else np.asarray(weights, dtype=np.float64)
    hom = np.hstack([ctrl * w[:, None], w[:, None]])
    t0, t1 = knots[p], knots[len(knots) - p - 1]
    out = []
    for i in range(n + 1):
        t = t0 + (t1 - t0) * i / n
        k = p
        while k < len(knots) - p - 2 and knots[k + 1] <= t:
            k += 1
        d = [hom[j + k - p].copy() for j in range(p + 1)]
        for r in range(1, p + 1):
            for j in range(p, r - 1, -1):
                den = knots[j + 1 + k - r] - knots[j + k - p]
                alpha = 0.0 if den == 0.0 else (t - knots[j + k - p]) / den
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j]
        out.append(d[p][:3] / d[p][3])
    return out


def _arc_points(centre, radius, a0, a1, ocs, closed):
    ax, ay, n = ocs
    sweep = 360.0 if closed else (a1 - a0) % 360.0 or 360.0
    k = max(4, int(math.ceil(sweep / 360.0 * _CIRCLE_SEGMENTS)))
    count = k if closed else k + 1
    pts = []
    for i in range(count):
        a = (a0 + sweep * i / k) * math.pi / 180.0
        lx = centre[0] + radius * math.cos(a)
        ly = centre[1] + radius * math.sin(a)
        pts.append(tuple(lx * ax[r] + ly * ay[r] + centre[2] * n[r] for r in range(3)))
    return pts


def read(filename):
    """Read the LINE, ARC, CIRCLE, LWPOLYLINE and SPLINE entities of a DXF file
    as a line mesh, with ``dxf:entity`` (the entity's index) and ``dxf:type``
    (1 line, 2 arc, 3 circle, 4 polyline, 5 spline) cell data."""
    if hasattr(filename, "read"):
        text = filename.read()
    else:
        with open(filename, encoding="utf-8", errors="replace") as f:
            text = f.read()
    pairs = list(_pairs(text))
    points, index, segs, ent_of, type_of = [], {}, [], [], []

    def point(p):
        key = tuple(float(v) for v in p)
        if key not in index:
            index[key] = len(points)
            points.append(key)
        return index[key]

    def polyline(pts, closed, ent, kind):
        ids = [point(p) for p in pts]
        if closed:
            ids.append(ids[0])
        for a, b in zip(ids[:-1], ids[1:]):
            if a != b:
                segs.append((a, b))
                ent_of.append(ent)
                type_of.append(kind)

    units = 0
    k = 0
    section = None
    entity = 0
    while k < len(pairs):
        code, value = pairs[k]
        if code == 0 and value == "SECTION" and k + 1 < len(pairs):
            section = pairs[k + 1][1]
            k += 2
            continue
        if code == 0 and value == "ENDSEC":
            section = None
            k += 1
            continue
        if (
            section == "HEADER"
            and code == 9
            and value == "$INSUNITS"
            and k + 1 < len(pairs)
        ):
            units = int(_float(pairs[k + 1][1], 70))
            k += 2
            continue
        if section != "ENTITIES" or code != 0:
            k += 1
            continue
        kind = value
        if kind not in _READ:
            raise ReadError(
                f"meshio++: dxf: the entity {kind} is not read (only LINE, ARC, CIRCLE, "
                "LWPOLYLINE and SPLINE)"
            )
        k += 1
        groups = []
        while k < len(pairs) and pairs[k][0] != 0:
            groups.append(pairs[k])
            k += 1
        g = {}
        for c, v in groups:
            g.setdefault(c, []).append(v)

        def num(c, default=0.0, g=g):
            return _float(g[c][0], c) if c in g else default

        n = (num(210), num(220), num(230, 1.0))
        if not any(n):
            raise ReadError(f"meshio++: dxf: {kind}: zero extrusion direction")
        ocs = _ocs(n)
        if kind == "LINE":
            polyline(
                [(num(10), num(20), num(30)), (num(11), num(21), num(31))],
                False,
                entity,
                1,
            )
        elif kind in ("ARC", "CIRCLE"):
            r = num(40)
            if not r > 0.0:
                raise ReadError(f"meshio++: dxf: {kind}: radius must be positive")
            c = (num(10), num(20), num(30))
            closed = kind == "CIRCLE"
            a0, a1 = (0.0, 360.0) if closed else (num(50), num(51))
            polyline(
                _arc_points(c, r, a0, a1, ocs, closed),
                closed,
                entity,
                2 if kind == "ARC" else 3,
            )
        elif kind == "LWPOLYLINE":
            xs, ys = g.get(10, []), g.get(20, [])
            if len(xs) != len(ys):
                raise ReadError(
                    "meshio++: dxf: LWPOLYLINE: unpaired vertex coordinates"
                )
            bulges = [0.0] * len(xs)
            seen = 0
            for c, v in groups:
                if c == 10:
                    seen += 1
                elif c == 42 and seen > 0:
                    bulges[seen - 1] = _float(v, 42)
            elev = num(38)
            closed = int(num(70)) & 1 == 1
            verts = [(_float(x, 10), _float(y, 20)) for x, y in zip(xs, ys)]
            local = []
            count = len(verts)
            for i in range(count if closed else count - 1):
                p0, p1 = verts[i], verts[(i + 1) % count]
                local.append(p0)
                b = bulges[i]
                if b != 0.0:
                    theta = 4.0 * math.atan(b)
                    chord = math.hypot(p1[0] - p0[0], p1[1] - p0[1])
                    if chord > 0.0:
                        r = chord / (2.0 * math.sin(abs(theta) / 2.0))
                        mx, my = (p0[0] + p1[0]) / 2.0, (p0[1] + p1[1]) / 2.0
                        h = r * math.cos(theta / 2.0) * (1 if b > 0 else -1)
                        dx, dy = (p1[0] - p0[0]) / chord, (p1[1] - p0[1]) / chord
                        cx, cy = mx - h * dy, my + h * dx
                        a0 = math.atan2(p0[1] - cy, p0[0] - cx)
                        steps = max(
                            2,
                            int(
                                math.ceil(abs(theta) / (2 * math.pi) * _CIRCLE_SEGMENTS)
                            ),
                        )
                        for s in range(1, steps):
                            a = a0 + theta * s / steps
                            local.append((cx + r * math.cos(a), cy + r * math.sin(a)))
            if not closed and count:
                local.append(verts[-1])
            world = [
                tuple(
                    x * ocs[0][r] + y * ocs[1][r] + elev * ocs[2][r] for r in range(3)
                )
                for x, y in local
            ]
            polyline(world, closed, entity, 4)
        else:  # SPLINE
            degree = int(num(71, 3))
            knots = [_float(v, 40) for v in g.get(40, [])]
            ctrl = list(zip(g.get(10, []), g.get(20, []), g.get(30, [])))
            ctrl = [(_float(a, 10), _float(b, 20), _float(c, 30)) for a, b, c in ctrl]
            weights = [_float(v, 41) for v in g.get(41, [])] or None
            if (
                degree < 1
                or len(ctrl) <= degree
                or len(knots) != len(ctrl) + degree + 1
            ):
                raise ReadError(
                    "meshio++: dxf: SPLINE: inconsistent degree, knots and points"
                )
            if weights is not None and len(weights) != len(ctrl):
                raise ReadError("meshio++: dxf: SPLINE: one weight per control point")
            polyline(
                _bspline(degree, knots, ctrl, weights, 8 * len(ctrl)), False, entity, 5
            )
        entity += 1
    pts = np.array(points, dtype=np.float64).reshape(-1, 3)
    cells = (
        [CellBlock("line", np.array(segs, dtype=np.int64).reshape(-1, 2))]
        if segs
        else []
    )
    mesh = Mesh(pts, cells)
    if segs:
        mesh.cell_data["dxf:entity"] = [np.array(ent_of, dtype=np.int64)]
        mesh.cell_data["dxf:type"] = [np.array(type_of, dtype=np.int64)]
    for code, metres in _UNITS:
        if code == units:
            mesh.field_data[LENGTH_UNIT_KEY] = np.array([metres])
    return mesh
