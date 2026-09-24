"""I/O for the MSC Patran 2 neutral file (``.pat``/``.out``).

The pure-Python twin of ``src/cpp/src/formats/patran.cpp``: both engines read
the same meshes and write the same bytes.

The file is a sequence of packets. Each opens with a fixed-width header card
``(I2,8I8)`` -- ``IT, ID, IV, KC, N1..N5`` -- followed by ``KC`` data cards.

- ``01`` node: ``ID`` is the node id, card 1 its coordinates ``(3E16.9)``.
- ``02`` element: ``ID`` is the element id, ``IV`` its shape (2 bar, 3 tri,
  4 quad, 5 tet, 6 pyramid, 7 wedge, 8 hex); the node count on card 1 tells
  linear from quadratic, and the property id there becomes the
  ``patran:property`` cell data. Hex20/wedge15 use the ``"patran"`` tables of
  :mod:`meshioplusplus._node_order`.
- ``21`` named component: type 5 (node) entries become a point region, element
  entries (6 bar ... 12 hex) a cell region of the same name, tagged with the
  component number.
- Loads and boundary conditions, per load or constraint set ``<set>``:
  ``06`` distributed loads as the ``patran:distributed_load`` (flags) and
  ``patran:distributed_load_values`` field data tables; ``07`` node forces and
  ``08`` node displacements as ``patran:force:<set>`` and
  ``patran:displacement:<set>`` point data (six components, NaN where not
  given; a nonzero coordinate frame in ``..._frame:<set>``); ``10`` node and
  ``11`` element temperatures as ``patran:temperature:<set>`` point data and
  ``patran:element_temperature:<set>`` cell data.
- ``25``/``26`` (title, summary) and every other packet are skipped; ``99``
  ends the file.

Patran 2.5 result files (``.nod``/``.dis`` nodal, ``.els`` element; text or
binary) are read onto the mesh with ``read(filename, results={name: path})``.

Elements no component names are grouped by property into ``property_<pid>``
cell regions (tag = pid).
"""

import re

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import open_file
from .._mesh import Mesh, topological_dimension
from .._node_order import node_order
from .._regions import Region

__all__ = ["read", "write"]

_MAX_ID = 99999999
_MAX_NAME = 12
_EXPONENT_WITHOUT_LETTER = re.compile(r"^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d+)$")

# (shape code, node count) -> (meshio++ type, packet 21 entity code)
_SHAPES = {
    (2, 2): ("line", 6),
    (2, 3): ("line3", 6),
    (3, 3): ("triangle", 7),
    (3, 6): ("triangle6", 7),
    (3, 7): ("triangle7", 7),
    (4, 4): ("quad", 8),
    (4, 8): ("quad8", 8),
    (4, 9): ("quad9", 8),
    (5, 4): ("tetra", 9),
    (5, 10): ("tetra10", 9),
    (6, 5): ("pyramid", 10),
    (6, 13): ("pyramid13", 10),
    (7, 6): ("wedge", 11),
    (7, 15): ("wedge15", 11),
    (8, 8): ("hexahedron", 12),
    (8, 20): ("hexahedron20", 12),
}
# meshio++ type -> (shape code, node count, packet 21 entity code)
_SHAPE_OF = {t: (s, k, code) for (s, k), (t, code) in _SHAPES.items()}


def _fail(what, line):
    raise ReadError(f"Patran neutral: {what} (line {line})")


def _int(text, line):
    text = text.strip()
    if not text:
        return 0
    try:
        return int(text)
    except ValueError:
        raise ReadError(
            f"Patran neutral: invalid integer field {text!r} (line {line})"
        ) from None


def _real(text, line):
    text = text.strip()
    if not text:
        return 0.0
    s = text.replace("D", "E").replace("d", "e")
    m = _EXPONENT_WITHOUT_LETTER.match(s)  # Fortran's 1.5-3
    if m:
        s = m.group(1) + "e" + m.group(2)
    try:
        return float(s)
    except ValueError:
        raise ReadError(
            f"Patran neutral: invalid real field {text!r} (line {line})"
        ) from None


def _fields(line, widths):
    out = []
    pos = 0
    for w in widths:
        if pos >= len(line):
            break
        out.append(line[pos : pos + w])
        pos += w
    return out


_HEADER = (2,) + (8,) * 8
_LOAD_FLAGS = (1,) * 17 + (2,)  # packet 06 data card 1: (3I1,6I1,8I1,I2)
_FRAME_FLAGS = (8,) + (1,) * 6  # packets 07/08 data card 1: (I8,6I1)
_REALS = (16,) * 5
_XYZ = (16, 16, 16)
_ELEM = (8, 8, 8, 8, 16, 16, 16)
_INTS = (8,) * 10


def _header(line, line_no):
    f = [_int(t, line_no) for t in _fields(line, _HEADER)]
    f += [0] * (9 - len(f))
    if f[3] < 0:
        _fail(f"negative card count {f[3]}", line_no)
    return f


def _int_cards(lines, first, count):
    out = []
    for c in range(count):
        f = [_int(t, first + c + 1) for t in _fields(lines[first + c], _INTS)]
        out += f + [0] * (10 - len(f))
    return out


def _reals(lines, first, count, n, line_no):
    """The first ``n`` ``E16.9`` fields of the cards ``lines[first:first+count]``."""
    out = []
    for c in range(count):
        out += [_real(t, first + c + 1) for t in _fields(lines[first + c], _REALS)]
    if len(out) < n:
        _fail(f"expected {n} values, found {len(out)}", line_no)
    return out[:n]


def read(filename, results=None):
    with open_file(filename, "rb") as f:
        raw = f.read()
    text = raw.decode("latin-1") if isinstance(raw, bytes) else raw
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    lines = [ln[:-1] if ln.endswith("\r") else ln for ln in lines]

    node_ids = []
    coords = []
    elements = []  # (id, type, pid, node ids, line)
    components = []  # (number, name, [(code, id)])
    loads = []  # (packet, id, set, fields, values, line)
    warned_shapes = set()
    saw_end = False

    i = 0
    n = len(lines)
    while i < n:
        if not lines[i].strip(" \t"):
            i += 1
            continue
        head_line = i + 1
        it, ident, iv, kc, *_ = _header(lines[i], head_line)
        i += 1
        if it == 99:
            saw_end = True
            break
        if i + kc > n:
            _fail(
                f"packet {it} announces {kc} data cards but the file ends first",
                head_line,
            )
        if it == 1:
            if kc < 1:
                _fail(f"node {ident} has no coordinate card", head_line)
            f = _fields(lines[i], _XYZ)
            node_ids.append(ident)
            coords.append([_real(f[d], i + 1) if d < len(f) else 0.0 for d in range(3)])
        elif it == 2:
            if kc < 1:
                _fail(f"element {ident} has no data card", head_line)
            f = _fields(lines[i], _ELEM)
            nodes = _int(f[0], i + 1) if f else 0
            pid = _int(f[2], i + 1) if len(f) > 2 else 0
            if nodes < 0:
                _fail(f"element {ident} has a negative node count", i + 1)
            node_cards = (nodes + 9) // 10
            if 1 + node_cards > kc:
                _fail(
                    f"element {ident} lists {nodes} nodes but has {kc} data cards",
                    head_line,
                )
            shape = _SHAPES.get((iv, nodes))
            if shape is None:
                if (iv, nodes) not in warned_shapes:
                    warned_shapes.add((iv, nodes))
                    warn(
                        f"Patran neutral: skipping elements of shape {iv} with {nodes} "
                        f"nodes (no meshio++ equivalent; first one is element {ident})"
                    )
            else:
                ids = _int_cards(lines, i + 1, node_cards)[:nodes]
                elements.append((ident, shape[0], pid, ids, head_line))
        elif it == 21:
            if kc < 1:
                _fail(f"component {ident} has no name card", head_line)
            name = lines[i].strip(" \t")
            values = _int_cards(lines, i + 1, kc - 1)
            count = min(max(iv, 0), len(values))
            pairs = [(values[k], values[k + 1]) for k in range(0, count - 1, 2)]
            components.append((ident, name, pairs))
        elif it in (7, 8):
            if kc < 1:
                _fail(f"packet {it} of node {ident} has no data card", head_line)
            f = [_int(t, i + 1) for t in _fields(lines[i], _FRAME_FLAGS)]
            f += [0] * (7 - len(f))
            given = sum(1 for v in f[1:] if v)
            values = _reals(lines, i + 1, kc - 1, given, head_line)
            loads.append((it, ident, iv, f, values, head_line))
        elif it == 6:
            if kc < 1:
                _fail(
                    f"distributed load on element {ident} has no data card", head_line
                )
            f = [_int(t, i + 1) for t in _fields(lines[i], _LOAD_FLAGS)]
            f += [0] * (18 - len(f))
            nc = sum(1 for v in f[3:9] if v)
            nn = sum(1 for v in f[9:17] if v)
            npv = nc * ((1 if f[1] else 0) + nn * (1 if f[2] else 0))
            values = _reals(lines, i + 1, kc - 1, npv, head_line)
            loads.append((it, ident, iv, f, values, head_line))
        elif it in (10, 11):
            if kc < 1:
                _fail(f"temperature packet {it} of {ident} has no data card", head_line)
            header = _header(lines[head_line - 1], head_line)
            value = _real(lines[i][:16], i + 1)
            loads.append((it, ident, iv, [header[4]], [value], head_line))
        i += kc
    if not saw_end:
        warn(
            f"Patran neutral: '{filename}' has no end packet (99); reading it to the end"
        )

    node_index = {}
    for p, ident in enumerate(node_ids):
        if ident in node_index:
            raise ReadError(f"Patran neutral: node {ident} is defined twice")
        node_index[ident] = p
    points = np.array(coords, dtype=np.float64).reshape(len(node_ids), 3)

    order_of_types = []
    by_type = {}
    for e, el in enumerate(elements):
        if el[1] not in by_type:
            by_type[el[1]] = []
            order_of_types.append(el[1])
        by_type[el[1]].append(e)

    cells = []
    pid_blocks = []
    element_index = {}
    cell_pid = []
    cell_dim = []
    for cell_type in order_of_types:
        members = by_type[cell_type]
        k = _SHAPE_OF[cell_type][1]
        order = node_order("patran", cell_type)
        to_meshio = order.to_meshio if order else range(k)
        dim = topological_dimension[cell_type]
        conn = np.empty((len(members), k), dtype=np.int64)
        pids = np.empty(len(members), dtype=np.int64)
        for r, e in enumerate(members):
            ident, _, pid, ids, line = elements[e]
            for j, src in enumerate(to_meshio):
                idx = node_index.get(ids[src])
                if idx is None:
                    _fail(f"element {ident} names undefined node {ids[src]}", line)
                conn[r, j] = idx
            if ident in element_index:
                _fail(f"element {ident} is defined twice", line)
            element_index[ident] = len(cell_pid)
            pids[r] = pid
            cell_pid.append(pid)
            cell_dim.append(dim)
        cells.append((cell_type, conn))
        pid_blocks.append(pids)

    mesh = Mesh(points, cells)
    if pid_blocks:
        mesh.cell_data["patran:property"] = pid_blocks

    regions = []
    named = [False] * len(cell_pid)
    warned_codes = set()
    for number, name, pairs in components:
        pts, cls = [], []
        seen_pts, seen_cls = set(), set()
        dim = -1
        missing = 0
        for code, ident in pairs:
            if code == 5:
                idx = node_index.get(ident)
                if idx is None:
                    missing += 1
                elif idx not in seen_pts:
                    seen_pts.add(idx)
                    pts.append(idx)
            elif 6 <= code <= 12:
                idx = element_index.get(ident)
                if idx is None:
                    missing += 1
                elif idx not in seen_cls:
                    seen_cls.add(idx)
                    cls.append(idx)
                    named[idx] = True
                    dim = max(dim, cell_dim[idx])
            elif code not in warned_codes:
                warned_codes.add(code)
                warn(
                    f"Patran neutral: component '{name}' lists entities of type {code}; "
                    "skipped"
                )
        if missing:
            warn(
                f"Patran neutral: component '{name}' names {missing} undefined node(s) "
                "or element(s)"
            )
        if cls or not pts:
            regions.append(
                Region(name, "cell", np.array(cls, dtype=np.int64), dim, number)
            )
        if pts:
            regions.append(
                Region(name, "point", np.array(pts, dtype=np.int64), -1, number)
            )

    by_pid = {}
    pid_dim = {}
    for g, pid in enumerate(cell_pid):
        if named[g]:
            continue
        by_pid.setdefault(pid, []).append(g)
        pid_dim[pid] = max(pid_dim.get(pid, -1), cell_dim[g])
    for pid in sorted(by_pid):
        regions.append(
            Region(
                f"property_{pid}",
                "cell",
                np.array(by_pid[pid], dtype=np.int64),
                pid_dim[pid],
                pid,
            )
        )
    mesh.regions = regions
    _attach_loads(mesh, loads, node_index, element_index, cells)
    for name, path in _result_items(results):
        _attach_result(mesh, name, path, node_index, element_index, cells)
    return mesh


def _attach_loads(mesh, loads, node_index, element_index, cells):
    """Packets 06, 07, 08, 10 and 11 as point, cell and field data."""
    npts = len(mesh.points)
    ncells = len(element_index)
    offsets = np.cumsum([0] + [len(c) for _, c in cells])
    point = {}  # key -> array
    frames = {}
    element_temps = {}
    dist_flags, dist_values = [], []
    missing = 0
    for packet, ident, iv, f, values, _ in loads:
        if packet == 6:
            idx = element_index.get(ident)
            if idx is None:
                missing += 1
                continue
            dist_flags.append([idx, iv] + list(f))
            dist_values.append(values + [np.nan] * (54 - len(values)))
            continue
        if packet == 11:
            idx = element_index.get(ident)
            if idx is None:
                missing += 1
                continue
            arr = element_temps.setdefault(iv, np.full(ncells, np.nan))
            arr[idx] = values[0] if f[0] else np.nan
            continue
        idx = node_index.get(ident)
        if idx is None:
            missing += 1
            continue
        if packet == 10:
            arr = point.setdefault(f"patran:temperature:{iv}", np.full(npts, np.nan))
            arr[idx] = values[0] if f[0] else np.nan
            continue
        kind = "force" if packet == 7 else "displacement"
        arr = point.setdefault(f"patran:{kind}:{iv}", np.full((npts, 6), np.nan))
        k = 0
        for c in range(6):
            if f[1 + c]:
                arr[idx, c] = values[k]
                k += 1
        if f[0]:
            frame = frames.setdefault(
                f"patran:{kind}_frame:{iv}", np.zeros(npts, dtype=np.int64)
            )
            frame[idx] = f[0]
    if missing:
        warn(
            f"Patran neutral: {missing} load or boundary condition record(s) name an "
            "undefined node or element; skipped"
        )
    for key in sorted(point):
        mesh.point_data[key] = point[key]
    for key in sorted(frames):
        mesh.point_data[key] = frames[key]
    for s_id in sorted(element_temps):
        a = element_temps[s_id]
        mesh.cell_data[f"patran:element_temperature:{s_id}"] = [
            a[offsets[b] : offsets[b + 1]] for b in range(len(cells))
        ]
    if dist_flags:
        mesh.field_data["patran:distributed_load"] = np.array(
            dist_flags, dtype=np.int64
        )
        mesh.field_data["patran:distributed_load_values"] = np.array(
            dist_values, dtype=np.float64
        )


def _result_items(results):
    if not results:
        return []
    if isinstance(results, dict):
        return [(str(k), str(v)) for k, v in results.items()]
    import pathlib

    return [(pathlib.Path(p).stem, str(p)) for p in results]


def _parse_result(path):
    """A Patran 2.5 result file: ``(kind, width, [(id, values)])``, ``kind``
    ``"nodal"`` (``.nod``/``.dis``: NODID and NWIDTH values) or ``"element"``
    (``.els``: ID, shape and NWIDTH values); text or Fortran unformatted
    binary (4-byte words, 4- or 8-byte reals, either byte order)."""
    with open_file(path, "rb") as f:
        raw = f.read()
    for order in ("<", ">"):
        if len(raw) >= 4:
            head = int(np.frombuffer(raw[:4], dtype=order + "i4")[0])
            if head in (340, 324):
                return _parse_binary_result(raw, order, head == 340, path)
    text = raw.decode("latin-1").replace("\r", "")
    lines = text.split("\n")
    if len(lines) < 4:
        raise ReadError(f"Patran results: {path} is too short")
    head = lines[1].split()
    nodal = len(head) >= 3
    width = _int(head[-1] if nodal else lines[1][:5], 2)
    if width < 1:
        raise ReadError(f"Patran results: {path} has {width} columns")
    rows = []
    i = 4
    n = len(lines)
    while i < n:
        line = lines[i]
        if not line.strip():
            i += 1
            continue
        if nodal:
            ident = _int(line[:8], i + 1)
            vals = [_real(line[8 + 13 * k : 21 + 13 * k], i + 1) for k in range(5)]
            vals = [
                v for k, v in enumerate(vals) if line[8 + 13 * k : 21 + 13 * k].strip()
            ]
            per_line = 5
        else:
            f = _fields(line, (8, 8))
            ident = _int(f[0], i + 1)
            if ident == 0:
                break
            vals = []
            per_line = 6
        i += 1
        while len(vals) < width:
            if i >= n:
                raise ReadError(f"Patran results: {path} ends inside record {ident}")
            cont = lines[i]
            fields = [cont[13 * k : 13 * k + 13] for k in range(per_line)]
            vals += [_real(t, i + 1) for t in fields if t.strip()]
            i += 1
        rows.append((ident, vals[:width]))
    return ("nodal" if nodal else "element"), width, rows


def _parse_binary_result(raw, order, nodal, path):
    i4 = np.dtype(order + "i4")
    pos = 0
    records = []
    while pos + 4 <= len(raw):
        size = int(np.frombuffer(raw, dtype=i4, count=1, offset=pos)[0])
        if size < 0 or pos + 8 + size > len(raw):
            raise ReadError(f"Patran results: {path} has a truncated record")
        records.append(raw[pos + 4 : pos + 4 + size])
        pos += 8 + size
    if len(records) < 3:
        raise ReadError(f"Patran results: {path} is too short")
    first = records[0]
    width = int(np.frombuffer(first, dtype=i4, count=1, offset=len(first) - 4)[0])
    lead = 1 if nodal else 2
    rows = []
    for rec in records[3:]:
        ident = int(np.frombuffer(rec, dtype=i4, count=1)[0])
        if ident == 0:
            break
        real = len(rec) - 4 * lead
        if real == 4 * width:
            dt = np.dtype(order + "f4")
        elif real == 8 * width:
            dt = np.dtype(order + "f8")
        else:
            raise ReadError(
                f"Patran results: a record of {path} holds {len(rec)} bytes for "
                f"{width} columns"
            )
        vals = np.frombuffer(rec, dtype=dt, count=width, offset=4 * lead)
        rows.append((ident, [float(v) for v in vals]))
    return ("nodal" if nodal else "element"), width, rows


def _attach_result(mesh, name, path, node_index, element_index, cells):
    kind, width, rows = _parse_result(path)
    index = node_index if kind == "nodal" else element_index
    size = len(mesh.points) if kind == "nodal" else len(element_index)
    data = np.full((size, width), np.nan)
    missing = 0
    for ident, vals in rows:
        idx = index.get(ident)
        if idx is None:
            missing += 1
            continue
        data[idx] = vals
    if missing:
        warn(
            f"Patran results: {missing} record(s) of {path} name an undefined "
            f"{'node' if kind == 'nodal' else 'element'}; skipped"
        )
    if width == 1:
        data = data[:, 0].copy()
    if kind == "nodal":
        mesh.point_data[name] = data
        return
    offsets = np.cumsum([0] + [len(c) for _, c in cells])
    mesh.cell_data[name] = [
        data[offsets[b] : offsets[b + 1]].copy() for b in range(len(cells))
    ]


def _header_card(it, ident, iv, kc, n1=0, n2=0, n3=0, n4=0, n5=0):
    return f"{it:2d}{ident:8d}{iv:8d}{kc:8d}{n1:8d}{n2:8d}{n3:8d}{n4:8d}{n5:8d}\n"


def _int_lines(values):
    out = []
    for k in range(0, len(values), 10):
        out.append("".join(f"{v:8d}" for v in values[k : k + 10]) + "\n")
    return "".join(out)


def _component_name(name, taken):
    clean = (name or "COMPONENT").replace("\n", " ").replace("\r", " ")
    clean = clean[:_MAX_NAME]
    out = clean
    k = 1
    while out in taken:
        suffix = f"_{k}"
        out = clean[: min(len(clean), _MAX_NAME - len(suffix))] + suffix
        k += 1
    taken.add(out)
    if out != name:
        warn(f"Patran neutral: component '{name}' is written as '{out}'")
    return out


_LOAD_KEY = re.compile(
    r"^patran:(force|displacement|temperature|element_temperature)"
    r"(_frame)?:(-?\d+)$"
)


def _load_arrays(mesh, npts):
    """The load and boundary-condition arrays the writer turns into packets
    06, 07, 08, 10 and 11, keyed by (kind, set), and how many there are."""
    out = {"point": {}, "frame": {}, "cell": {}, "count": 0, "dist": None}
    for key, a in mesh.point_data.items():
        m = _LOAD_KEY.match(key)
        a = np.asarray(a)
        if not m or m.group(1) == "element_temperature":
            continue
        kind, frame, set_id = m.group(1), m.group(2), int(m.group(3))
        want = (npts,) if kind == "temperature" or frame else (npts, 6)
        if a.shape != want:
            continue
        out["frame" if frame else "point"][(kind, set_id)] = a
        out["count"] += 1
    for key, blocks in mesh.cell_data.items():
        m = _LOAD_KEY.match(key)
        if m and m.group(1) == "element_temperature" and not m.group(2):
            out["cell"][int(m.group(3))] = [np.asarray(b).ravel() for b in blocks]
            out["count"] += 1
    flags = mesh.field_data.get("patran:distributed_load")
    values = mesh.field_data.get("patran:distributed_load_values")
    if flags is not None and values is not None:
        flags, values = np.asarray(flags), np.asarray(values)
        if flags.ndim == 2 and flags.shape[1] == 20 and len(values) == len(flags):
            out["dist"] = (flags, values)
            out["count"] += 2
    return out


def _real_lines(values):
    out = []
    for k in range(0, len(values), 5):
        out.append("".join("%16.9E" % v for v in values[k : k + 5]) + "\n")
    return "".join(out)


def _load_packets(mesh, loads, cell_label):
    """Packets 06 (distributed loads), 07 (node forces), 08 (node
    displacements), 10 (node temperatures) and 11 (element temperatures)."""
    out = []
    if loads["dist"] is not None:
        flags, values = loads["dist"]
        for row, vals in zip(flags, values):
            cell = int(row[0])
            if not 0 <= cell < len(cell_label) or not cell_label[cell]:
                continue
            f = [int(v) for v in row[2:]]
            nc = sum(1 for v in f[3:9] if v)
            nn = sum(1 for v in f[9:17] if v)
            npv = nc * ((1 if f[1] else 0) + nn * (1 if f[2] else 0))
            vals = [float(v) for v in vals[:npv]]
            out.append(
                _header_card(6, cell_label[cell], int(row[1]), 1 + (npv + 4) // 5)
            )
            out.append("".join(str(v) for v in f[:17]) + f"{f[17]:2d}\n")
            out.append(_real_lines(vals))
    for packet, kind in ((7, "force"), (8, "displacement")):
        for (k, set_id), a in sorted(loads["point"].items()):
            if k != kind:
                continue
            frame = loads["frame"].get((kind, set_id))
            for p in range(len(a)):
                given = [not np.isnan(v) for v in a[p]]
                if not any(given):
                    continue
                vals = [float(v) for v, g in zip(a[p], given) if g]
                out.append(
                    _header_card(packet, p + 1, set_id, 1 + (len(vals) + 4) // 5)
                )
                cid = int(frame[p]) if frame is not None else 0
                out.append(
                    f"{cid:8d}" + "".join("1" if g else "0" for g in given) + "\n"
                )
                out.append(_real_lines(vals))
    for (k, set_id), a in sorted(loads["point"].items()):
        if k != "temperature":
            continue
        for p, v in enumerate(a):
            if not np.isnan(v):
                out.append(_header_card(10, p + 1, set_id, 1, 1))
                out.append("%16.9E\n" % float(v))
    for set_id, blocks in sorted(loads["cell"].items()):
        g = 0
        for a in blocks:
            for v in a:
                if not np.isnan(v) and g < len(cell_label) and cell_label[g]:
                    out.append(_header_card(11, cell_label[g], set_id, 1, 1))
                    out.append("%16.9E\n" % float(v))
                g += 1
    return out


def write(filename, mesh):
    points = np.asarray(mesh.points)
    pdim = points.shape[1] if points.ndim == 2 else 0
    if pdim > 3:
        raise WriteError(
            f"Patran neutral writer: points of dimension {pdim} (at most 3)"
        )
    npts = len(points)

    shapes = []
    cell_label = []
    written = 0
    dropped_types = set()
    for block in mesh.cells:
        shape = (
            _SHAPE_OF.get(block.type) if isinstance(block.data, np.ndarray) else None
        )
        shapes.append(shape)
        if shape is None:
            dropped_types.add(block.type)
        for _ in range(len(block.data)):
            if shape is None:
                cell_label.append(0)
            else:
                written += 1
                cell_label.append(written)
    if npts > _MAX_ID or written > _MAX_ID:
        raise WriteError(
            "Patran neutral writer: more than 99,999,999 nodes or elements do not fit "
            "the I8 id fields"
        )

    for t in sorted(dropped_types):
        warn(f"Patran neutral has no '{t}' element; those cells are dropped")
        _provenance.note("cells-dropped", f"Patran neutral has no '{t}' element")
    pid_data = mesh.cell_data.get("patran:property")
    regions = sorted(getattr(mesh, "regions", []) or [], key=lambda r: r.key)
    side_regions = sum(1 for r in regions if r.kind == "side")
    if side_regions:
        warn(
            f"Patran neutral has no facet components; {side_regions} side region(s) "
            "dropped"
        )
        _provenance.note(
            "regions-dropped",
            f"{side_regions} side region(s) have no Patran neutral component",
        )
    loads = _load_arrays(mesh, npts)
    other = (
        len(mesh.point_data)
        + len(mesh.field_data)
        + len(mesh.cell_data)
        - (1 if pid_data is not None else 0)
        - loads["count"]
    )
    if other:
        warn(
            "Patran neutral holds no data arrays but the element property and the "
            "loads; the other point, cell and field data are dropped"
        )
        _provenance.note(
            "data-dropped",
            "a Patran neutral file holds no data arrays besides the element property "
            "and its loads and boundary conditions",
        )

    out = []
    title = _provenance.lines(_provenance.SlotTier.BOUNDED)[0][:80]
    out.append(_header_card(25, 0, 0, 1))
    out.append(title + "\n")

    if pid_data is not None:
        pids = set()
        for b, block in enumerate(mesh.cells):
            if shapes[b] is not None:
                pids.update(int(v) for v in np.asarray(pid_data[b]).ravel())
    else:
        pids = {1} if written else set()
    out.append(_header_card(26, 0, 0, 1, npts, written, 0, len(pids), 0))
    out.append(" " * 40 + "3.0\n")

    pts = points.astype(np.float64, copy=False)
    for p in range(npts):
        out.append(_header_card(1, p + 1, 0, 2))
        xyz = [float(pts[p, d]) if d < pdim else 0.0 for d in range(3)]
        out.append("".join("%16.9E" % v for v in xyz))
        out.append("\n1G       6       0       0  000000\n")

    zeros = "%16.9E" % 0.0
    for b, block in enumerate(mesh.cells):
        shape = shapes[b]
        if shape is None:
            continue
        code, k, _ = shape
        order = node_order("patran", block.type)
        data = np.asarray(block.data, dtype=np.int64)
        if order is not None:
            data = data[:, list(order.from_meshio)]
        pid = None if pid_data is None else np.asarray(pid_data[b]).ravel()
        kc = 1 + (k + 9) // 10
        g0 = sum(len(mesh.cells[q].data) for q in range(b))
        for r in range(len(data)):
            out.append(_header_card(2, cell_label[g0 + r], code, kc))
            p = 1 if pid is None else int(pid[r])
            out.append(f"{k:8d}{0:8d}{p:8d}{0:8d}" + zeros * 3 + "\n")
            out.append(_int_lines([int(v) + 1 for v in data[r]]))

    out += _load_packets(mesh, loads, cell_label)

    cell_code = []
    for b, block in enumerate(mesh.cells):
        c = shapes[b][2] if shapes[b] is not None else 0
        cell_code += [c] * len(block.data)
    names = []
    comps = {}
    tags = {}
    for region in regions:
        if region.kind == "side":
            continue
        if region.name not in comps:
            comps[region.name] = []
            names.append(region.name)
            tags[region.name] = region.tag
        entries = comps[region.name]
        for e in np.asarray(region.entries, dtype=np.int64).ravel():
            e = int(e)
            if region.kind == "point":
                entries.append((5, e + 1))
            elif e < len(cell_label) and cell_label[e]:
                entries.append((cell_code[e], cell_label[e]))
    used = {tags[nm] for nm in names if tags[nm] > 0}
    assigned = set()
    taken = set()
    nxt = 1
    for name in names:
        number = tags[name]
        if number <= 0 or number in assigned:
            while nxt in used or nxt in assigned:
                nxt += 1
            number = nxt
        assigned.add(number)
        flat = [v for pair in comps[name] for v in pair]
        iv = len(flat)
        out.append(_header_card(21, number, iv, 1 + (iv + 9) // 10))
        out.append(_component_name(name, taken) + "\n")
        out.append(_int_lines(flat))
    out.append(_header_card(99, 0, 0, 1))

    with open_file(filename, "w", newline="\n") as f:
        f.write("".join(out))
