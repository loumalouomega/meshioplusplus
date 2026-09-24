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
- ``25``/``26`` (title, summary) and every other packet are skipped; ``99``
  ends the file.

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


def read(filename):
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
    return mesh


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
    other = (
        len(mesh.point_data)
        + len(mesh.field_data)
        + len(mesh.cell_data)
        - (1 if pid_data is not None else 0)
    )
    if other:
        warn("Patran neutral holds no data arrays; point, cell and field data dropped")
        _provenance.note(
            "data-dropped",
            "a Patran neutral file holds no data arrays besides the element property",
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
