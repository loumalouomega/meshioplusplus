"""I/O for the Femap neutral file (``.neu``).

The pure-Python twin of ``src/cpp/src/formats/femap.cpp``: both engines read the
same meshes and write the same bytes.

A neutral file is a sequence of blocks, each opened by a ``-1`` line and the
block id and closed by the next ``-1`` line; records are comma-separated. Their
layout changes with the Femap version, so every record is read by position and
length:

- ``403`` nodes (x, y, z are fields 11-13 in every version);
- ``404`` elements: seven lines, then (from 4.5) one node list per non-zero list
  flag. The topology code picks the cell type, the nodes sit in a 20-slot
  "degenerate brick" layout, and the property and element type become the
  ``femap:property`` and ``femap:type`` cell data;
- ``402`` properties name the ``property_<id>`` cell regions by their titles;
- ``408`` groups become point and cell regions;
- ``450`` output sets are steps (``time_step``); the ``451``/``1051`` output
  vectors of the selected set become point or cell data, NaN where missing.
"""

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import is_buffer, open_file
from .._mesh import Mesh, topological_dimension
from .._regions import Region

__all__ = ["read", "write", "time_values", "SeriesWriter"]

# topology code -> (meshio++ type, slot of each meshio++ node, default element type)
_TOPOLOGIES = {
    0: ("line", [0, 1], 1),
    1: ("line3", [0, 1, 2], 1),
    2: ("triangle", [0, 1, 2], 17),
    3: ("triangle6", [0, 1, 2, 4, 5, 6], 18),
    4: ("quad", [0, 1, 2, 3], 17),
    5: ("quad8", [0, 1, 2, 3, 4, 5, 6, 7], 18),
    6: ("tetra", [0, 1, 2, 4], 25),
    7: ("wedge", [0, 1, 2, 4, 5, 6], 25),
    8: ("hexahedron", [0, 1, 2, 3, 4, 5, 6, 7], 25),
    9: ("vertex", [0], 27),
    10: ("tetra10", [0, 1, 2, 4, 8, 9, 10, 12, 13, 14], 26),
    11: ("wedge15", [0, 1, 2, 4, 5, 6, 8, 9, 10, 16, 17, 18, 12, 13, 14], 26),
    12: ("hexahedron20", list(range(12)) + [16, 17, 18, 19, 12, 13, 14, 15], 26),
    14: ("pyramid", [0, 1, 2, 3, 4], 25),
}
_TOPOLOGY_OF = {t: (code, slots, et) for code, (t, slots, et) in _TOPOLOGIES.items()}
_SKIPPED_NAMES = {
    13: "rigid",
    15: "multi-list",
    16: "contact",
    17: "weld",
    18: "rigid",
    19: "pyramid13",
}


def _fields(line):
    out = [f.strip(" \t") for f in line.split(",")]
    while out and out[-1] == "":
        out.pop()
    return out


def _parse_int(text):
    body = text[1:] if text[:1] in "+-" else text
    if not body or not body.isascii() or not body.isdigit():
        return None
    return int(text)


def _parse_real(text):
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


class _Cursor:
    def __init__(self, bid, first_line, lines):
        self.bid = bid
        self.first = first_line
        self.lines = lines
        self.pos = 0

    def at_end(self):
        return self.pos >= len(self.lines)

    def remaining(self):
        return len(self.lines) - self.pos

    def line(self):
        return self.first + self.pos

    def peek(self, ahead=0):
        return self.lines[self.pos + ahead][1]

    def fail(self, why):
        raise ReadError(f"Femap neutral: {why} (line {self.line()})")

    def next(self, what):
        if self.at_end():
            self.fail(f"block {self.bid} ends inside {what}")
        text = self.lines[self.pos][1]
        self.pos += 1
        return text

    def fields(self, what):
        return _fields(self.next(what))

    def int(self, f, k, what):
        v = _parse_int(f[k]) if k < len(f) else None
        if v is None:
            self.fail(f"bad {what}" + (f" '{f[k]}'" if k < len(f) else ""))
        return v

    def real(self, f, k, what):
        v = _parse_real(f[k]) if k < len(f) else None
        if v is None:
            self.fail(f"bad {what}" + (f" '{f[k]}'" if k < len(f) else ""))
        return v


def _title(line):
    t = line.strip(" \t")
    return "" if t == "<NULL>" else t


def _blocks(text):
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    lines = [ln[:-1] if ln.endswith("\r") else ln for ln in lines]
    blocks = []
    i = 0
    n = len(lines)
    while i < n:
        if lines[i].strip(" \t") != "-1":
            i += 1
            continue
        if i + 1 >= n:
            break
        found = lines[i + 1].strip(" \t")
        bid = _parse_int(found)
        if bid is None:
            raise ReadError(
                f"Femap neutral: expected a block id after '-1', found '{found}' "
                f"(line {i + 2})"
            )
        body = []
        j = i + 2
        while j < n and lines[j].strip(" \t") != "-1":
            if lines[j].strip(" \t"):
                body.append((j + 1, lines[j]))
            j += 1
        if j >= n:
            warn(f"Femap neutral: block {bid} (line {i + 1}) has no closing -1")
        blocks.append((bid, i + 3, body))
        i = j + 1
    return blocks


class _File:
    def __init__(self):
        self.version = 0.0
        self.node_ids = []
        self.coords = []
        self.elements = []  # (id, prop, type, topology code, slots, line)
        self.properties = {}
        self.groups = []  # (id, title, nodes, elements)
        self.sets = []  # (id, title, value)
        self.vectors = []  # (set, id, title, entity, [(id, value)])


class _BlockCursor(_Cursor):
    def __init__(self, bid, first, body):
        super().__init__(bid, first, body)

    def line(self):
        if self.pos < len(self.lines):
            return self.lines[self.pos][0]
        return self.lines[-1][0] + 1 if self.lines else self.first


def _read_nodes(c, f):
    while not c.at_end():
        fl = c.fields("a node")
        if len(fl) < 14:
            c.fail(f"a node record with {len(fl)} fields (x, y, z are 11-13)")
        f.node_ids.append(c.int(fl, 0, "node id"))
        f.coords.append([c.real(fl, k, "coordinate") for k in (11, 12, 13)])


def _skip_list(c):
    while True:
        fl = c.fields("an element node list")
        if fl and fl[0] == "-1":
            return


def _read_elements(c, f, warned):
    while not c.at_end():
        line = c.line()
        head = c.fields("an element")
        eid = c.int(head, 0, "element id")
        prop = c.int(head, 2, "element property")
        etype = c.int(head, 3, "element type")
        topology = c.int(head, 4, "element topology")
        slots = [0] * 20
        for part in range(2):
            fl = c.fields("an element's nodes")
            for k in range(10):
                if k < len(fl):
                    v = _parse_int(fl[k])
                    if v is None:
                        c.fail(f"bad node id '{fl[k]}'")
                    slots[part * 10 + k] = v
        for _ in range(3):
            c.next("an element record")
        last = c.fields("an element record")
        for k in range(12, min(16, len(last))):
            flag = _parse_int(last[k])
            if flag and f.version >= 4.5 - 1e-9:
                _skip_list(c)
        if topology not in _TOPOLOGIES:
            if topology not in warned:
                warned.add(topology)
                warn(
                    f"Femap neutral: skipping {_SKIPPED_NAMES.get(topology, 'unknown')} "
                    f"elements (topology {topology}; first one is {eid})"
                )
            continue
        f.elements.append((eid, prop, etype, topology, slots, line))


def _read_properties(c, f):
    while not c.at_end():
        head = c.fields("a property")
        pid = c.int(head, 0, "property id")
        f.properties[pid] = _title(c.next("a property title"))
        c.next("property flags")
        values = 0
        for per_line, what in ((8, "laminate count"), (5, "property value count")):
            count = c.int(c.fields(what), 0, what)
            if count < 0:
                c.fail(f"negative {what}")
            for _ in range((count + per_line - 1) // per_line):
                c.next(what)
            values = count
        # Femap 2401 follows the values with as many integers, five to a line
        # (function references): a repeat of the value count before a line of
        # several integers.
        if c.remaining() >= 2 and values > 0:
            fl = _fields(c.peek())
            nxt = _fields(c.peek(1))
            if (
                len(fl) == 1
                and _parse_int(fl[0]) == values
                and len(nxt) > 1
                and all(_parse_int(v) is not None for v in nxt)
            ):
                c.next("a function count")
                for _ in range((values + 4) // 5):
                    c.next("a function reference")
        while not c.at_end():
            fl = _fields(c.peek())
            count = _parse_int(fl[0]) if len(fl) == 1 else None
            if count is None or count < 0:
                break
            c.next("an outline count")
            for _ in range(count):
                c.next("an outline point")


def _read_groups(c, f):
    try:
        while not c.at_end():
            head = c.fields("a group")
            gid = c.int(head, 0, "group id")
            title = _title(c.next("a group title"))
            nodes, elements = [], []
            for _ in range(3 + 18):
                c.next("a group record")
            c.fields("a rule count")
            while True:
                if c.int(c.fields("a group rule"), 0, "rule type") == -1:
                    break
                while True:
                    if c.int(c.fields("a group rule entry"), 0, "rule entry") == -1:
                        break
            c.fields("a list count")
            while True:
                ltype = c.int(c.fields("a group list"), 0, "list type")
                if ltype == -1:
                    break
                while True:
                    eid = c.int(c.fields("a group list entry"), 0, "list entry")
                    if eid == -1:
                        break
                    if ltype == 7:
                        nodes.append(eid)
                    elif ltype == 8:
                        elements.append(eid)
            f.groups.append((gid, title, nodes, elements))
    except ReadError as e:
        warn(f"Femap neutral: groups (block 408) past the ones read are skipped: {e}")


def _is_int_line(line, lo, hi, positive):
    fl = _fields(line)
    if not lo <= len(fl) <= hi:
        return False
    for s in fl:
        v = _parse_int(s)
        if v is None or (positive and v <= 0):
            return False
    return True


def _is_real_line(line):
    fl = _fields(line)
    return len(fl) == 1 and _parse_real(fl[0]) is not None


def _set_starts(c, ahead=0):
    if c.remaining() < ahead + 5:
        return False
    return (
        _is_int_line(c.peek(ahead), 1, 1, True)
        and _is_int_line(c.peek(ahead + 2), 2, 4, False)
        and _is_real_line(c.peek(ahead + 3))
        and _is_int_line(c.peek(ahead + 4), 1, 1, False)
    )


def _read_sets(c, f):
    while not c.at_end():
        if not _set_starts(c):
            c.fail("expected an output set record")
        sid = c.int(c.fields("an output set"), 0, "output set id")
        title = _title(c.next("an output set title"))
        c.next("an output set record")
        value = c.real(c.fields("an output set value"), 0, "output set value")
        notes = c.int(c.fields("a note count"), 0, "note count")
        for _ in range(notes):
            c.next("an output set note")
        while not c.at_end() and not _set_starts(c):
            c.next("an output set record")
        f.sets.append((sid, title, value))


def _read_vectors(c, f, ranges):
    while not c.at_end():
        head = c.fields("an output vector")
        sid = c.int(head, 0, "output set id")
        vid = c.int(head, 1, "output vector id")
        title = _title(c.next("an output vector title"))
        c.next("an output vector range")
        c.next("output vector components")
        c.next("output vector components")
        fl = c.fields("an output vector record")
        if len(fl) == 1:
            fl = c.fields("an output vector record")
        entity = c.int(fl, 3, "output vector entity type")
        c.next("an output vector record")
        values = []
        while True:
            d = c.fields("output vector data")
            if d and d[0] == "-1":
                break
            # A 1051 record is a range when its second field is an integer (the
            # range end); Femap 11 also writes plain ``id,value`` records there.
            if not ranges or len(d) < 2 or _parse_int(d[1]) is None:
                values.append((c.int(d, 0, "entity id"), c.real(d, 1, "value")))
                continue
            start = c.int(d, 0, "range start")
            stop = c.int(d, 1, "range end")
            if stop < start:
                c.fail(f"range {start}..{stop}")
            k = 2
            for eid in range(start, stop + 1):
                if k >= len(d):
                    d = c.fields("output vector data")
                    k = 0
                values.append((eid, c.real(d, k, "value")))
                k += 1
        f.vectors.append((sid, vid, title, entity, values))


def _parse(filename):
    with open_file(filename, "rb") as fh:
        raw = fh.read()
    text = raw.decode("latin-1") if isinstance(raw, bytes) else raw
    blocks = _blocks(text)
    f = _File()
    for bid, first, body in blocks:
        if bid == 100 and len(body) >= 2:
            c = _BlockCursor(bid, first, body)
            c.next("a title")
            f.version = c.real(c.fields("a version"), 0, "version")
    warned = set()
    readers = {
        403: _read_nodes,
        404: lambda c, f: _read_elements(c, f, warned),
        402: _read_properties,
        408: _read_groups,
        450: _read_sets,
        451: lambda c, f: _read_vectors(c, f, False),
        1051: lambda c, f: _read_vectors(c, f, True),
    }
    for bid, first, body in blocks:
        if bid in readers:
            readers[bid](_BlockCursor(bid, first, body), f)
    return f


def _resolve_step(time_step, n):
    k = time_step + n if time_step < 0 else time_step
    if n == 0:
        if time_step not in (0, -1):
            raise ReadError(
                f"time_step {time_step} requested but the file has no steps"
            )
        return None
    if not 0 <= k < n:
        raise ReadError(f"time_step {time_step} out of range: the file has {n} step(s)")
    return k


def time_values(filename):
    return [value for _, _, value in _parse(filename).sets]


def read(filename, points_only=False, arrays=None, time_step=0):
    f = _parse(filename)
    if not f.node_ids:
        raise ReadError(
            f"Femap neutral: '{filename}' holds no nodes (block 403): a geometry-only or "
            "results-only file has no mesh to read"
        )
    node_index = {}
    for p, nid in enumerate(f.node_ids):
        if nid in node_index:
            raise ReadError(f"Femap neutral: node {nid} is defined twice")
        node_index[nid] = p
    npts = len(f.node_ids)
    points = np.array(f.coords, dtype=np.float64).reshape(npts, 3)

    order = []
    by_topology = {}
    for e, el in enumerate(f.elements):
        if el[3] not in by_topology:
            by_topology[el[3]] = []
            order.append(el[3])
        by_topology[el[3]].append(e)
    cells = []
    prop_blocks, type_blocks = [], []
    element_index = {}
    cell_property = []
    cell_dim = []
    starts = [0]
    for code in order:
        cell_type, slots, _ = _TOPOLOGIES[code]
        dim = topological_dimension[cell_type]
        members = by_topology[code]
        conn = np.empty((len(members), len(slots)), dtype=np.int64)
        props = np.empty(len(members), dtype=np.int64)
        types = np.empty(len(members), dtype=np.int64)
        for r, e in enumerate(members):
            eid, prop, etype, _, el_slots, line = f.elements[e]
            for j, s in enumerate(slots):
                nid = el_slots[s]
                idx = node_index.get(nid)
                if idx is None:
                    what = (
                        f" has no node in slot {s} its topology needs"
                        if nid == 0
                        else f" names undefined node {nid}"
                    )
                    raise ReadError(f"Femap neutral: element {eid}{what} (line {line})")
                conn[r, j] = idx
            if eid in element_index:
                raise ReadError(
                    f"Femap neutral: element {eid} is defined twice (line {line})"
                )
            element_index[eid] = len(cell_property)
            props[r] = prop
            types[r] = etype
            cell_property.append(prop)
            cell_dim.append(dim)
        cells.append((cell_type, conn))
        prop_blocks.append(props)
        type_blocks.append(types)
        starts.append(starts[-1] + len(members))
    mesh = Mesh(points, cells)
    if prop_blocks:
        mesh.cell_data["femap:property"] = prop_blocks
        mesh.cell_data["femap:type"] = type_blocks
    ncells = len(cell_property)

    regions = []
    by_property = {}
    for g, pid in enumerate(cell_property):
        by_property.setdefault(pid, []).append(g)
    for pid in sorted(by_property):
        ids = by_property[pid]
        name = f.properties.get(pid) or f"property_{pid}"
        dim = max(cell_dim[g] for g in ids)
        regions.append(Region(name, "cell", np.array(ids, dtype=np.int64), dim, pid))
    missing = 0
    for gid, title, nodes, elements in f.groups:
        name = title or f"group_{gid}"
        pts, cls = [], []
        dim = -1
        for nid in nodes:
            idx = node_index.get(nid)
            if idx is None:
                missing += 1
            else:
                pts.append(idx)
        for eid in elements:
            idx = element_index.get(eid)
            if idx is None:
                missing += 1
            else:
                cls.append(idx)
                dim = max(dim, cell_dim[idx])
        if cls or not pts:
            regions.append(
                Region(name, "cell", np.array(cls, dtype=np.int64), dim, gid)
            )
        if pts:
            regions.append(
                Region(name, "point", np.array(pts, dtype=np.int64), -1, gid)
            )
    if missing:
        warn(
            f"Femap neutral: groups name {missing} undefined or skipped node(s) or "
            "element(s)"
        )
    mesh.regions = regions
    mesh.time_values = [value for _, _, value in f.sets]

    if not f.sets:
        if f.vectors:
            warn(
                "Femap neutral: output vectors without an output set (block 450) are "
                "skipped"
            )
        _resolve_step(time_step, 0)
        return mesh
    sid, _, value = f.sets[_resolve_step(time_step, len(f.sets))]
    mesh.field_data["meshio:time"] = np.array(value, dtype=np.float64)
    mesh.field_data["femap:set"] = np.array(sid, dtype=np.int64)
    if points_only or (arrays is not None and len(arrays) == 0):
        return mesh
    used = set()
    skipped = 0
    for vsid, vid, title, entity, values in f.vectors:
        if vsid != sid:
            continue
        if entity not in (7, 8):
            skipped += 1
            continue
        name = title or f"vector_{vid}"
        if name in used:
            name += f" ({vid})"
        used.add(name)
        if arrays is not None and name not in arrays:
            continue
        if entity == 7:
            data = np.full(npts, np.nan)
            for eid, v in values:
                idx = node_index.get(eid)
                if idx is not None:
                    data[idx] = v
            mesh.point_data[name] = data
        else:
            per_cell = np.full(ncells, np.nan)
            for eid, v in values:
                idx = element_index.get(eid)
                if idx is not None:
                    per_cell[idx] = v
            mesh.cell_data[name] = [
                per_cell[starts[b] : starts[b + 1]].copy() for b in range(len(order))
            ]
    if skipped:
        warn(
            f"Femap neutral: {skipped} output vector(s) on neither nodes nor elements "
            "skipped"
        )
    return mesh


def _fmt(v):
    return "%.17g" % float(v)


def _clean_title(title):
    t = title.replace("\n", " ").replace("\r", " ")
    return t or "<NULL>"


class _Plan:
    """What the writer makes of a mesh: each cell's topology, property, type
    and label (0 when it is dropped), and the regions in their canonical order."""

    def __init__(self, mesh):
        points = np.asarray(mesh.points)
        self.pdim = points.shape[1] if points.ndim == 2 else 0
        if self.pdim > 3:
            raise WriteError(
                f"Femap neutral writer: points of dimension {self.pdim} (at most 3)"
            )
        self.points = points
        self.npts = len(points)
        self.tops = []
        self.starts = [0]
        dropped = set()
        for block in mesh.cells:
            top = (
                _TOPOLOGY_OF.get(block.type)
                if isinstance(block.data, np.ndarray)
                else None
            )
            self.tops.append(top)
            if top is None and len(block.data):
                dropped.add(block.type)
            self.starts.append(self.starts[-1] + len(block.data))
        self.ncells = self.starts[-1]
        pdata = mesh.cell_data.get("femap:property")
        tdata = mesh.cell_data.get("femap:type")
        self.prop = [1] * self.ncells
        self.etype = [0] * self.ncells
        self.label = [0] * self.ncells
        self.written = 0
        for b, block in enumerate(mesh.cells):
            for r in range(len(block.data)):
                g = self.starts[b] + r
                if pdata is not None:
                    self.prop[g] = int(np.asarray(pdata[b]).ravel()[r])
                if tdata is not None:
                    self.etype[g] = int(np.asarray(tdata[b]).ravel()[r])
                elif self.tops[b] is not None:
                    self.etype[g] = self.tops[b][2]
                if self.tops[b] is not None:
                    self.written += 1
                    self.label[g] = self.written

        for t in sorted(dropped):
            warn(f"Femap neutral has no '{t}' topology; those cells are dropped")
            _provenance.note("cells-dropped", f"Femap neutral has no '{t}' topology")
        self.regions = sorted(getattr(mesh, "regions", []) or [], key=lambda r: r.key)
        side_regions = sum(1 for r in self.regions if r.kind == "side")
        if side_regions:
            warn(
                f"Femap neutral has no facet groups; {side_regions} side region(s) "
                "dropped"
            )
            _provenance.note(
                "regions-dropped", f"{side_regions} side region(s) have no Femap group"
            )


def _vectors(mesh, plan):
    """The output vectors of a mesh's data, ``(title, entity, [(id, value)])``: a
    point array per component as nodal vectors (entity 7), a cell array per
    component as elemental ones (8); and the arrays that have none."""
    npts, starts, label = plan.npts, plan.starts, plan.label
    vectors = []
    unwritable = []
    for name in sorted(mesh.point_data):
        a = np.asarray(mesh.point_data[name])
        if a.ndim not in (1, 2) or len(a) != npts or a.dtype.kind not in "biuf":
            unwritable.append(name)
            continue
        a = a.reshape(npts, -1).astype(np.float64)
        nc = a.shape[1]
        for c in range(nc):
            values = [
                (p + 1, float(a[p, c])) for p in range(npts) if not np.isnan(a[p, c])
            ]
            vectors.append((name if nc == 1 else f"{name}_{c}", 7, values))
    for name in sorted(mesh.cell_data):
        if name.startswith("femap:"):
            continue
        arrays = [np.asarray(x) for x in mesh.cell_data[name]]
        widths = {x.shape[1] if x.ndim == 2 else 1 for x in arrays}
        ok = (
            len(arrays) == len(mesh.cells)
            and len(widths) == 1
            and all(
                x.ndim in (1, 2)
                and len(x) == len(block.data)
                and x.dtype.kind in "biuf"
                for x, block in zip(arrays, mesh.cells)
            )
        )
        if not ok:
            unwritable.append(name)
            continue
        nc = widths.pop()
        for c in range(nc):
            values = []
            for b, x in enumerate(arrays):
                x = x.reshape(len(x), -1).astype(np.float64)
                for r in range(len(x)):
                    g = starts[b] + r
                    if label[g] and not np.isnan(x[r, c]):
                        values.append((label[g], float(x[r, c])))
            vectors.append((name if nc == 1 else f"{name}_{c}", 8, values))
    return vectors, unwritable


def _set_fields(mesh, unwritable):
    """The output set's ``femap:set`` and ``meshio:time`` (None when absent);
    other field data goes to ``unwritable``."""
    set_id = set_value = None
    for name in sorted(mesh.field_data):
        a = np.asarray(mesh.field_data[name]).ravel()
        if name == "femap:set" and a.size == 1:
            set_id = int(a[0])
        elif name == "meshio:time" and a.size == 1:
            set_value = float(a[0])
        else:
            unwritable.append(name)
    return set_id, set_value


def _note_unwritable(unwritable):
    if unwritable:
        listed = ", ".join(unwritable)
        warn(
            f"Femap neutral writer: arrays with no Femap output vector dropped: {listed}"
        )
        _provenance.note(
            "data-dropped", f"arrays with no Femap output vector: {listed}"
        )


def _block_open(out, bid):
    out.append(f"   -1\n{bid:6d}\n")


def _block_close(out):
    out.append("   -1\n")


def _mesh_blocks(mesh, plan):
    """Blocks 100 (header), 402 (properties), 403 (nodes), 404 (elements) and
    408 (groups)."""
    ncells, label, prop, etype = plan.ncells, plan.label, plan.prop, plan.etype
    regions = plan.regions
    prop_cells = {}
    prop_type = {}
    for g in range(ncells):
        if label[g]:
            prop_cells.setdefault(prop[g], []).append(g)
            prop_type.setdefault(prop[g], etype[g])
    prop_title = {}
    property_regions = set()
    for k, reg in enumerate(regions):
        if reg.kind != "cell" or reg.tag not in prop_cells or reg.tag in prop_title:
            continue
        if [int(v) for v in np.asarray(reg.entries).ravel()] == prop_cells[reg.tag]:
            prop_title[reg.tag] = reg.name
            property_regions.add(k)

    out = []
    _block_open(out, 100)
    out.append(
        _clean_title(_provenance.lines(_provenance.SlotTier.SINGLE_LINE)[0])
        + "\n8.2,\n"
    )
    _block_close(out)

    def zero_lines(count, per_line, zero):
        for k in range(0, count, per_line):
            out.append(f"{zero}," * (min(count, k + per_line) - k) + "\n")

    if prop_cells:
        _block_open(out, 402)
        for pid in sorted(prop_cells):
            out.append(f"{pid},110,0,{prop_type[pid]},1,0,\n")
            out.append(
                _clean_title(prop_title.get(pid, f"property_{pid}"))
                + "\n0,0,0,0,\n90,\n"
            )
            zero_lines(90, 8, "0")
            out.append("190,\n")
            zero_lines(190, 5, "0.")
            out.append("0,\n0,\n")
        _block_close(out)

    _block_open(out, 403)
    pts = plan.points.astype(np.float64, copy=False)
    for p in range(plan.npts):
        xyz = [pts[p, d] if d < plan.pdim else 0.0 for d in range(3)]
        out.append(
            f"{p + 1},0,0,1,46,0,0,0,0,0,0,"
            + "".join(f"{_fmt(v)}," for v in xyz)
            + "0,\n"
        )
    _block_close(out)

    if plan.written:
        _block_open(out, 404)
        for b, block in enumerate(mesh.cells):
            if plan.tops[b] is None:
                continue
            code, slots, _ = plan.tops[b]
            data = np.asarray(block.data, dtype=np.int64)
            for r in range(len(data)):
                g = plan.starts[b] + r
                sl = [0] * 20
                for j, s in enumerate(slots):
                    sl[s] = int(data[r, j]) + 1
                out.append(
                    f"{label[g]},124,{prop[g]},{etype[g]},{code},1,0,0,0,0,0,0,0,\n"
                )
                out.append("".join(f"{v}," for v in sl[:10]) + "\n")
                out.append("".join(f"{v}," for v in sl[10:]) + "\n")
                out.append(
                    "0.,0.,0.,\n0.,0.,0.,\n0.,0.,0.,\n0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,\n"
                )
        _block_close(out)

    names = []
    groups = {}
    tags = {}
    for k, reg in enumerate(regions):
        if reg.kind == "side" or k in property_regions:
            continue
        if reg.name not in groups:
            groups[reg.name] = ([], [])
            names.append(reg.name)
            tags[reg.name] = reg.tag
        nodes, elements = groups[reg.name]
        for e in np.asarray(reg.entries, dtype=np.int64).ravel():
            e = int(e)
            if reg.kind == "point":
                nodes.append(e + 1)
            elif e < ncells and label[e]:
                elements.append(label[e])
    if names:
        used = {tags[n] for n in names if tags[n] > 0}
        assigned = set()
        nxt = 1
        _block_open(out, 408)
        for name in names:
            gid = tags[name]
            if gid <= 0 or gid in assigned:
                while nxt in used or nxt in assigned:
                    nxt += 1
                gid = nxt
            assigned.add(gid)
            out.append(f"{gid},0,0,\n{_clean_title(name)}\n")
            out.append("0,0,0,\n0,0,0,0,0.,0.,\n0,0,\n")
            out.append("0,0,\n0.,0.,0.,\n0.,0.,0.,\n" * 6)
            out.append("91,\n-1,\n23,\n")
            nodes, elements = groups[name]
            if nodes:
                out.append("7,\n" + "".join(f"{v},\n" for v in nodes) + "-1,\n")
            if elements:
                out.append("8,\n" + "".join(f"{v},\n" for v in elements) + "-1,\n")
            out.append("-1,\n")
        _block_close(out)
    return out


def _set_blocks(set_id, set_value, vectors, ids):
    """One output set (450) and its vectors (451); ``ids`` gives each vector's id."""
    out = []
    _block_open(out, 450)
    out.append(f"{set_id},\nmeshio++\n0,1,\n{_fmt(set_value)},\n0,\n")
    _block_close(out)
    if not vectors:
        return out
    _block_open(out, 451)
    for (title, entity, values), vid_ in zip(vectors, ids):
        lo = hi = absmax = 0.0
        id_lo = id_hi = 0
        for j, (vid, x) in enumerate(values):
            if j == 0 or x < lo:
                lo, id_lo = x, vid
            if j == 0 or x > hi:
                hi, id_hi = x, vid
            absmax = max(absmax, abs(x))
        out.append(f"{set_id},{vid_},1,\n{_clean_title(title)}\n")
        out.append(f"{_fmt(lo)},{_fmt(hi)},{_fmt(absmax)},\n")
        out.append("0,0,0,0,0,0,0,0,0,0,\n" * 2)
        out.append(f"{id_lo},{id_hi},0,{entity},\n0,0,1,\n")
        out.append("".join(f"{vid},{_fmt(x)},\n" for vid, x in values))
        out.append("-1,0.,\n")
    _block_close(out)
    return out


def write(filename, mesh):
    plan = _Plan(mesh)
    vectors, unwritable = _vectors(mesh, plan)
    set_id, set_value = _set_fields(mesh, unwritable)
    _note_unwritable(unwritable)
    out = _mesh_blocks(mesh, plan)
    if vectors:
        out += _set_blocks(
            max(1, set_id or 0),
            0.0 if set_value is None else set_value,
            vectors,
            range(1, len(vectors) + 1),
        )
    with open_file(filename, "w", newline="\n") as fh:
        fh.write("".join(out))


def _fingerprint(mesh):
    """What must not change between the steps of a series, as digests (no copy
    of a step is kept): the cells, and separately the points."""
    import hashlib

    cells = hashlib.blake2b(digest_size=16)
    for block in mesh.cells:
        data = np.asarray(block.data, dtype=np.int64)
        cells.update(f"{block.type}:{data.shape};".encode())
        cells.update(np.ascontiguousarray(data).tobytes())
    points = hashlib.blake2b(
        np.ascontiguousarray(np.asarray(mesh.points, dtype=np.float64)).tobytes(),
        digest_size=16,
    )
    return cells.digest(), points.digest()


class SeriesWriter:
    """Write a time series into one neutral file: the mesh once, then an output
    set (450) with its vectors (451) per step, the steps ``read(time_step=k)``
    reads back. Every step must have the first step's cells; a step whose points
    moved is written with the first step's, with a warning. Use as a context
    manager, or call :meth:`close`."""

    def __init__(self, filename):
        # A path is opened as `write` opens it; a caller's buffer is written
        # into and left open.
        self._owned = not is_buffer(filename, "w")
        self._fh = open(filename, "w", newline="\n") if self._owned else filename
        self._cells = self._points = None
        self._plan = None
        self._set_ids = set()
        self._next_set = 1
        self._vector_ids = {}
        self._moved = False
        self._step = 0

    def write(self, time, mesh):
        if self._plan is None:
            self._plan = _Plan(mesh)
            vectors, unwritable = _vectors(mesh, self._plan)
            set_id, _ = _set_fields(mesh, unwritable)
            _note_unwritable(unwritable)
            self._fh.write("".join(_mesh_blocks(mesh, self._plan)))
            self._cells, self._points = _fingerprint(mesh)
        else:
            cells, points = _fingerprint(mesh)
            if cells != self._cells:
                raise WriteError(
                    f"Femap series: step {self._step}'s cells differ from the first "
                    "step's; a neutral file holds one mesh -- write one file per step "
                    "with '{step}'"
                )
            if points != self._points and not self._moved:
                self._moved = True
                warn(
                    "Femap series: points moved after the first step; written as the first step's"
                )
                _provenance.note(
                    "points-moved", "every output set shares the first step's points"
                )
            vectors, unwritable = _vectors(mesh, self._plan)
            set_id, _ = _set_fields(mesh, unwritable)
        if set_id is None or set_id <= 0 or set_id in self._set_ids:
            while self._next_set in self._set_ids:
                self._next_set += 1
            set_id = self._next_set
        self._set_ids.add(set_id)
        ids = []
        for title, entity, _ in vectors:
            key = (title, entity)
            if key not in self._vector_ids:
                self._vector_ids[key] = len(self._vector_ids) + 1
            ids.append(self._vector_ids[key])
        self._fh.write("".join(_set_blocks(set_id, float(time), vectors, ids)))
        self._step += 1

    def close(self):
        if self._fh is not None:
            if self._owned:
                self._fh.close()
            self._fh = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
