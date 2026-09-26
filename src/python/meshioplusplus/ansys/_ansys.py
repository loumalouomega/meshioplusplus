"""
I/O for ANSYS Fluent's ``.msh`` format (also TGrid and GAMBIT meshes).

<https://romeo.univ-reims.fr/documents/fluent/tgrid/ug/appb.pdf>

A Fluent mesh stores points, faces and cell *zones*: a face row lists its nodes
and the two cells it separates, ``n0 .. nk c0 c1`` (hexadecimal in ASCII), and a
cell zone normally declares only a range of cell ids. Volume cells are rebuilt
from their faces: the right-hand normal of ``n0 .. nk`` points into ``c0`` (in
2-D, walking ``n0 -> n1`` leaves ``c0`` on the left), so a face is outward for
``c1`` and reversed for ``c0``. Boundary faces are kept as surface cells, every
cell carries its zone id in ``cell_data["ansys:zone"]``, and every zone is a
cell region named from its ``(39 ...)``/``(45 ...)`` declaration.

A cell section with a connectivity body is meshio's old (legacy) layout, which
Fluent itself cannot read; such a file reads as before, cells only.
:func:`write` writes faces: every face once, with ``c0``/``c1``, zone by zone
(the ``ansys:zone`` values, else one per block), named from the regions.
"""

import re
from collections import defaultdict

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._face_cells import polygon_from_edges, reconstruct_cell
from .._files import open_file
from .._mesh import CellBlock, Mesh
from .._regions import Region
from .._skin import _CELL_FACES

# Cell element types (zone header field 4) with a fixed node count.
_CELL_TYPES = {
    1: ("triangle", 3),
    2: ("tetra", 4),
    3: ("quad", 4),
    4: ("hexahedron", 8),
    5: ("pyramid", 5),
    6: ("wedge", 6),
}
# Face element types with a fixed node count; 0 (mixed) and 5 (polygonal) rows
# lead with their node count.
_FACE_NODES = {2: 2, 3: 3, 4: 4}
_INTERIOR = 2  # the face-zone type of interior faces

_SECTION = re.compile(rb"\(\s*(\d+)")
_ZONE_NAME = re.compile(rb"\(\s*(?:39|45)\s*\(\s*(\d+)\s+([^\s()]+)\s+([^\s()]+)")


class _Reader:
    def __init__(self, data):
        self.d = data
        self.p = 0

    def skip_ws(self):
        d, n = self.d, len(self.d)
        while self.p < n and d[self.p] in b" \t\r\n":
            self.p += 1

    def skip_balanced(self, depth):
        """Advance past ``depth`` unmatched closing brackets; quoted text is
        skipped whole."""
        d, n = self.d, len(self.d)
        while depth > 0 and self.p < n:
            c = d[self.p]
            self.p += 1
            if c == 0x28:  # (
                depth += 1
            elif c == 0x29:  # )
                depth -= 1
            elif c == 0x22:  # "
                q = d.find(b'"', self.p)
                self.p = n if q < 0 else q + 1

    def skip_binary_trailer(self):
        """Past ``End of Binary Section NNNN)`` if it follows, else past the
        section's own closing bracket."""
        self.skip_ws()
        if self.d.startswith(b"End of Binary Section", self.p):
            q = self.d.find(b")", self.p)
            self.p = len(self.d) if q < 0 else q + 1
        elif self.p < len(self.d) and self.d[self.p] == 0x29:
            self.p += 1

    def header(self):
        """The ``(a b c ...)`` group after a section index, as hex integers,
        and whether a body follows."""
        self.skip_ws()
        if self.p >= len(self.d) or self.d[self.p] != 0x28:
            raise ReadError("Fluent: expected a section header")
        q = self.d.find(b")", self.p)
        if q < 0:
            raise ReadError("Fluent: unterminated section header")
        values = [int(t, 16) for t in self.d[self.p + 1 : q].split()]
        self.p = q + 1
        # A body opens with '(' after nothing but blanks; anything else (a
        # closing bracket, `End of Binary Section`) makes a declaration.
        start = self.p
        self.skip_ws()
        if self.p < len(self.d) and self.d[self.p] == 0x28:
            newline_before = b"\n" in self.d[start : self.p]
            self.p += 1
            # meshio's writer puts the body's '(' last on the header line and
            # starts the data on the next line.
            if not newline_before:
                if self.d.startswith(b"\r\n", self.p):
                    self.p += 2
                elif self.d.startswith(b"\n", self.p):
                    self.p += 1
            return values, True
        self.p = start
        self.skip_binary_trailer()
        return values, False

    def ascii_body(self):
        q = self.d.find(b")", self.p)
        if q < 0:
            raise ReadError("Fluent: unterminated section body")
        tokens = self.d[self.p : q].split()
        self.p = q + 1
        self.skip_ws()
        if self.p < len(self.d) and self.d[self.p] == 0x29:
            self.p += 1
        return tokens

    def binary(self, count, dtype):
        size = count * np.dtype(dtype).itemsize
        if self.p + size > len(self.d):
            raise ReadError("Fluent: binary section runs past the end of the file")
        out = np.frombuffer(self.d, dtype=dtype, count=count, offset=self.p)
        self.p += size
        # the body's ')' then `End of Binary Section NNNN)` or the section's ')'
        self.skip_ws()
        if self.p < len(self.d) and self.d[self.p] == 0x29:
            self.p += 1
        self.skip_binary_trailer()
        return out


def _int_dtype(prefix):
    return np.dtype("<i4") if prefix == "20" else np.dtype("<i8")


def _face_rows(values, count, face_type):
    """Split a face body into (nodes, c0, c1) rows."""
    k = _FACE_NODES.get(face_type)
    if k is not None:
        rows = np.asarray(values).reshape(count, k + 2)
        return [list(r[:k]) for r in rows], rows[:, k], rows[:, k + 1]
    nodes, c0, c1 = [], [], []
    i = 0
    for _ in range(count):
        n = int(values[i])
        nodes.append(list(values[i + 1 : i + 1 + n]))
        c0.append(values[i + 1 + n])
        c1.append(values[i + 2 + n])
        i += n + 3
    return nodes, np.asarray(c0), np.asarray(c1)


def _binary_mixed_faces(rd, count, dtype):
    """A mixed binary face body: each row leads with its node count."""
    item = np.dtype(dtype).itemsize
    d = rd.d
    values = []
    for _ in range(count):
        if rd.p + item > len(d):
            raise ReadError("Fluent: binary section runs past the end of the file")
        n = int(np.frombuffer(d, dtype=dtype, count=1, offset=rd.p)[0])
        row = np.frombuffer(d, dtype=dtype, count=n + 3, offset=rd.p)
        values.extend(int(v) for v in row)
        rd.p += (n + 3) * item
    rd.skip_ws()
    if rd.p < len(d) and d[rd.p] == 0x29:
        rd.p += 1
    rd.skip_binary_trailer()
    return values


def _refuse_early(head):
    """Refuse a file that does not open with a Fluent section from its first
    bytes, with the message the full read would give -- before reading it
    whole. Every ``.msh`` (Gmsh, FreeFEM) is offered to this reader first."""
    p, n = 0, len(head)
    while p < n and head[p] in b" \t\r\n":
        p += 1
    if p == n:
        return  # undecided: all blanks so far
    q = p + 1
    while q < n and head[q] in b" \t\r\n":
        q += 1
    if head[p] != 0x28 or (q < n and not 0x30 <= head[q] <= 0x39):
        raise ReadError(f"Fluent: expected a section at byte {p}")


def read(filename):  # noqa: C901
    with open_file(filename, "rb") as f:
        head = f.read(256)
        if isinstance(head, str):
            head = head.encode()
        _refuse_early(head)
        if f.seekable():
            f.seek(0)
            data = f.read()
        else:
            rest = f.read()
            data = head + (rest.encode() if isinstance(rest, str) else rest)
    if isinstance(data, str):
        data = data.encode()
    rd = _Reader(data)

    dim = None
    node_zones = []  # (first, last, points)
    faces = []  # (zone, bc_type, nodes, c0, c1) per face, file order
    cell_zones = []  # (zone, first, last, zone_type, element_type)
    legacy = []  # (meshio type, connectivity) of cell sections with bodies
    names = {}  # zone id -> name

    while True:
        rd.skip_ws()
        if rd.p >= len(data):
            break
        m = _SECTION.match(data, rd.p)
        if m is None:
            raise ReadError(f"Fluent: expected a section at byte {rd.p}")
        index = m.group(1).decode()
        start = rd.p
        rd.p = m.end()
        prefix, core = ("", index) if len(index) <= 2 else (index[:-2], index[-2:])

        if index == "2":
            rd.skip_ws()
            q = data.find(b")", rd.p)
            dim = int(data[rd.p : q].split()[0])
            rd.p = q + 1
        elif index in ("39", "45"):
            z = _ZONE_NAME.match(data, start)
            if z is not None:
                names[int(z.group(1))] = z.group(3).decode()
            rd.skip_balanced(1)
        elif core in ("10", "12", "13") and prefix in ("", "20", "30"):
            head, has_body = rd.header()
            if len(head) < 4:
                continue
            zone, first, last = head[0], head[1], head[2]
            count = last - first + 1
            if core == "10":
                if not has_body:
                    continue
                nd = head[4] if len(head) > 4 else (dim or 3)
                if prefix == "":
                    pts = np.array(rd.ascii_body(), dtype=float)
                else:
                    pts = rd.binary(count * nd, "<f4" if prefix == "20" else "<f8")
                if pts.size != count * nd:
                    raise ReadError("Fluent: node section size mismatch")
                node_zones.append((first, last, pts.reshape(count, nd).astype(float)))
            elif core == "13":
                if not has_body or zone == 0:
                    continue
                bc_type, face_type = head[3], head[4] if len(head) > 4 else 0
                if prefix == "":
                    values = [int(t, 16) for t in rd.ascii_body()]
                elif face_type in _FACE_NODES:
                    values = rd.binary(
                        count * (_FACE_NODES[face_type] + 2), _int_dtype(prefix)
                    )
                else:
                    values = _binary_mixed_faces(rd, count, _int_dtype(prefix))
                nodes, c0, c1 = _face_rows(values, count, face_type)
                for f, a, b in zip(nodes, c0, c1):
                    faces.append((zone, bc_type, f, int(a), int(b)))
            else:
                zone_type = head[3]
                element_type = head[4] if len(head) > 4 else 0
                body = None
                if has_body:
                    if prefix == "":
                        body = [int(t, 16) for t in rd.ascii_body()]
                    else:
                        npc = _CELL_TYPES.get(element_type, (None, 1))[1]
                        body = rd.binary(count * npc, _int_dtype(prefix))
                if zone == 0 or zone_type == 0:
                    continue  # the global declaration, or a dead zone
                if body is not None and element_type in _CELL_TYPES:
                    key, npc = _CELL_TYPES[element_type]
                    legacy.append(
                        (key, np.asarray(body, dtype=np.int64).reshape(count, npc))
                    )
                else:
                    cell_zones.append((zone, first, last, zone_type, element_type))
        else:
            if prefix in ("20", "30"):
                # An unknown binary section: its body is not bracket-balanced.
                q = data.find(b"End of Binary Section", rd.p)
                q = data.find(b")", q) if q >= 0 else -1
                rd.p = len(data) if q < 0 else q + 1
            else:
                rd.skip_balanced(1)

    if not node_zones:
        raise ReadError("Fluent: no nodes")
    base = min(z[0] for z in node_zones)
    top = max(z[1] for z in node_zones)
    nd = max(z[2].shape[1] for z in node_zones)
    points = np.zeros((top - base + 1, nd))
    for first, last, pts in node_zones:
        points[first - base : last - base + 1, : pts.shape[1]] = pts

    if legacy:
        if faces:
            warn(
                "Fluent: cells with connectivity bodies; the face sections are ignored"
            )
        return Mesh(points, [(k, c - base) for k, c in legacy])

    dim = dim or nd
    return _from_faces(points, base, dim, faces, cell_zones, names)


def _from_faces(points, base, dim, faces, cell_zones, names):  # noqa: C901
    zone_of = {}
    for zone, first, last, _, _ in cell_zones:
        for c in range(first, last + 1):
            zone_of[c] = zone
    if not cell_zones:
        # No cell zone declared: every referenced cell is in zone 0.
        for _, _, _, a, b in faces:
            for c in (a, b):
                if c:
                    zone_of.setdefault(c, 0)

    P = (
        points
        if points.shape[1] == 3
        else np.column_stack([points, np.zeros(len(points))])
    )
    # Outward faces (3-D) or directed edges (2-D) per cell, in file order.
    per_cell = defaultdict(list)
    for _, _, nodes, c0, c1 in faces:
        f = [n - base for n in nodes]
        if c1 in zone_of:
            per_cell[c1].append(f)
        if c0 in zone_of:
            per_cell[c0].append(f[::-1])

    volume = {}  # (zone, type) -> list of connectivities, first-seen order
    skipped = 0
    for cid in sorted(zone_of):
        cf = per_cell.get(cid)
        if not cf:
            skipped += 1
            continue
        if dim == 2:
            ring = polygon_from_edges([(f[0], f[-1]) for f in cf], P)
            if ring is None:
                skipped += 1
                continue
            t = {3: "triangle", 4: "quad"}.get(len(ring), f"polygon{len(ring)}")
            volume.setdefault((zone_of[cid], t), []).append(ring)
            continue
        t, conn = reconstruct_cell(cf, P)
        if conn is None:
            skipped += 1
            continue
        if t == "polyhedron":
            t = f"polyhedron{len(set().union(*conn))}"
            conn = [np.array(f, dtype=int) for f in conn]
        volume.setdefault((zone_of[cid], t), []).append(conn)
    if skipped:
        warn(f"Fluent: {skipped} cell(s) with no usable faces skipped")

    surface = {}
    for zone, bc_type, nodes, c0, c1 in faces:
        if bc_type == _INTERIOR and c0 in zone_of and c1 in zone_of:
            continue
        f = [n - base for n in nodes]
        # Outward from the domain: the normal points into c0.
        if c0 in zone_of and c1 not in zone_of:
            f = f[::-1]
        t = {2: "line", 3: "triangle", 4: "quad"}.get(len(f), f"polygon{len(f)}")
        surface.setdefault((zone, t), []).append(f)

    cells, zone_data, regions = [], [], defaultdict(list)
    offset = 0
    for group, cdim in ((volume, dim), (surface, dim - 1)):
        for (zone, t), conn in group.items():
            if t.startswith("polyhedron"):
                data = np.empty(len(conn), dtype=object)
                for i, c in enumerate(conn):
                    data[i] = c
            else:
                data = np.array(conn, dtype=int)
            cells.append(CellBlock(t, data))
            zone_data.append(np.full(len(conn), zone, dtype=int))
            regions[(zone, cdim)].append(np.arange(offset, offset + len(conn)))
            offset += len(conn)

    mesh = Mesh(points, cells, cell_data={"ansys:zone": zone_data} if cells else {})
    mesh.regions = [
        Region(names.get(zone, f"zone_{zone}"), "cell", np.concatenate(ids), cdim, zone)
        for (zone, cdim), ids in regions.items()
    ]
    return mesh


def _element_type(cell_type):
    """Fluent element type of a cell zone (12): 7 is a polyhedral (3-D) or
    polygonal (2-D) cell, defined by its faces."""
    for prefix, code in (
        ("triangle", 1),
        ("tetra", 2),
        ("quad", 3),
        ("hexahedron", 4),
        ("pyramid", 5),
        ("wedge", 6),
    ):
        if cell_type.startswith(prefix):
            return code
    return 7


_LINEAR = {
    "triangle",
    "tetra",
    "quad",
    "hexahedron",
    "pyramid",
    "wedge",
    "polygon",
    "line",
}


def _is_linear(cell_type):
    return (
        cell_type in _LINEAR
        or cell_type.startswith("polygon")
        or cell_type.startswith("polyhedron")
    )


def _block_dim(block):
    if block.type.startswith("polyhedron"):
        return 3
    if block.type.startswith("polygon"):
        return 2
    return block.dim


def _ring(block, cell):
    """The corner ring of a 2-D cell (a surface facet in 3-D, a cell in 2-D)."""
    row = [int(v) for v in block.data[cell]]
    t = block.type
    if t.startswith("triangle"):
        return row[:3]
    if t.startswith("quad"):
        return row[:4]
    if t.startswith("line"):
        return row[:2]
    return row


def _signed_volume(rings, points):
    """Signed volume enclosed by ``rings`` (positive when they wind outward)."""
    nodes = sorted({n for r in rings for n in r})
    c = points[nodes].mean(axis=0)
    vol = 0.0
    for r in rings:
        p = points[r] - c
        f = p.mean(axis=0)
        for i in range(len(r)):
            a, b = p[i], p[(i + 1) % len(r)]
            vol += float(np.dot(np.cross(a, b), f))
    return vol / 6.0


def _orient(rings, points):
    """``orient_rings``: wind the rings of one cell consistently (a BFS over
    shared edges from face 0), then all outward. Left as given when the rings
    are not a closed, orientable surface."""
    uses = defaultdict(list)
    for f, r in enumerate(rings):
        if len(r) < 3:
            return rings
        for i in range(len(r)):
            a, b = r[i], r[(i + 1) % len(r)]
            if a == b:
                return rings
            uses[(min(a, b), max(a, b))].append((f, a < b))
    if any(len(u) != 2 for u in uses.values()):
        return rings
    flip = [-1] * len(rings)
    flip[0] = 0
    stack = [0]
    visited = 0
    while stack:
        f = stack.pop()
        visited += 1
        r = rings[f]
        for i in range(len(r)):
            a, b = r[i], r[(i + 1) % len(r)]
            u = uses[(min(a, b), max(a, b))]
            mine = u[0] if u[0][0] == f and u[0][1] == (a < b) else u[1]
            other = u[1] if mine is u[0] else u[0]
            if other[0] == f:
                return rings
            mine_fwd = mine[1] != (flip[f] == 1)
            want = 1 if other[1] == mine_fwd else 0
            if flip[other[0]] == -1:
                flip[other[0]] = want
                stack.append(other[0])
            elif flip[other[0]] != want:
                return rings
    if visited != len(rings):
        return rings
    out = [r[::-1] if flip[k] == 1 else r for k, r in enumerate(rings)]
    if _signed_volume(out, points) < 0.0:
        out = [r[::-1] for r in out]
    return out


def _cell_rings(block, cell):
    """The outward (on the reference element) face rings of a volume cell."""
    if block.type.startswith("polyhedron"):
        return [[int(v) for v in face] for face in block.data[cell]]
    row = block.data[cell]
    return [
        [int(row[i]) for i in local[:ncorner]]
        for _, ncorner, local in _CELL_FACES[block.type]
    ]


def _name(name):
    name = "".join("_" if c in ' \t()"' else c for c in name)
    return name or "zone"


def write(filename, mesh, binary=True):
    """Write a Fluent mesh (``write_ansys`` in ansys.cpp, byte for byte).

    The cells are the blocks of the mesh's highest dimension (2 or 3), written
    as faces with ``c0``/``c1``; blocks one dimension lower name boundary zones;
    anything else is dropped with a warning.
    """
    points = np.asarray(mesh.points, dtype=np.float64)
    npoints, pdim = points.shape
    if pdim not in (2, 3):
        raise WriteError("Fluent: can only write points of dimension 2 or 3")
    dims = [_block_dim(b) for b in mesh.cells]
    dim = max(dims, default=0)
    if dim < 2:
        raise WriteError("Fluent: the mesh has no 2-D or 3-D cells")
    if dim == 3 and pdim != 3:
        raise WriteError("Fluent: 3-D cells need 3-D points")
    # a 2-D Fluent mesh lies in the xy plane: z is dropped only when it is 0
    if dim == 2 and pdim == 3:
        nonzero = np.flatnonzero(points[:, 2] != 0.0)
        if len(nonzero):
            raise WriteError(
                f"Fluent: a 2-D mesh must lie in the z = 0 plane (point {nonzero[0]} "
                "has z != 0); Fluent has no 3-D surface meshes"
            )
    bases = np.concatenate([[0], np.cumsum([len(b) for b in mesh.cells])]).astype(
        np.int64
    )
    zones_data = mesh.cell_data.get("ansys:zone")

    def zone_value(b, i):
        if zones_data is None or i >= len(zones_data[b]):
            return 0
        return max(int(zones_data[b][i]), 0)

    faces = []  # [nodes, c0, c1] with compact cell ids, -1 for none
    cell_to_global = []
    cell_block = []
    surface_blocks = []
    dropped = 0
    if dim == 3:
        seen = {}
        for b, block in enumerate(mesh.cells):
            volume = block.type.startswith("polyhedron") or block.type in _CELL_FACES
            if not volume:
                if dims[b] == 2:
                    surface_blocks.append(b)
                else:
                    dropped += 1
                continue
            for i in range(len(block)):
                c = len(cell_to_global)
                cell_to_global.append(int(bases[b]) + i)
                cell_block.append(b)
                for ring in _orient(_cell_rings(block, i), points):
                    key = tuple(sorted(ring))
                    f = seen.get(key)
                    if f is None:
                        seen[key] = len(faces)
                        # stored outward from the owner; Fluent's normal points into c0
                        faces.append([ring[::-1], c, -1])
                    elif faces[f][2] < 0:
                        faces[f][2] = c
    else:
        edge_of = {}
        for b, block in enumerate(mesh.cells):
            if dims[b] != 2:
                if dims[b] == 1:
                    surface_blocks.append(b)
                else:
                    dropped += 1
                continue
            for i in range(len(block)):
                ring = _ring(block, i)
                area = 0.0
                for k in range(len(ring)):
                    p, q = ring[k], ring[(k + 1) % len(ring)]
                    area += points[p, 0] * points[q, 1] - points[q, 0] * points[p, 1]
                if area < 0.0:
                    ring = ring[::-1]
                c = len(cell_to_global)
                cell_to_global.append(int(bases[b]) + i)
                cell_block.append(b)
                # counter-clockwise: an edge a -> b has the cell on its left (c0)
                for k in range(len(ring)):
                    a, e = ring[k], ring[(k + 1) % len(ring)]
                    key = (min(a, e), max(a, e))
                    f = edge_of.get(key)
                    if f is None:
                        edge_of[key] = len(faces)
                        faces.append([[a, e], c, -1])
                    elif faces[f][2] < 0:
                        faces[f][2] = c
    if not cell_to_global:
        raise WriteError("Fluent: the mesh has no cells Fluent can hold")
    if dropped:
        warn(f"Fluent: {dropped} cell block(s) of other dimensions are not written")
    quadratic = sum(
        1 for b in sorted(set(cell_block)) if not _is_linear(mesh.cells[b].type)
    )
    if quadratic:
        warn(
            f"Fluent: cells are linear; the mid-side nodes of {quadratic} block(s) are not used"
        )

    # zones: ansys:zone values are kept; the rest get fresh ids
    used = set()
    cell_zone_id = []
    for c, g in enumerate(cell_to_global):
        z = zone_value(cell_block[c], g - int(bases[cell_block[c]]))
        cell_zone_id.append(z)
        if z:
            used.add(z)
    cell_ids_used = set(used)
    surface_explicit = {}
    for b in surface_blocks:
        for i in range(len(mesh.cells[b])):
            z = zone_value(b, i)
            if z and z not in cell_ids_used:
                surface_explicit[(b, i)] = z
                used.add(z)
    state = {"next": max(used) + 1 if used else 1}

    def fresh():
        while state["next"] in used:
            state["next"] += 1
        used.add(state["next"])
        state["next"] += 1
        return state["next"] - 1

    block_zone = {}
    for c in range(len(cell_to_global)):
        if cell_zone_id[c] == 0:
            if cell_block[c] not in block_zone:
                block_zone[cell_block[c]] = fresh()
            cell_zone_id[c] = block_zone[cell_block[c]]

    cell_zones = []  # [id, type, global cells, dim, members]
    pos = {}
    for c, z in enumerate(cell_zone_id):
        if z not in pos:
            pos[z] = len(cell_zones)
            cell_zones.append([z, "fluid", [], dim, []])
        cell_zones[pos[z]][4].append(c)
        cell_zones[pos[z]][2].append(cell_to_global[c])
    fluent_cell = [0] * len(cell_to_global)
    n = 1
    for zone in cell_zones:
        for c in zone[4]:
            fluent_cell[c] = n
            n += 1

    face_zone = [0] * len(faces)
    face_zones = [[fresh(), "interior", [], dim - 1, []]]
    face_zone_pos = {}
    by_key = {}
    for f, face in enumerate(faces):
        if face[2] < 0:
            by_key.setdefault(tuple(sorted(face[0])), f)
    unmatched = 0
    for b in surface_blocks:
        block_id = 0
        for i in range(len(mesh.cells[b])):
            f = by_key.get(tuple(sorted(_ring(mesh.cells[b], i))))
            if f is None:
                unmatched += 1
                continue
            if face_zone[f]:
                continue
            z = surface_explicit.get((b, i))
            if z is None:
                if not block_id:
                    block_id = fresh()
                z = block_id
            face_zone[f] = z
            if z not in face_zone_pos:
                face_zone_pos[z] = len(face_zones)
                face_zones.append([z, "wall", [], dim - 1, []])
            face_zones[face_zone_pos[z]][2].append(int(bases[b]) + i)
    if unmatched:
        warn(
            f"Fluent: {unmatched} facet cell(s) are not on the boundary of the cells "
            "and are not written"
        )
    default_wall = 0
    for f, face in enumerate(faces):
        if face[2] >= 0:
            face_zones[0][4].append(f)
            continue
        if not face_zone[f]:
            if not default_wall:
                default_wall = fresh()
                face_zone_pos[default_wall] = len(face_zones)
                face_zones.append([default_wall, "wall", [], dim - 1, []])
            face_zone[f] = default_wall
        face_zones[face_zone_pos[face_zone[f]]][4].append(f)
    node_zone = fresh()

    def zone_name(zone):
        zid, ztype, members, zdim, _ = zone
        for r in mesh.regions:
            if r.kind == "cell" and r.tag == zid and r.dim == zdim:
                return _name(r.name)
        want = sorted(members)
        if want:
            for r in mesh.regions:
                if r.kind == "cell" and len(r.entries) == len(want):
                    if np.array_equal(np.asarray(r.entries), want):
                        return _name(r.name)
        return f"{ztype}_{zid}"

    other_cell_data = [k for k in mesh.cell_data if k != "ansys:zone"]
    if mesh.point_data or other_cell_data or mesh.field_data:
        _provenance.note("data-dropped", "a Fluent mesh file holds no data arrays")

    def hx(v):
        return f"{v:x}"

    with open_file(filename, "wb") as fh:
        out = []
        out.append(f'(1 "{_provenance.lines(_provenance.SlotTier.SINGLE_LINE)[0]}")\n')
        out.append(f"(2 {dim})\n")
        out.append(f"(10 (0 1 {hx(npoints)} 0 {dim}))\n")
        out.append(f"(13 (0 1 {hx(len(faces))} 0))\n")
        out.append(f"(12 (0 1 {hx(len(cell_to_global))} 0))\n")
        out.append(
            f"({'3010' if binary else '10'} ({hx(node_zone)} 1 {hx(npoints)} 1 {dim})"
        )
        fh.write("".join(out).encode())
        if binary:
            fh.write(b"\n(")
            fh.write(np.ascontiguousarray(points[:, :dim], dtype="<f8").tobytes())
            fh.write(b")\nEnd of Binary Section 3010)\n")
        else:
            fh.write(b"(\n")
            fh.write(
                "".join(
                    " ".join(f"{v:.16e}" for v in p) + "\n" for p in points[:, :dim]
                ).encode()
            )
            fh.write(b"))\n")

        def write_ints(key, head, rows):
            if binary:
                fh.write(f"(20{key} ({head})\n(".encode())
                flat = [v for row in rows for v in row]
                if flat and max(flat) > np.iinfo(np.int32).max:
                    raise WriteError(
                        "Fluent: an id does not fit a binary 32-bit section"
                    )
                fh.write(np.asarray(flat, dtype="<i4").tobytes())
                fh.write(f")\nEnd of Binary Section 20{key})\n".encode())
                return
            text = [f"({key} ({head})(\n"]
            text.extend(" ".join(hx(v) for v in row) + "\n" for row in rows)
            text.append("))\n")
            fh.write("".join(text).encode())

        first = 1
        for zid, _, _, _, members in cell_zones:
            types = [_element_type(mesh.cells[cell_block[c]].type) for c in members]
            mixed = any(t != types[0] for t in types)
            last = first + len(members) - 1
            head = f"{hx(zid)} {hx(first)} {hx(last)} 1 {hx(0 if mixed else types[0])}"
            if mixed:
                write_ints("12", head, [[t] for t in types])
            else:
                fh.write(f"(12 ({head}))\n".encode())
            first = last + 1

        first = 1
        for zid, ztype, _, _, members in face_zones:
            if not members:
                continue
            sizes = [len(faces[f][0]) for f in members]
            uniform = all(s == sizes[0] for s in sizes)
            if uniform and 2 <= sizes[0] <= 4:
                ftype = sizes[0]
            else:
                ftype = 5 if max(sizes) > 4 else 0
            last = first + len(members) - 1
            bc = "2" if ztype == "interior" else "3"
            head = f"{hx(zid)} {hx(first)} {hx(last)} {bc} {hx(ftype)}"
            rows = []
            for f in members:
                nodes, c0, c1 = faces[f]
                row = [len(nodes)] if ftype in (0, 5) else []
                row.extend(v + 1 for v in nodes)
                row.append(fluent_cell[c0])
                row.append(0 if c1 < 0 else fluent_cell[c1])
                rows.append(row)
            write_ints("13", head, rows)
            first = last + 1

        names = []
        for zone in cell_zones:
            names.append(f"(45 ({zone[0]} fluid {zone_name(zone)})())\n")
        for zone in face_zones:
            if zone[4]:
                names.append(f"(45 ({zone[0]} {zone[1]} {zone_name(zone)})())\n")
        fh.write("".join(names).encode())
