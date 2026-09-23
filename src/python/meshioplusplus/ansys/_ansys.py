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

A cell section with a connectivity body is meshio's own (legacy) layout, the
one :func:`write` produces; such a file reads as before, cells only.
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


def read(filename):  # noqa: C901
    with open_file(filename, "rb") as f:
        data = f.read()
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


def write(filename, mesh, binary=True):
    with open_file(filename, "wb") as fh:
        # header
        fh.write(
            f'(1 "{_provenance.lines(_provenance.SlotTier.SINGLE_LINE)[0]}")\n'.encode()
        )

        # dimension
        num_points, dim = mesh.points.shape
        if dim not in [2, 3]:
            raise WriteError(f"Can only write dimension 2, 3, got {dim}.")
        fh.write((f"(2 {dim})\n").encode())

        # total number of nodes
        first_node_index = 1
        fh.write((f"(10 (0 {first_node_index:x} {num_points:x} 0))\n").encode())

        # total number of cells
        total_num_cells = sum(len(c) for c in mesh.cells)
        fh.write((f"(12 (0 1 {total_num_cells:x} 0))\n").encode())

        # Write nodes
        key = "3010" if binary else "10"
        fh.write(
            f"({key} (1 {first_node_index:x} {num_points:x} 1 {dim:x})(\n".encode()
        )
        if binary:
            mesh.points.tofile(fh)
            fh.write(b"\n)")
            fh.write(b"End of Binary Section 3010)\n")
        else:
            np.savetxt(fh, mesh.points, fmt="%.16e")
            fh.write(b"))\n")

        # Write cells
        meshio_to_ansys_type = {
            # "mixed": 0,
            "triangle": 1,
            "tetra": 2,
            "quad": 3,
            "hexahedron": 4,
            "pyramid": 5,
            "wedge": 6,
            # "polyhedral": 7,
        }
        first_index = 0
        binary_dtypes = {
            # np.int16 is not allowed
            np.dtype("int32"): "2012",
            np.dtype("int64"): "3012",
        }
        for cell_block in mesh.cells:
            cell_type = cell_block.type
            values = cell_block.data
            key = binary_dtypes[values.dtype] if binary else "12"
            last_index = first_index + len(values) - 1
            try:
                ansys_cell_type = meshio_to_ansys_type[cell_type]
            except KeyError:
                legal_keys = ", ".join(meshio_to_ansys_type.keys())
                raise KeyError(
                    f"Illegal ANSYS cell type '{cell_type}'. (legal: {legal_keys})"
                )
            fh.write(
                f"({key} (1 {first_index:x} {last_index:x} 1 {ansys_cell_type})(\n".encode()
            )
            if binary:
                (values + first_node_index).tofile(fh)
                fh.write(b"\n)")
                fh.write((f"End of Binary Section {key})\n").encode())
            else:
                np.savetxt(fh, values + first_node_index, fmt="%x")
                fh.write(b"))\n")
            first_index = last_index + 1
