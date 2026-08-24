"""
I/O for EnSight Gold (.case/.geo), geometry only, c.f.
<https://vis.lbl.gov/archive/NERSC/Software/ensight/doc/OnlineHelp/UM-C11.pdf>

Handles the FORMAT/GEOMETRY sections of the .case file and the Gold geometry
file in ASCII and C-binary form (32-bit ints/floats, 80-char string records;
a foreign byte order is auto-detected from the plausibility of the
part-number/node-count records). Per the Gold specification connectivity is
positional (1-based index into the part's coordinate list), so "node id
given/ignore" id arrays are skipped. Multi-part files are concatenated into
one point array with the owning part recorded as the ``ensight:part``
cell_data field (only when there are two or more parts). The writer emits a
single part with ``node id assign`` / ``element id assign`` and drops
point/cell/field data (mesh-only scope).
"""

import pathlib

import numpy as np

from .._common import warn
from .._exceptions import ReadError, WriteError
from .._mesh import CellBlock, Mesh
from .._provenance import TAG as _PROVENANCE_TAG

ensight_to_meshio_type = {
    "point": "vertex",
    "bar2": "line",
    "bar3": "line3",
    "tria3": "triangle",
    "tria6": "triangle6",
    "quad4": "quad",
    "quad8": "quad8",
    "tetra4": "tetra",
    "tetra10": "tetra10",
    "pyramid5": "pyramid",
    "pyramid13": "pyramid13",
    "penta6": "wedge",
    "penta15": "wedge15",
    "hexa8": "hexahedron",
    "hexa20": "hexahedron20",
}
meshio_to_ensight_type = {v: k for k, v in ensight_to_meshio_type.items()}
ensight_type_num_nodes = {
    "point": 1,
    "bar2": 2,
    "bar3": 3,
    "tria3": 3,
    "tria6": 6,
    "quad4": 4,
    "quad8": 8,
    "tetra4": 4,
    "tetra10": 10,
    "pyramid5": 5,
    "pyramid13": 13,
    "penta6": 6,
    "penta15": 15,
    "hexa8": 8,
    "hexa20": 20,
}

# meshio (== VTK) <-> EnSight node order differs only in the prism triangle
# winding of penta15 vs wedge15 — the same involution VTK's EnSight readers
# apply, so one table serves both directions.
_wedge15_permutation = [0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13]

_PLAUSIBLE_MAX = 100_000_000


def _permute(cell_type, conn):
    if cell_type == "wedge15":
        return conn[:, _wedge15_permutation]
    return conn


def _case_geo_paths(filename):
    filename = pathlib.Path(filename)
    if filename.suffix == ".case":
        return filename, filename.parent / (filename.stem + ".geo")
    if filename.suffix == ".geo":
        return filename.parent / (filename.stem + ".case"), filename
    return None, None


def _parse_case(case_path):
    """Parse the .case file and return the resolved geometry file path."""
    section = None
    format_type = ""
    model_value = ""
    with open(case_path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line in (
                "FORMAT",
                "GEOMETRY",
                "VARIABLE",
                "TIME",
                "FILE",
                "MATERIAL",
                "SCRIPTS",
            ):
                section = line
                continue
            if section == "FORMAT" and line.startswith("type:"):
                format_type = line[5:].strip()
            elif section == "GEOMETRY" and line.startswith("model:"):
                model_value = line[6:].strip()

    if "ensight gold" not in format_type:
        raise ReadError(
            f"EnSight: case file is not 'type: ensight gold' (got '{format_type}')"
        )
    if not model_value:
        raise ReadError("EnSight: case file has no GEOMETRY 'model:' entry")

    # model: [ts] [fs] filename [change_coords_only] — drop leading integer
    # timeset/fileset tokens, take the first remaining token as the filename.
    tokens = model_value.split()
    while tokens:
        try:
            int(tokens[0])
        except ValueError:
            break
        tokens.pop(0)
    if not tokens:
        raise ReadError("EnSight: malformed 'model:' line in case file")
    geo_name = tokens[0]
    if "*" in geo_name:
        raise ReadError("EnSight: transient (wildcard) geometry is not supported")
    return pathlib.Path(case_path).parent / geo_name


class _AsciiCursor:
    """Record/number stream over an ASCII Gold geometry file."""

    def __init__(self, text):
        self.lines = text.splitlines()
        self.li = 0
        self.toks = []

    def at_end(self):
        if self.toks:
            return False
        for line in self.lines[self.li :]:
            if line.strip():
                return False
        return True

    def next_record(self):
        # A string record always starts on a fresh line.
        self.toks = []
        while self.li < len(self.lines):
            line = self.lines[self.li].strip()
            self.li += 1
            if line:
                return line
        raise ReadError("EnSight: unexpected end of geometry file")

    def peek_record(self):
        if self.at_end():
            return ""
        saved_li, saved_toks = self.li, self.toks
        rec = self.next_record()
        self.li, self.toks = saved_li, saved_toks
        return rec

    def _next_token(self):
        while not self.toks:
            if self.li >= len(self.lines):
                raise ReadError("EnSight: unexpected end of geometry file")
            self.toks = self.lines[self.li].split()
            self.li += 1
        return self.toks.pop(0)

    def next_int(self):
        try:
            return int(self._next_token())
        except ValueError:
            raise ReadError("EnSight: expected an integer in geometry file")

    def read_ints(self, n):
        return np.array([self.next_int() for _ in range(n)], dtype=np.int64)

    def read_floats(self, n):
        try:
            return np.array(
                [float(self._next_token()) for _ in range(n)], dtype=np.float64
            )
        except ValueError:
            raise ReadError("EnSight: expected a number in geometry file")

    def skip_ints(self, n):
        for _ in range(n):
            self._next_token()

    def check_swap(self, maximum, prefer_smaller=False):
        pass


class _BinaryCursor:
    """Record/number stream over a C-binary Gold geometry file."""

    def __init__(self, data):
        self.data = data
        self.pos = 80  # past the leading "C Binary" record
        self.i4 = np.dtype(np.int32)
        self.f4 = np.dtype(np.float32)

    def at_end(self):
        return self.pos >= len(self.data)

    def next_record(self):
        if self.pos + 80 > len(self.data):
            raise ReadError("EnSight: truncated binary geometry file")
        rec = self.data[self.pos : self.pos + 80]
        self.pos += 80
        return rec.split(b"\0", 1)[0].decode("ascii", errors="replace").strip()

    def peek_record(self):
        if self.pos + 80 > len(self.data):
            return ""
        saved = self.pos
        rec = self.next_record()
        self.pos = saved
        return rec

    def next_int(self):
        return int(self.read_ints(1)[0])

    def read_ints(self, n):
        if self.pos + 4 * n > len(self.data):
            raise ReadError("EnSight: truncated binary geometry file")
        out = np.frombuffer(self.data, dtype=self.i4, count=n, offset=self.pos)
        self.pos += 4 * n
        return out.astype(np.int64)

    def read_floats(self, n):
        if self.pos + 4 * n > len(self.data):
            raise ReadError("EnSight: truncated binary geometry file")
        out = np.frombuffer(self.data, dtype=self.f4, count=n, offset=self.pos)
        self.pos += 4 * n
        return out.astype(np.float64)

    def skip_ints(self, n):
        if self.pos + 4 * n > len(self.data):
            raise ReadError("EnSight: truncated binary geometry file")
        self.pos += 4 * n

    def check_swap(self, maximum, prefer_smaller=False):
        """Enable byte-swapping when the next int32 is implausible as-is but
        plausible with the opposite byte order. With ``prefer_smaller`` (used
        at the tiny part-number record, where e.g. bswap(1) = 16777216 is
        still "plausible"), the smaller of two plausible interpretations
        wins."""
        if self.pos + 4 > len(self.data):
            return
        v = int(np.frombuffer(self.data, dtype=self.i4, count=1, offset=self.pos)[0])
        swapped = self.i4.newbyteorder("S")
        s = int(np.frombuffer(self.data, dtype=swapped, count=1, offset=self.pos)[0])
        v_ok = 0 <= v <= maximum
        s_ok = 0 <= s <= maximum
        if (not v_ok and s_ok) or (prefer_smaller and v_ok and s_ok and s < v):
            self.i4 = swapped
            self.f4 = self.f4.newbyteorder("S")


def _ids_in_file(record, what):
    # "given" and "ignore" both put id arrays in the file; only the presence
    # matters — Gold connectivity is positional, so ids are always skipped.
    mode = record.split()[-1]
    if mode in ("given", "ignore"):
        return True
    if mode in ("off", "assign"):
        return False
    raise ReadError(f"EnSight: malformed '{what} id' record: {record}")


def _parse_geo(cur):
    cur.next_record()  # description line 1
    cur.next_record()  # description line 2
    node_id_rec = cur.next_record()
    if not node_id_rec.startswith("node id"):
        raise ReadError(f"EnSight: expected 'node id' record, got: {node_id_rec}")
    elem_id_rec = cur.next_record()
    if not elem_id_rec.startswith("element id"):
        raise ReadError(f"EnSight: expected 'element id' record, got: {elem_id_rec}")
    node_ids_in_file = _ids_in_file(node_id_rec, "node")
    elem_ids_in_file = _ids_in_file(elem_id_rec, "element")

    if cur.peek_record().startswith("extents"):
        cur.next_record()
        cur.read_floats(6)

    all_points = []
    blocks = []  # (cell_type, data, num_cells, part_id)
    num_parts = 0

    while not cur.at_end():
        rec = cur.next_record()
        if not rec.startswith("part"):
            raise ReadError(f"EnSight: expected 'part' record, got: {rec}")
        cur.check_swap(_PLAUSIBLE_MAX, prefer_smaller=True)
        part_id = cur.next_int()
        num_parts += 1
        cur.next_record()  # part description

        rec = cur.next_record()
        if not rec.startswith("coordinates"):
            raise ReadError(f"EnSight: expected 'coordinates' record, got: {rec}")
        cur.check_swap(_PLAUSIBLE_MAX)
        nn = cur.next_int()
        if nn < 0:
            raise ReadError("EnSight: negative node count")
        point_offset = sum(p.shape[0] for p in all_points)

        if node_ids_in_file:
            cur.skip_ints(nn)
        x = cur.read_floats(nn)
        y = cur.read_floats(nn)
        z = cur.read_floats(nn)
        all_points.append(np.column_stack([x, y, z]))

        def resolve(conn, nn=nn, point_offset=point_offset):
            # 1-based positional index within this part -> global 0-based.
            out = np.asarray(conn, dtype=np.int64) - 1
            if out.size > 0 and (out.min() < 0 or out.max() >= nn):
                raise ReadError("EnSight: connectivity index out of range")
            return out + point_offset

        while not cur.at_end():
            kw = cur.peek_record()
            if not kw or kw.startswith("part"):
                break
            cur.next_record()
            if kw.startswith("g_"):  # ghost cells carry the same data
                kw = kw[2:]

            cur.check_swap(_PLAUSIBLE_MAX)
            ne = cur.next_int()
            if ne < 0:
                raise ReadError("EnSight: negative element count")
            if elem_ids_in_file:
                cur.skip_ints(ne)

            if kw == "nsided":
                sizes = cur.read_ints(ne)
                flat = resolve(cur.read_ints(int(sizes.sum())))
                offsets = np.concatenate([[0], np.cumsum(sizes)])
                rows = [flat[offsets[c] : offsets[c + 1]] for c in range(ne)]
                blocks.append(("polygon", rows, ne, part_id))
            elif kw == "nfaced":
                nfaces = cur.read_ints(ne)
                fsizes = cur.read_ints(int(nfaces.sum()))
                flat = resolve(cur.read_ints(int(fsizes.sum())))
                cells = []
                face_at = 0
                node_at = 0
                for c in range(ne):
                    faces = []
                    for _ in range(int(nfaces[c])):
                        size = int(fsizes[face_at])
                        face_at += 1
                        faces.append(flat[node_at : node_at + size])
                        node_at += size
                    cells.append(faces)
                # Group by unique node count into "polyhedron<N>" blocks,
                # preserving first-seen order (the openfoam convention).
                groups = {}
                order = []
                for cell in cells:
                    n = np.unique(np.concatenate(cell)).size
                    if n not in groups:
                        groups[n] = []
                        order.append(n)
                    groups[n].append(cell)
                for n in order:
                    blocks.append(
                        (f"polyhedron{n}", groups[n], len(groups[n]), part_id)
                    )
            else:
                if kw not in ensight_to_meshio_type:
                    raise ReadError(f"EnSight: unsupported element keyword: {kw}")
                npc = ensight_type_num_nodes[kw]
                cell_type = ensight_to_meshio_type[kw]
                conn = cur.read_ints(ne * npc).reshape(ne, npc)
                conn = resolve(_permute(cell_type, conn))
                blocks.append((cell_type, conn, ne, part_id))

    points = (
        np.concatenate(all_points) if all_points else np.empty((0, 3), dtype=np.float64)
    )
    cells = [CellBlock(cell_type, data) for cell_type, data, _, _ in blocks]
    cell_data = {}
    if num_parts >= 2:
        cell_data["ensight:part"] = [
            np.full(ne, part_id, dtype=np.int64) for _, _, ne, part_id in blocks
        ]
    return Mesh(points, cells, cell_data=cell_data)


def read(filename):
    filename = pathlib.Path(filename)
    geo_path = _parse_case(filename) if filename.suffix == ".case" else filename
    with open(geo_path, "rb") as f:
        data = f.read()
    if data.startswith(b"Fortran Binary"):
        raise ReadError("EnSight: Fortran-binary geometry files are not supported")
    if len(data) >= 80 and data.startswith(b"C Binary"):
        return _parse_geo(_BinaryCursor(data))
    return _parse_geo(_AsciiCursor(data.decode("utf-8", errors="replace")))


def _str80(s):
    return s.encode("ascii")[:79].ljust(80, b"\0")


def _write_geo_ascii(fh, points, cells):
    fh.write(b"EnSight Gold Geometry File\n")
    fh.write(f"{_PROVENANCE_TAG}\n".encode())
    fh.write(b"node id assign\n")
    fh.write(b"element id assign\n")
    fh.write(b"part\n")
    fh.write(f"{1:10d}\n".encode())
    fh.write(b"Mesh\n")
    fh.write(b"coordinates\n")
    fh.write(f"{points.shape[0]:10d}\n".encode())
    for c in range(3):
        fh.write("".join(f"{v:12.5e}\n" for v in points[:, c]).encode())
    for cell_block in cells:
        conn = _permute(cell_block.type, np.asarray(cell_block.data)) + 1
        fh.write(f"{meshio_to_ensight_type[cell_block.type]}\n".encode())
        fh.write(f"{conn.shape[0]:10d}\n".encode())
        fh.write(
            "".join("".join(f"{v:10d}" for v in row) + "\n" for row in conn).encode()
        )


def _write_geo_binary(fh, points, cells):
    int32_max = np.iinfo(np.int32).max
    if points.shape[0] > int32_max:
        raise WriteError("EnSight: mesh too large for 32-bit binary EnSight output")
    fh.write(_str80("C Binary"))
    fh.write(_str80("EnSight Gold Geometry File"))
    fh.write(_str80(_PROVENANCE_TAG))
    fh.write(_str80("node id assign"))
    fh.write(_str80("element id assign"))
    fh.write(_str80("part"))
    fh.write(np.int32(1).tobytes())
    fh.write(_str80("Mesh"))
    fh.write(_str80("coordinates"))
    fh.write(np.int32(points.shape[0]).tobytes())
    fh.write(np.ascontiguousarray(points.T, dtype=np.float32).tobytes())
    for cell_block in cells:
        conn = _permute(cell_block.type, np.asarray(cell_block.data)) + 1
        if conn.shape[0] > int32_max:
            raise WriteError("EnSight: mesh too large for 32-bit binary EnSight output")
        fh.write(_str80(meshio_to_ensight_type[cell_block.type]))
        fh.write(np.int32(conn.shape[0]).tobytes())
        fh.write(np.ascontiguousarray(conn, dtype=np.int32).tobytes())


def write(filename, mesh, binary=True):
    case_path, geo_path = _case_geo_paths(filename)
    if case_path is None:
        raise WriteError(f"EnSight: must specify a .case or .geo file. Got {filename}.")

    if mesh.points.shape[1] > 3:
        raise WriteError("EnSight: points must have at most three components")
    for cell_block in mesh.cells:
        if cell_block.type.startswith(("polygon", "polyhedron")):
            raise WriteError(
                "EnSight: writing nsided/nfaced (ragged) blocks is not supported"
            )
        if cell_block.type not in meshio_to_ensight_type:
            raise WriteError(
                f"EnSight: cell type '{cell_block.type}' has no EnSight keyword"
            )
    if mesh.point_data or mesh.cell_data or mesh.field_data:
        warn("EnSight writer is geometry-only; point/cell/field data are dropped.")

    points = np.asarray(mesh.points, dtype=np.float64)
    if points.shape[1] < 3:
        points = np.column_stack(
            [points, np.zeros((points.shape[0], 3 - points.shape[1]))]
        )

    with open(case_path, "wb") as fh:
        fh.write(
            (
                "FORMAT\n"
                "type: ensight gold\n"
                "\n"
                "GEOMETRY\n"
                f"model: {geo_path.name}\n"
            ).encode()
        )

    with open(geo_path, "wb") as fh:
        if binary:
            _write_geo_binary(fh, points, mesh.cells)
        else:
            _write_geo_ascii(fh, points, mesh.cells)
