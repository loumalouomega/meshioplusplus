"""
I/O for EnSight Gold (.case/.geo), c.f.
<https://vis.lbl.gov/archive/NERSC/Software/ensight/doc/OnlineHelp/UM-C11.pdf>

Handles the FORMAT/GEOMETRY/TIME/VARIABLE sections of the .case file and the
Gold geometry file in ASCII, C-binary and Fortran-binary form (32-bit ints/floats, 80-char
string records, framed as Fortran sequential records in Fortran binary;
a foreign byte order is auto-detected from the plausibility of the
part-number/node-count records). Per the Gold specification connectivity is
positional (1-based index into the part's coordinate list), so "node id
given/ignore" id arrays are skipped. Multi-part files are concatenated into
one point array with the owning part recorded as the ``ensight:part``
cell_data field (only when there are two or more parts). The writer emits a
single part with ``node id assign`` / ``element id assign``, one variable file
per point/cell array and a ``constant per case`` line per scalar field_data
entry. ``VARIABLE`` entries are read back (one file per step, the ``TIME``
section selecting the step) as ``point_data``/``cell_data``/``field_data``. ``polygon`` blocks are written as
``nsided`` and ``polyhedron<N>`` blocks as ``nfaced``, the exact inverse of
what the reader parses (counts, then for ``nfaced`` the face sizes, then the
node ids).
"""

import pathlib

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._fortran_records import fortran_records, sniff_fortran_records
from .._mesh import CellBlock, Mesh

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


_CASE_SECTIONS = (
    "FORMAT",
    "GEOMETRY",
    "VARIABLE",
    "TIME",
    "FILE",
    "MATERIAL",
    "SCRIPTS",
)
_VARIABLE_KINDS = (
    "scalar per node:",
    "vector per node:",
    "tensor symm per node:",
    "tensor asym per node:",
    "scalar per element:",
    "vector per element:",
    "tensor symm per element:",
    "tensor asym per element:",
)


def _tokens_after_leading_ints(value):
    """The tokens of a record without its leading ``[ts] [fs]`` integers."""
    tokens = value.split()
    while tokens:
        try:
            int(tokens[0])
        except ValueError:
            break
        tokens.pop(0)
    return tokens


def _leading_int(text):
    """``strtol``-style: the integer a string starts with, else 0."""
    digits = ""
    text = text.strip()
    for i, ch in enumerate(text):
        if ch.isdigit() or (i == 0 and ch in "+-"):
            digits += ch
        else:
            break
    try:
        return int(digits)
    except ValueError:
        return 0


def _parse_case_info(case_path):
    """The .case file's GEOMETRY, TIME and VARIABLE sections."""
    section = None
    format_type = ""
    model_value = ""
    info = {
        "time_values": [],
        "start": 0,
        "increment": 1,
        "variables": [],
        "constants": [],
    }
    in_time_values = False
    have_time_set = False
    num_steps = -1
    with open(case_path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                in_time_values = False
                continue
            if line in _CASE_SECTIONS:
                section = line
                in_time_values = False
                continue
            if section == "FORMAT" and line.startswith("type:"):
                format_type = line[5:].strip()
            elif section == "GEOMETRY" and line.startswith("model:"):
                model_value = line[6:].strip()
            elif section == "VARIABLE" and line.startswith("constant per case:"):
                toks = _tokens_after_leading_ints(line[len("constant per case:") :])
                if len(toks) < 2:
                    raise ReadError(
                        f"EnSight: malformed 'constant per case' line: {line}"
                    )
                try:
                    value = float(toks[-1])
                except ValueError:
                    raise ReadError(
                        f"EnSight: malformed 'constant per case' line: {line}"
                    )
                info["constants"].append((" ".join(toks[:-1]), value))
            elif section == "VARIABLE":
                for kind in _VARIABLE_KINDS:
                    if not line.startswith(kind):
                        continue
                    toks = _tokens_after_leading_ints(line[len(kind) :])
                    if len(toks) < 2:
                        raise ReadError(
                            f"EnSight: malformed '{kind[:-1]}' line: {line}"
                        )
                    info["variables"].append(
                        {
                            "kind": kind[:-1],
                            "name": " ".join(toks[:-1]),
                            "pattern": toks[-1],
                        }
                    )
                    break
            elif section == "TIME":
                if line.startswith("time set:"):
                    if have_time_set:
                        section = None  # only the first time set is honoured
                        continue
                    have_time_set = True
                elif line.startswith("number of steps:"):
                    num_steps = _leading_int(line[len("number of steps:") :])
                elif line.startswith("filename start number:"):
                    info["start"] = _leading_int(line[len("filename start number:") :])
                elif line.startswith("filename increment:"):
                    info["increment"] = _leading_int(line[len("filename increment:") :])
                elif line.startswith("time values:"):
                    in_time_values = True
                    info["time_values"] += _floats_prefix(line[len("time values:") :])
                elif in_time_values:
                    info["time_values"] += _floats_prefix(line)
                    if num_steps >= 0 and len(info["time_values"]) >= num_steps:
                        in_time_values = False

    if "ensight gold" not in format_type:
        raise ReadError(
            f"EnSight: case file is not 'type: ensight gold' (got '{format_type}')"
        )
    if not model_value:
        raise ReadError("EnSight: case file has no GEOMETRY 'model:' entry")

    # model: [ts] [fs] filename [change_coords_only] — drop leading integer
    # timeset/fileset tokens, take the first remaining token as the filename.
    tokens = _tokens_after_leading_ints(model_value)
    if not tokens:
        raise ReadError("EnSight: malformed 'model:' line in case file")
    geo_name = tokens[0]
    if "*" in geo_name:
        raise ReadError("EnSight: transient (wildcard) geometry is not supported")
    info["geo_path"] = pathlib.Path(case_path).parent / geo_name
    return info


def _floats_prefix(text):
    """The leading run of numbers of a line (``>>`` into a double stops at a non-number)."""
    out = []
    for tok in text.split():
        try:
            out.append(float(tok))
        except ValueError:
            break
    return out


def _parse_case(case_path):
    """Parse the .case file and return the resolved geometry file path."""
    return _parse_case_info(case_path)["geo_path"]


def _resolve_wildcard(pattern, number):
    """Replace the run of ``*`` with ``number``, zero-padded to the run's width."""
    star = pattern.find("*")
    if star < 0:
        return pattern
    width = 0
    while star + width < len(pattern) and pattern[star + width] == "*":
        width += 1
    digits = f"{number:0{width}d}"
    if len(digits) > width:
        digits = digits[-width:]
    return pattern[:star] + digits + pattern[star + width :]


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
    """Record/number stream over a C-binary Gold geometry file (or a
    Fortran-binary one with its record markers stripped)."""

    def __init__(self, data, order="=", start=80):
        # `order` seeds the byte order when the file says it (Fortran markers).
        self.data = data
        self.pos = start  # past the leading "C Binary"/"Fortran Binary" record
        self.i4 = np.dtype(np.int32).newbyteorder(order)
        self.f4 = np.dtype(np.float32).newbyteorder(order)

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


def _parse_geo(cur, layout=None):
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
    parts = []  # (part_id, point_offset, num_points), in file order
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
        parts.append((part_id, point_offset, nn))

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
    if layout is not None:
        for part_id, offset, count in parts:
            layout.append(
                {
                    "id": part_id,
                    "offset": offset,
                    "count": count,
                    "blocks": [
                        (i, ne)
                        for i, b in enumerate(blocks)
                        if b[3] == part_id
                        for ne in (b[2],)
                    ],
                }
            )
    cell_data = {}
    if num_parts >= 2:
        cell_data["ensight:part"] = [
            np.full(ne, part_id, dtype=np.int64) for _, _, ne, part_id in blocks
        ]
    return Mesh(points, cells, cell_data=cell_data)


def _unframe_fortran(data):
    """The joined payloads of a Fortran-binary geometry file and their byte order
    (``"<"`` or ``">"``), or ``None``.

    Fortran binary is the C-binary stream with every WRITE framed as a Fortran
    sequential unformatted record; its first record is the 80-byte ``Fortran
    Binary`` string. Joining the payloads gives the C-binary stream back.
    """
    layout = sniff_fortran_records(data)
    if layout is None:
        return None
    width = layout[0]
    if not data[width : width + 80].startswith(b"Fortran Binary"):
        return None
    records = fortran_records(data, layout, "EnSight geometry file")
    return b"".join(data[off : off + size] for off, size in records), layout[1]


def _unframe_variable_fortran(data):
    """Like :func:`_unframe_fortran` for a variable file, which has no format record."""
    layout = sniff_fortran_records(data)
    if layout is None:
        return None
    records = fortran_records(data, layout, "EnSight variable file")
    return b"".join(data[off : off + size] for off, size in records), layout[1]


def _variable_ncomp(kind):
    if kind.startswith("vector"):
        return 3
    if kind.startswith("tensor symm"):
        return 6
    if kind.startswith("tensor asym"):
        return 9
    return 1


def _read_variable_file(cur, per_node, ncomp, layout, total_points):
    """One variable file against the geometry's part/block layout."""
    cur.next_record()  # description line
    point_out = None
    if per_node:
        shape = (total_points,) if ncomp == 1 else (total_points, ncomp)
        point_out = np.zeros(shape, dtype=np.float64)
    block_out = {}
    for part in layout:
        rec = cur.next_record()
        if not rec.startswith("part"):
            raise ReadError(
                f"EnSight: expected 'part' record in variable file, got: {rec}"
            )
        cur.check_swap(_PLAUSIBLE_MAX, prefer_smaller=True)
        if cur.next_int() != part["id"]:
            raise ReadError(
                "EnSight: variable file's part sequence does not match the geometry's"
            )
        if per_node:
            rec = cur.next_record()
            if not rec.startswith("coordinates"):
                raise ReadError(
                    "EnSight: expected 'coordinates' record in variable file, "
                    f"got: {rec}"
                )
            lo, n = part["offset"], part["count"]
            for c in range(ncomp):
                values = cur.read_floats(n)
                if ncomp == 1:
                    point_out[lo : lo + n] = values
                else:
                    point_out[lo : lo + n, c] = values
            continue
        for block_index, num_cells in part["blocks"]:
            cur.next_record()  # element type; not re-validated by name
            shape = (num_cells,) if ncomp == 1 else (num_cells, ncomp)
            block = np.zeros(shape, dtype=np.float64)
            for c in range(ncomp):
                values = cur.read_floats(num_cells)
                if ncomp == 1:
                    block[:] = values
                else:
                    block[:, c] = values
            block_out[block_index] = block
    return point_out if per_node else block_out


def _read_variable_file_auto(path, encoding, per_node, ncomp, layout, total_points):
    data = pathlib.Path(path).read_bytes()
    binary_layout = (
        len(data) >= 160 and data[80:].startswith(b"part") and b"\n" not in data[:80]
    ) or data.startswith(b"C Binary")
    if encoding != "ascii" and not binary_layout:
        head = data[:1024]
        eol = head.find(b"\n")
        if b"\0" not in head and 0 <= eol < 81:
            encoding = "ascii"
    if encoding == "fortran":
        unframed = _unframe_variable_fortran(data)
        if unframed is None:
            raise ReadError(
                f"EnSight: variable file '{path}' is not Fortran binary like its "
                "geometry file"
            )
        cur = _BinaryCursor(unframed[0], unframed[1], start=0)
    elif encoding == "cbinary":
        standard = len(data) >= 160 and data[80:].startswith(b"part")
        start = 80 if (not standard and data.startswith(b"C Binary")) else 0
        cur = _BinaryCursor(data, start=start)
    else:
        cur = _AsciiCursor(data.decode("utf-8", errors="replace"))
    return _read_variable_file(cur, per_node, ncomp, layout, total_points)


def _swap_tensor_symm(arr):
    """EnSight's ``11 22 33 12 13 23`` <-> meshio++'s ``xx yy zz xy yz zx``."""
    out = np.array(arr, copy=True)
    if out.size:
        out[:, [4, 5]] = out[:, [5, 4]]
    return out


def read(filename, time_step=0):
    filename = pathlib.Path(filename)
    case_info = None
    if filename.suffix == ".case":
        case_info = _parse_case_info(filename)
        geo_path = case_info["geo_path"]
    else:
        geo_path = filename
    with open(geo_path, "rb") as f:
        data = f.read()
    layout = []
    unframed = _unframe_fortran(data)
    if unframed is not None:
        encoding = "fortran"
        mesh = _parse_geo(_BinaryCursor(*unframed), layout)
    elif len(data) >= 80 and data.startswith(b"C Binary"):
        encoding = "cbinary"
        mesh = _parse_geo(_BinaryCursor(data), layout)
    else:
        encoding = "ascii"
        mesh = _parse_geo(_AsciiCursor(data.decode("utf-8", errors="replace")), layout)

    if case_info is None or not (case_info["variables"] or case_info["constants"]):
        return mesh

    for name, value in case_info["constants"]:
        mesh.field_data[name] = np.array([value], dtype=np.float64)
    if not case_info["variables"]:
        return mesh

    times = case_info["time_values"]
    step = 0
    if times:
        step = time_step
        if step < 0:
            step += len(times)
        if not 0 <= step < len(times):
            raise ReadError(
                f"EnSight: time step {time_step} out of range for {len(times)} step(s)"
            )
    elif time_step != 0:
        raise ReadError(
            "EnSight: mTimeStep requested but the case file has no TIME section to "
            "resolve it against"
        )
    file_number = case_info["start"] + step * case_info["increment"]

    for var in case_info["variables"]:
        per_node = "per node" in var["kind"]
        per_element = "per element" in var["kind"]
        if not per_node and not per_element:
            continue
        ncomp = _variable_ncomp(var["kind"])
        tensor_symm = var["kind"].startswith("tensor symm")
        resolved = _resolve_wildcard(var["pattern"], file_number)
        var_path = filename.parent / resolved
        values = _read_variable_file_auto(
            var_path, encoding, per_node, ncomp, layout, len(mesh.points)
        )
        if per_node:
            mesh.point_data[var["name"]] = (
                _swap_tensor_symm(values) if tensor_symm else values
            )
        else:
            blocks = []
            for i in range(len(mesh.cells)):
                blk = values.get(i)
                if blk is None:
                    blk = np.zeros(
                        (len(mesh.cells[i]),)
                        if ncomp == 1
                        else (len(mesh.cells[i]), ncomp)
                    )
                blocks.append(
                    _swap_tensor_symm(blk) if tensor_symm and blk.size else blk
                )
            mesh.cell_data[var["name"]] = blocks
    return mesh


def _str80(s):
    return s.encode("ascii")[:79].ljust(80, b"\0")


def _ragged_layout(cell_block):
    """The ``nsided`` / ``nfaced`` runs of a ``polygon`` / ``polyhedron`` block.

    Returns ``(keyword, counts, sizes, rings)``: the per-cell node count
    (``nsided``) or face count (``nfaced``), the per-face node count (``nfaced``
    only, else empty), and the 0-based node-id arrays in file order (one per
    polygon, or one per face).
    """
    if cell_block.type.startswith("polyhedron"):
        rings = [
            np.asarray(face, dtype=np.int64)
            for cell in cell_block.data
            for face in cell
        ]
        counts = [len(cell) for cell in cell_block.data]
        return "nfaced", counts, [ring.size for ring in rings], rings
    rings = [np.asarray(row, dtype=np.int64) for row in cell_block.data]
    return "nsided", [ring.size for ring in rings], [], rings


def _write_geo_ascii(fh, points, cells):
    fh.write(b"EnSight Gold Geometry File\n")
    fh.write(f"{_provenance.lines(_provenance.SlotTier.BOUNDED)[0]}\n".encode())
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
        if cell_block.type.startswith(("polygon", "polyhedron")):
            keyword, counts, sizes, rings = _ragged_layout(cell_block)
            fh.write(f"{keyword}\n".encode())
            fh.write(f"{len(counts):10d}\n".encode())
            fh.write("".join(f"{n:10d}\n" for n in counts + sizes).encode())
            fh.write(
                "".join(
                    "".join(f"{v + 1:10d}" for v in ring) + "\n" for ring in rings
                ).encode()
            )
            continue
        conn = _permute(cell_block.type, np.asarray(cell_block.data)) + 1
        fh.write(f"{meshio_to_ensight_type[cell_block.type]}\n".encode())
        fh.write(f"{conn.shape[0]:10d}\n".encode())
        fh.write(
            "".join("".join(f"{v:10d}" for v in row) + "\n" for row in conn).encode()
        )


def _write_geo_binary(fh, points, cells, fortran=False):
    # One record per string, count or array: Fortran binary frames each as a
    # Fortran sequential unformatted record (native 4-byte length markers), the
    # layout vtkEnSightGoldBinaryReader reads.
    def record(payload):
        if fortran:
            if len(payload) > np.iinfo(np.int32).max:
                raise WriteError(
                    "EnSight: a Fortran-binary record over 2 GiB cannot be written"
                )
            marker = np.int32(len(payload)).tobytes()
            fh.write(marker + payload + marker)
        else:
            fh.write(payload)

    int32_max = np.iinfo(np.int32).max
    if points.shape[0] > int32_max:
        raise WriteError("EnSight: mesh too large for 32-bit binary EnSight output")
    record(_str80("Fortran Binary" if fortran else "C Binary"))
    record(_str80("EnSight Gold Geometry File"))
    record(_str80(_provenance.lines(_provenance.SlotTier.BOUNDED)[0]))
    record(_str80("node id assign"))
    record(_str80("element id assign"))
    record(_str80("part"))
    record(np.int32(1).tobytes())
    record(_str80("Mesh"))
    record(_str80("coordinates"))
    record(np.int32(points.shape[0]).tobytes())
    for c in range(3):
        record(np.ascontiguousarray(points[:, c], dtype=np.float32).tobytes())
    for cell_block in cells:
        if cell_block.type.startswith(("polygon", "polyhedron")):
            keyword, counts, sizes, rings = _ragged_layout(cell_block)
            if len(counts) > int32_max:
                raise WriteError(
                    "EnSight: mesh too large for 32-bit binary EnSight output"
                )
            record(_str80(keyword))
            record(np.int32(len(counts)).tobytes())
            nodes = np.concatenate(rings) + 1 if rings else np.empty(0, dtype=np.int64)
            # Each non-empty run is one record, as the C++ writer emits them.
            for run in (counts, sizes, nodes):
                if len(run):
                    record(np.asarray(run, dtype=np.int32).tobytes())
            continue
        conn = _permute(cell_block.type, np.asarray(cell_block.data)) + 1
        if conn.shape[0] > int32_max:
            raise WriteError("EnSight: mesh too large for 32-bit binary EnSight output")
        record(_str80(meshio_to_ensight_type[cell_block.type]))
        record(np.int32(conn.shape[0]).tobytes())
        record(np.ascontiguousarray(conn, dtype=np.int32).tobytes())


def _block_keyword(cell_block):
    if cell_block.type.startswith("polygon"):
        return "nsided"
    if cell_block.type.startswith("polyhedron"):
        return "nfaced"
    return meshio_to_ensight_type[cell_block.type]


def _kind_for_ncomp(ncomp):
    """``(kind word, written component count)`` or ``None`` for an unsupported count."""
    return {
        1: ("scalar", 1),
        2: ("vector", 3),
        3: ("vector", 3),
        6: ("tensor symm", 6),
        9: ("tensor asym", 9),
    }.get(ncomp)


def _variable_extension(kind, per_node):
    ext = {
        "scalar": "scl",
        "vector": "vec",
        "tensor symm": "tsym",
        "tensor asym": "tasym",
    }[kind]
    return ext if per_node else "e" + ext


def _ncomp(arr):
    arr = np.asarray(arr)
    return arr.shape[1] if arr.ndim >= 2 else 1


def _collect_variables(mesh):
    """The arrays this writer can express: (variables, constants)."""
    variables = []
    constants = []
    for name in sorted(mesh.point_data):
        ncomp = _ncomp(mesh.point_data[name])
        kind = _kind_for_ncomp(ncomp)
        if kind is None:
            warn(
                f"EnSight: skipping point data '{name}' with {ncomp} components "
                "(1, 2, 3, 6 or 9 are supported)"
            )
            continue
        variables.append((name, kind[0], kind[1], True))
    for name in sorted(mesh.cell_data):
        blocks = mesh.cell_data[name]
        if len(blocks) != len(mesh.cells):
            warn(
                f"EnSight: skipping cell data '{name}': not present on every cell block"
            )
            continue
        ncomp = _ncomp(blocks[0]) if len(blocks) else 1
        kind = _kind_for_ncomp(ncomp)
        if kind is None:
            warn(
                f"EnSight: skipping cell data '{name}' with {ncomp} components "
                "(1, 2, 3, 6 or 9 are supported)"
            )
            continue
        variables.append((name, kind[0], kind[1], False))
    for name in sorted(mesh.field_data):
        arr = np.asarray(mesh.field_data[name])
        if arr.size != 1:
            warn(
                f"EnSight: skipping field data '{name}': only a single scalar can "
                "become a 'constant per case'"
            )
            continue
        constants.append((name, float(arr.reshape(-1)[0])))
    return variables, constants


def _variable_columns(array, kind, written):
    """The ``written`` components of ``array`` (n,) or (n, k), in EnSight file order."""
    array = np.asarray(array, dtype=np.float64)
    stored = array.shape[1] if array.ndim >= 2 else 1
    columns = []
    for comp in range(written):
        source = comp
        if kind == "tensor symm" and comp in (4, 5):
            source = 9 - comp
        if source < stored:
            columns.append(array if array.ndim < 2 else array[:, source])
        else:
            columns.append(np.zeros(array.shape[0]))
    return columns


def _write_variable(path, mesh, variable, binary, fortran):
    name, kind, written, per_node = variable
    if per_node:
        sections = [
            ("coordinates", _variable_columns(mesh.point_data[name], kind, written))
        ]
    else:
        sections = [
            (_block_keyword(block), _variable_columns(data, kind, written))
            for block, data in zip(mesh.cells, mesh.cell_data[name])
        ]
    with open(path, "wb") as fh:
        if binary:

            def record(payload):
                if fortran:
                    marker = np.int32(len(payload)).tobytes()
                    fh.write(marker + payload + marker)
                else:
                    fh.write(payload)

            record(_str80("variable"))
            record(_str80("part"))
            record(np.int32(1).tobytes())
            for keyword, columns in sections:
                record(_str80(keyword))
                for column in columns:
                    record(np.ascontiguousarray(column, dtype=np.float32).tobytes())
        else:
            fh.write(b"variable\npart\n" + f"{1:10d}\n".encode())
            for keyword, columns in sections:
                fh.write(f"{keyword}\n".encode())
                for column in columns:
                    fh.write("".join(f"{v:12.5e}\n" for v in column).encode())


def write(filename, mesh, binary=True, fortran=False):
    if fortran and not binary:
        raise WriteError("EnSight: Fortran binary is a binary encoding (binary=True)")
    case_path, geo_path = _case_geo_paths(filename)
    if case_path is None:
        raise WriteError(f"EnSight: must specify a .case or .geo file. Got {filename}.")

    if mesh.points.shape[1] > 3:
        raise WriteError("EnSight: points must have at most three components")
    for cell_block in mesh.cells:
        if cell_block.type.startswith(("polygon", "polyhedron")):
            continue  # nsided / nfaced
        if cell_block.type not in meshio_to_ensight_type:
            raise WriteError(
                f"EnSight: cell type '{cell_block.type}' has no EnSight keyword"
            )
    variables, constants = _collect_variables(mesh)

    points = np.asarray(mesh.points, dtype=np.float64)
    if points.shape[1] < 3:
        points = np.column_stack(
            [points, np.zeros((points.shape[0], 3 - points.shape[1]))]
        )

    text = (
        "FORMAT\n" "type: ensight gold\n" "\n" "GEOMETRY\n" f"model: {geo_path.name}\n"
    )
    if variables or constants:
        text += "\nVARIABLE\n"
        for name, kind, _, per_node in variables:
            ext = _variable_extension(kind, per_node)
            text += (
                f"{kind} per {'node' if per_node else 'element'}: "
                f"{name} {name}.{ext}\n"
            )
        for name, value in constants:
            text += f"constant per case: {name} {format(value, '.17g')}\n"
    with open(case_path, "wb") as fh:
        fh.write(text.encode())

    with open(geo_path, "wb") as fh:
        if binary:
            _write_geo_binary(fh, points, mesh.cells, fortran)
        else:
            _write_geo_ascii(fh, points, mesh.cells)

    for variable in variables:
        name, kind, _, per_node = variable
        var_path = geo_path.parent / f"{name}.{_variable_extension(kind, per_node)}"
        _write_variable(var_path, mesh, variable, binary, fortran)
