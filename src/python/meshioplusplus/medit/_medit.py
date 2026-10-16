"""
I/O for Medit's format/Gamma Mesh Format,
Latest official up-to-date documentation and a reference C implementation at
<https://github.com/LoicMarechal/libMeshb>
"""

import struct
from ctypes import c_double, c_float

import numpy as np

from .._common import _pick_first_int_data, warn
from .._exceptions import ReadError, WriteError
from .._files import open_file
from .._mesh import Mesh
from .._regions import Region
from ._medit_internal import medit_codes


def read(filename):
    if str(filename)[-1] == "b":
        with open_file(filename, "rb") as f:
            return read_binary_buffer(f)
    with open_file(filename) as f:
        return read_ascii_buffer(f)


def _produce_dtype(string_type, dim, itype, ftype):
    """
    convert a medit_code to a dtype appropriate for building a numpy array
    """
    res = ""
    c = 0
    while c < len(string_type):
        s = string_type[c]
        if s == "i":
            res += itype
        elif s == "r":
            res += ftype
        elif s == "d":
            res += str(dim)
            c += 1
            continue
        else:
            raise ReadError("Invalid string type")
        c += 1
        if c != len(string_type):
            res += ","
    return res


# libMeshb codes of the sections meshio++ reads: element sections and id lists.
_BINARY_CELLS = {
    5: ("line", 2),
    6: ("triangle", 3),
    7: ("quad", 4),
    8: ("tetra", 4),
    9: ("wedge", 6),
    10: ("hexahedron", 8),
    49: ("pyramid", 5),
}
_BINARY_ID_LISTS = {
    13: "Corners",
    14: "Ridges",
    15: "RequiredVertices",
    16: "RequiredEdges",
    17: "RequiredTriangles",
}


class _Cursor:
    """Bounds-checked reads from the bytes of a binary file."""

    def __init__(self, data, order):
        self.data = data
        self.at = 0
        self.order = order

    def left(self):
        return len(self.data) - self.at

    def array(self, count, dtype):
        dtype = np.dtype(self.order + dtype)
        if count < 0 or count * dtype.itemsize > self.left():
            raise ReadError("Medit: the file is truncated")
        out = np.frombuffer(self.data, dtype, count, self.at)
        self.at += count * dtype.itemsize
        return out

    def int(self, dtype):
        return int(self.array(1, dtype)[0])


def read_binary_buffer(f):
    data = f.read()
    cursor = _Cursor(data, "")
    code = cursor.int("i4")
    if code == 16777216:
        # swap endianness
        cursor.order = ">" if struct.unpack("=l", struct.pack("<l", 1))[0] == 1 else "<"
    elif code != 1:
        raise ReadError("Invalid code")

    version = cursor.int("i4")
    if version < 1 or version > 4:
        raise ReadError("Invalid version")
    itype = "i8" if version == 4 else "i4"
    ftype = "f4" if version == 1 else "f8"
    postype = "i8" if version >= 3 else "i4"

    field = cursor.int("i4")
    if field != 3:  # =  GmfDimension
        raise ReadError("Invalid dimension code : " + str(field) + " it should be 3")
    cursor.int(postype)
    dim = cursor.int("i4")
    if dim != 2 and dim != 3:
        raise ReadError("Invalid mesh dimension : " + str(dim))

    points = None
    cells = []
    point_data = {}
    cell_data = {"medit:ref": []}
    id_lists = {}
    while True:
        if cursor.left() < 4:
            warn("End-of-file reached before GmfEnd keyword")
            break
        field = cursor.int("i4")
        if field not in medit_codes:
            raise ReadError("Unsupported field")
        name = medit_codes[field][0]
        if name == "GmfEnd":
            break
        if name == "GmfReserved":
            continue
        following = cursor.int(postype)
        section_start = cursor.at

        if field == 4:  # GmfVertices
            n = cursor.int(itype)
            row = dim * np.dtype(ftype).itemsize + np.dtype(itype).itemsize
            if n < 0 or n * row > cursor.left():
                raise ReadError("Medit: the file is truncated")
            raw = np.frombuffer(
                data,
                np.dtype(
                    [("x", cursor.order + ftype, (dim,)), ("ref", cursor.order + itype)]
                ),
                n,
                cursor.at,
            )
            cursor.at += n * row
            points = raw["x"].astype(ftype)
            point_data["medit:ref"] = raw["ref"].astype(itype)
        elif field in _BINARY_CELLS:
            meshio_type, ncols = _BINARY_CELLS[field]
            n = cursor.int(itype)
            flat = cursor.array(n * (ncols + 1), itype).reshape(n, ncols + 1)
            cells.append((meshio_type, flat[:, :ncols].astype(itype) - 1))
            cell_data["medit:ref"].append(flat[:, -1].astype(itype))
        elif field in _BINARY_ID_LISTS:
            n = cursor.int(itype)
            id_lists[_BINARY_ID_LISTS[field]] = cursor.array(n, itype).astype(itype) - 1
        else:
            # A section meshio++ has no home for: skipped through the position
            # of the next one.
            if following <= section_start or following > len(data):
                raise ReadError(
                    f"Medit: cannot skip keyword {field}: its next-section position "
                    "is outside the file"
                )
            warn(f"meshio++ doesn't know Medit keyword {field}. Skipping.")
            cursor.at = following

    if points is None:
        raise ReadError("Medit: expected Vertices")
    regions = _regions_from_id_lists(id_lists, len(points), cells)
    return Mesh(
        points, cells, point_data=point_data, cell_data=cell_data, regions=regions
    )


# Keyword -> (region name, kind, the cell type its ids index). Corners and
# RequiredVertices list vertices; Ridges and RequiredEdges list entries of the
# Edges section, RequiredTriangles of the Triangles section.
_ID_LIST_KEYWORDS = {
    "Corners": ("medit:Corners", "point", None),
    "RequiredVertices": ("medit:RequiredVertices", "point", None),
    "Ridges": ("medit:Ridges", "cell", "line"),
    "RequiredEdges": ("medit:RequiredEdges", "cell", "line"),
    "RequiredTriangles": ("medit:RequiredTriangles", "cell", "triangle"),
}


def _regions_from_id_lists(id_lists, num_points, cells):
    """The regions the id-list keywords describe, ids checked against their section."""
    regions = []
    for keyword, ids in id_lists.items():
        name, kind, cell_type = _ID_LIST_KEYWORDS[keyword]
        if kind == "point":
            limit, base = num_points, 0
        else:
            base = 0
            limit = None
            for cell_type_i, data in cells:
                if cell_type_i == cell_type:
                    limit = len(data)
                    break
                base += len(data)
            if limit is None:
                warn(f"Medit: {keyword} lists {cell_type} cells the file does not have")
                continue
        if len(ids) and (ids.min() < 0 or ids.max() >= limit):
            raise ReadError(f"Medit: {keyword} names an entry the file does not have")
        regions.append(Region(name, kind, ids + base))
    return regions


def read_ascii_buffer(f):
    dim = 0
    cells = []
    point_data = {}
    cell_data = {"medit:ref": []}

    meshio_from_medit = {
        "Edges": ("line", 2),
        "Triangles": ("triangle", 3),
        "Quadrilaterals": ("quad", 4),
        "Tetrahedra": ("tetra", 4),
        "Prisms": ("wedge", 6),
        "Pyramids": ("pyramid", 5),
        "Hexahedra": ("hexahedron", 8),  # Frey
        "Hexaedra": ("hexahedron", 8),  # Dobrzynski
    }
    points = None
    dtype = None
    id_lists = {}  # keyword -> 0-based ids, kept as regions

    while True:
        line = f.readline()
        if not line:
            # EOF
            break

        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue

        items = line.split()

        if not items[0].isalpha():
            raise ReadError()

        if items[0] == "MeshVersionFormatted":
            version = items[1]
            dtype = {"0": c_float, "1": c_float, "2": c_double}[version]
        elif items[0] == "Dimension":
            if len(items) >= 2:
                dim = int(items[1])
            else:
                dim = int(
                    int(f.readline())
                )  # e.g. Dimension\n3, where the number of dimensions is on the next line
        elif items[0] == "Vertices":
            if dtype is None:
                raise ReadError("Expected `MeshVersionFormatted` before `Vertices`")
            num_verts = int(f.readline())
            first = []
            if dim <= 0:
                # No `Dimension` keyword (FEconv writes none): a vertex row
                # holds the coordinates and a reference.
                dim = 3
                if num_verts > 0:
                    line = ""
                    while not line.strip():
                        line = f.readline()
                        if not line:  # EOF: readline() returns "" for ever
                            raise ReadError("Medit: file ends inside `Vertices`")
                    first = [float(v) for v in line.split("#")[0].split()]
                    dim = len(first) - 1
                if dim < 1 or dim > 3:
                    raise ReadError(
                        "Medit: cannot tell the dimension from the Vertices rows"
                    )
            rest = np.fromfile(
                f,
                count=(num_verts - (1 if first else 0)) * (dim + 1),
                dtype=dtype,
                sep=" ",
            )
            out = np.concatenate([np.array(first, dtype=dtype), rest]).reshape(
                num_verts, dim + 1
            )
            points = out[:, :dim]
            point_data["medit:ref"] = out[:, dim].astype(int)
        elif items[0] in meshio_from_medit:
            meshio_type, points_per_cell = meshio_from_medit[items[0]]
            # The first value is the number of elements
            num_cells = int(f.readline())

            out = np.fromfile(
                f, count=num_cells * (points_per_cell + 1), dtype=int, sep=" "
            ).reshape(num_cells, points_per_cell + 1)

            # adapt for 0-base
            cells.append((meshio_type, out[:, :points_per_cell] - 1))
            cell_data["medit:ref"].append(out[:, -1])
        elif items[0] in _ID_LIST_KEYWORDS:
            # Vertex or element ids (1-based), kept as a region.
            count = int(f.readline())
            ids = np.fromfile(f, count=count, dtype=int, sep=" ")
            id_lists[items[0]] = ids - 1
        elif items[0] == "Normals":
            # those are just discarded
            num_normals = int(f.readline())
            np.fromfile(f, count=num_normals * dim, dtype=dtype, sep=" ").reshape(
                num_normals, dim
            )
        elif items[0] == "NormalAtVertices":
            # those are just discarded
            num_normal_at_vertices = int(f.readline())
            np.fromfile(
                f, count=num_normal_at_vertices * 2, dtype=int, sep=" "
            ).reshape(num_normal_at_vertices, 2)
        elif items[0] == "SubDomainFromMesh":
            # those are just discarded
            num_sub_domain_from_mesh = int(f.readline())
            np.fromfile(
                f, count=num_sub_domain_from_mesh * 4, dtype=int, sep=" "
            ).reshape(num_sub_domain_from_mesh, 4)
        elif items[0] == "VertexOnGeometricVertex":
            # those are just discarded
            num_vertex_on_geometric_vertex = int(f.readline())
            np.fromfile(
                f, count=num_vertex_on_geometric_vertex * 2, dtype=int, sep=" "
            ).reshape(num_vertex_on_geometric_vertex, 2)
        elif items[0] == "VertexOnGeometricEdge":
            # those are just discarded
            num_vertex_on_geometric_edge = int(f.readline())
            np.fromfile(
                f, count=num_vertex_on_geometric_edge * 3, dtype=float, sep=" "
            ).reshape(num_vertex_on_geometric_edge, 3)
        elif items[0] == "EdgeOnGeometricEdge":
            # those are just discarded
            num_edge_on_geometric_edge = int(f.readline())
            np.fromfile(
                f, count=num_edge_on_geometric_edge * 2, dtype=int, sep=" "
            ).reshape(num_edge_on_geometric_edge, 2)
        elif items[0] == "Identifier" or items[0] == "Geometry":
            f.readline()
        elif items[0] in [
            "TangentAtVertices",
            "Tangents",
        ]:
            msg = f"Meshio doesn't know keyword {items[0]}. Skipping."
            warn(msg)
            num_to_pass = int(f.readline())
            for _ in range(num_to_pass):
                f.readline()
        else:
            if items[0] != "End":
                raise ReadError(f"Unknown keyword '{items[0]}'.")

    if points is None:
        raise ReadError("Expected `Vertices`")
    regions = _regions_from_id_lists(id_lists, len(points), cells)
    return Mesh(
        points, cells, point_data=point_data, cell_data=cell_data, regions=regions
    )


def write(filename, mesh, float_fmt=".16e"):
    if str(filename)[-1] == "b":
        write_binary_file(filename, mesh)
    else:
        write_ascii_file(filename, mesh, float_fmt)


def write_ascii_file(filename, mesh, float_fmt=".16e"):
    with open_file(filename, "wb") as fh:
        version = {np.dtype(c_float): 1, np.dtype(c_double): 2}[mesh.points.dtype]
        # N. B.: PEP 461 Adding % formatting to bytes and bytearray
        fh.write(f"MeshVersionFormatted {version}\n".encode())

        n, d = mesh.points.shape

        fh.write(f"Dimension {d}\n".encode())

        # vertices
        fh.write(b"\nVertices\n")
        fh.write(f"{n}\n".encode())

        # pick out point data
        labels_key, other = _pick_first_int_data(mesh.point_data)
        if labels_key and other:
            string = ", ".join(other)
            warn(
                "Medit can only write one point data array. "
                f"Picking {labels_key}, skipping {string}."
            )
        labels = mesh.point_data[labels_key] if labels_key else np.ones(n, dtype=int)

        fmt = " ".join(["{:" + float_fmt + "}"] * d) + " {:d}\n"
        for x, label in zip(mesh.points, labels):
            fh.write(fmt.format(*x, label).encode())

        medit_from_meshio = {
            "line": ("Edges", 2),
            "triangle": ("Triangles", 3),
            "quad": ("Quadrilaterals", 4),
            "tetra": ("Tetrahedra", 4),
            "wedge": ("Prisms", 6),
            "pyramid": ("Pyramids", 5),
            "hexahedron": ("Hexahedra", 8),
        }

        # pick out cell_data
        labels_key, other = _pick_first_int_data(mesh.cell_data)
        if labels_key and other:
            string = ", ".join(other)
            warn(
                "Medit can only write one cell data array. "
                f"Picking {labels_key}, skipping {string}."
            )

        for k, cell_block in enumerate(mesh.cells):
            cell_type = cell_block.type
            data = cell_block.data
            try:
                medit_name, num = medit_from_meshio[cell_type]
            except KeyError:
                msg = f"MEDIT's mesh format doesn't know {cell_type} cells. Skipping."
                warn(msg)
                continue
            fh.write(b"\n")
            fh.write(f"{medit_name}\n".encode())
            fh.write(f"{len(data)}\n".encode())

            # pick out cell data
            labels = (
                mesh.cell_data[labels_key][k]
                if labels_key
                else np.ones(len(data), dtype=data.dtype)
            )

            fmt = " ".join(["{:d}"] * (num + 1)) + "\n"
            # adapt 1-base
            for d, label in zip(data + 1, labels):
                fh.write(fmt.format(*d, label).encode())

        _write_id_lists(fh, mesh)
        fh.write(b"\nEnd\n")


def _id_list_sections(mesh):
    """The regions named like the id-list keywords: ``(keyword, 0-based ids)``,
    cell regions' ids counted within the first block of the section's type."""
    sizes = [len(cb.data) for cb in mesh.cells]
    starts = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
    out = []
    for keyword, (name, kind, cell_type) in _ID_LIST_KEYWORDS.items():
        for region in getattr(mesh, "regions", []):
            if region.name != name or region.kind != kind:
                continue
            ids = np.asarray(region.entries, dtype=np.int64)
            if kind == "cell":
                block = next(
                    (i for i, cb in enumerate(mesh.cells) if cb.type == cell_type),
                    None,
                )
                if block is None:
                    warn(f"Medit: region '{name}' needs {cell_type} cells; not written")
                    continue
                lo, hi = starts[block], starts[block + 1]
                kept = ids[(ids >= lo) & (ids < hi)]
                if len(kept) != len(ids):
                    warn(
                        f"Medit: region '{name}' holds cells outside the first "
                        f"{cell_type} block; they are not written"
                    )
                ids = kept - lo
            out.append((keyword, ids))
            break
    return out


def _write_id_lists(fh, mesh):
    """The id-list regions as their keywords."""
    for keyword, ids in _id_list_sections(mesh):
        fh.write(f"\n{keyword}\n{len(ids)}\n".encode())
        fh.write("".join(f"{int(i) + 1}\n" for i in ids).encode())


def _first_int_key(data):
    """The first integer array's name in sorted order (the core's order), else ``None``."""
    for key in sorted(data):
        if data[key][0].dtype.kind in ("i", "u"):
            return key
    return None


# element type -> libMeshb code; the id lists' keyword -> code
_BINARY_CELL_CODES = {
    "line": 5,
    "triangle": 6,
    "quad": 7,
    "tetra": 8,
    "wedge": 9,
    "pyramid": 49,
    "hexahedron": 10,
}
_BINARY_ID_CODES = {
    "Corners": 13,
    "Ridges": 14,
    "RequiredVertices": 15,
    "RequiredEdges": 16,
    "RequiredTriangles": 17,
}


def write_binary_file(f, mesh):
    num_verts, dim = mesh.points.shape
    if dim not in (2, 3):
        raise WriteError(f"Medit: a binary file holds 2-D or 3-D points, not {dim}-D")

    point_key = _first_int_key(mesh.point_data)
    if point_key and len(mesh.point_data) > 1:
        warn(
            "Medit can only write one point data array. Picking the first integer one."
        )
    labels = (
        np.asarray(mesh.point_data[point_key], dtype=np.int64)
        if point_key
        else np.ones(num_verts, dtype=np.int64)
    )
    cell_key = _first_int_key(mesh.cell_data)
    if cell_key and len(mesh.cell_data) > 1:
        warn("Medit can only write one cell data array. Picking the first integer one.")

    blocks = []
    for k, cell_block in enumerate(mesh.cells):
        code = _BINARY_CELL_CODES.get(cell_block.type)
        if code is None:
            warn(
                f"MEDIT's mesh format doesn't know {cell_block.type} cells. "
                "Skipping."
            )
            continue
        data = np.asarray(cell_block.data, dtype=np.int64)
        refs = (
            np.asarray(mesh.cell_data[cell_key][k], dtype=np.int64)
            if cell_key
            else np.ones(len(data), dtype=np.int64)
        )
        blocks.append((code, data, refs))
    sections = _id_list_sections(mesh)

    # version 3 (32-bit integers) unless an index or reference needs 64 bits
    limit = np.iinfo(np.int32)
    values = [labels] + [d + 1 for _, d, _ in blocks] + [r for _, _, r in blocks]
    values += [ids + 1 for _, ids in sections]
    big = num_verts > limit.max or any(
        len(v) and (v.min() < limit.min or v.max() > limit.max) for v in values
    )
    version = 4 if big else 3
    isize = 8 if big else 4
    itype = "<i8" if big else "<i4"
    psize = 8
    pos = 4 * 4 + psize

    out = [np.array([1, version, 3], "<i4").tobytes(), np.array([pos], "<i8").tobytes()]
    out.append(np.array([dim], "<i4").tobytes())

    def header(code, count, nbytes):
        nonlocal pos
        pos += nbytes + 4 + psize + isize
        out.append(np.array([code], "<i4").tobytes())
        out.append(np.array([pos], "<i8").tobytes())
        out.append(np.array([count], itype).tobytes())

    header(4, num_verts, num_verts * (dim * 8 + isize))
    vertices = np.empty(num_verts, dtype=[("x", "<f8", (dim,)), ("ref", itype)])
    vertices["x"] = mesh.points
    vertices["ref"] = labels
    out.append(vertices.tobytes())
    for code, data, refs in blocks:
        num_cells, nodes = data.shape
        header(code, num_cells, num_cells * (nodes + 1) * isize)
        rows = np.empty((num_cells, nodes + 1), dtype=itype)
        rows[:, :nodes] = data + 1
        rows[:, nodes] = refs
        out.append(rows.tobytes())
    for keyword, ids in sections:
        header(_BINARY_ID_CODES[keyword], len(ids), len(ids) * isize)
        out.append((ids + 1).astype(itype).tobytes())
    out.append(np.array([54], "<i4").tobytes())
    out.append(np.array([0], "<i8").tobytes())
    with open_file(f, "wb") as fh:
        fh.write(b"".join(out))


# NOTE: format registration now lives in meshioplusplus/medit/__init__.py, which wraps the
# reader/writer above with the C++-backed fast paths (ascii .mesh).
