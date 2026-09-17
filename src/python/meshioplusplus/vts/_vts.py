"""VTK XML StructuredGrid, the pure-Python reference.

StructuredGrid states the same `nx * ny * nz` hexahedron topology ImageData
does -- points ordered x-fastest, the same `meshioplusplus._grid` numbering
:func:`meshioplusplus.grid` and :func:`meshioplusplus.voxelize` produce -- but
with an explicit ``<Points>`` array instead of ``Origin``/``Spacing``
attributes. That is the whole difference: connectivity is still implicit (the
index formula, never written), so this reader never requires the file's
points to sit on an even grid at all -- unlike :mod:`meshioplusplus.vti`, a
genuinely curved structured mesh reads correctly here.

The writer, unlike the reader, does require a dense lattice (via
:func:`meshioplusplus._grid.lattice_from_mesh`, exactly as ``vti``'s does):
the writer has no way to recover which ``(nx, ny, nz)`` a mesh's points were
meant to tile without one. Once confirmed, the mesh's own points are written
unchanged -- there is nothing to recompute, unlike ``vti``'s Origin/Spacing.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

import numpy as np

from .. import _provenance
from .._exceptions import ReadError, WriteError
from .._grid import _lattice_py, lattice_from_mesh
from .._mesh import Mesh
from ..vti._vti import _COMPRESSION_TO_ATTR, _ArrayReader, _encode_binary, _parse_n
from ..vtu._vtu import numpy_to_vtu_type


def _hex_conn(dims):
    """Hexahedron connectivity for a `dims`-cell grid -- the CSR the writer's
    own points already imply. Discards `_lattice_py`'s own points (which it
    would compute from a dummy origin/spacing): only `conn` depends on
    `dims` alone, so calling it this way is not wasted work, just a reuse of
    the one function that owns this index formula."""
    _, conn = _lattice_py(dims, np.zeros(3), np.ones(3))
    return conn


def read(filename):
    tree = ET.parse(str(filename))
    root = tree.getroot()
    if root.tag != "VTKFile":
        raise ReadError("Expected tag 'VTKFile'")
    if root.get("type") != "StructuredGrid":
        raise ReadError("Expected type StructuredGrid")

    compression = root.get("compressor")
    if compression == "vtkLZMADataCompressor":
        raise ReadError("lzma-compressed VTS is not supported")
    header_type = root.get("header_type", "UInt32")
    byte_order = root.get("byte_order")

    appended = root.find("AppendedData")
    appended_data = None
    if appended is not None:
        encoding = appended.get("encoding", "base64")
        if encoding != "base64":
            raise ReadError(f"VTS appended data encoding '{encoding}' is not supported")
        text = appended.text or ""
        appended_data = text.strip().lstrip("_")

    grid = root.find("StructuredGrid")
    if grid is None:
        raise ReadError("No StructuredGrid found")
    whole = _parse_n(grid.get("WholeExtent"), 6, int)
    if whole is None:
        raise ReadError("StructuredGrid has no readable WholeExtent")

    pieces = grid.findall("Piece")
    if not pieces:
        raise ReadError("No Piece found")
    if len(pieces) > 1:
        raise ReadError("multi-piece VTS is not supported")
    piece = pieces[0]
    piece_extent = _parse_n(piece.get("Extent"), 6, int)
    if piece_extent is not None and not np.array_equal(piece_extent, whole):
        raise ReadError(
            "VTS Piece Extent differs from WholeExtent; a partial piece "
            "is not supported"
        )

    dims = np.array([whole[2 * k + 1] - whole[2 * k] for k in range(3)], dtype=np.int64)
    if np.any(dims < 0):
        raise ReadError("VTS WholeExtent is inverted")
    num_points = int(np.prod(dims + 1))
    num_cells = int(np.prod(dims)) if np.all(dims > 0) else 0

    reader = _ArrayReader(header_type, byte_order, compression, appended_data)
    points_node = piece.find("Points")
    if points_node is None:
        raise ReadError("VTS Piece has no Points")
    points_da = points_node.find("DataArray")
    if points_da is None:
        raise ReadError("VTS Piece has no Points/DataArray")
    points = reader.read_data(points_da)
    if points.shape[0] != num_points:
        raise ReadError(
            f"VTS Points has {points.shape[0]} rows, but WholeExtent has "
            f"{num_points} points"
        )

    cells = [] if num_cells == 0 else [("hexahedron", _hex_conn(dims))]

    point_data = {}
    cell_data = {}
    for section, sink, expected, what in (
        ("PointData", point_data, num_points, "point"),
        ("CellData", cell_data, num_cells, "cell"),
    ):
        node = piece.find(section)
        if node is None:
            continue
        for da in node.findall("DataArray"):
            name = da.get("Name")
            arr = reader.read_data(da)
            if arr.size and arr.shape[0] != expected:
                raise ReadError(
                    f"VTS {what} array '{name}' has {arr.shape[0]} rows, but the "
                    f"extent has {expected} {what}s"
                )
            if what == "cell":
                if num_cells == 0:
                    continue
                sink[name] = [arr]
            else:
                sink[name] = arr

    return Mesh(points, cells, point_data=point_data, cell_data=cell_data)


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    spec = lattice_from_mesh(mesh)
    if spec is None:
        raise WriteError(
            "StructuredGrid is a lattice, and this mesh is not one: it needs "
            "exactly one hexahedron block whose points tile an axis-aligned box "
            "with uniform spacing (the writer recovers WholeExtent from it, then "
            "writes the mesh's own points unchanged). A partial grid "
            "(voxelize's 'surface'/'inside' fill, or an octree) cannot be "
            "written as .vts either -- write it as .vtu, which stores the "
            "cells explicitly."
        )
    dims, _origin, _spacing = spec
    if header_type is None:
        header_type = "UInt32"
    if compression not in (None, "zlib", "lzma", "lz4", "zstd"):
        raise WriteError(f"Unknown VTS compression '{compression}'")
    if not binary:
        compression = None

    ext = f"0 {int(dims[0])} 0 {int(dims[1])} 0 {int(dims[2])}"
    lines = ['<?xml version="1.0"?>']
    attrs = 'type="StructuredGrid" version="0.1" byte_order="LittleEndian"'
    if compression:
        attrs += f' compressor="{_COMPRESSION_TO_ATTR[compression]}"'
    if header_type != "UInt32":
        attrs += f' header_type="{header_type}"'
    lines.append(f"<VTKFile {attrs}>")
    lines.append(_provenance.render_xml_comment(_provenance.SlotTier.BLOCK))
    lines.append(f'<StructuredGrid WholeExtent="{ext}">')
    lines.append(f'<Piece Extent="{ext}">')

    def emit(section, items):
        if not items:
            return
        lines.append(f"<{section}>")
        for name, data in items:
            data = np.asarray(data)
            data = data.astype(data.dtype.newbyteorder("="), copy=False)
            vtu_type = numpy_to_vtu_type[data.dtype]
            head = f'<DataArray type="{vtu_type}" Name="{name}"'
            if data.ndim == 2:
                head += f' NumberOfComponents="{data.shape[1]}"'
            if binary:
                lines.append(head + ' format="binary">')
                lines.append(_encode_binary(data, compression, header_type))
            else:
                lines.append(head + ' format="ascii">')
                fmt = "{:.11e}" if vtu_type.startswith("Float") else "{:d}"
                lines.extend(fmt.format(v) for v in data.reshape(-1))
            lines.append("</DataArray>")
        lines.append(f"</{section}>")

    points = np.asarray(mesh.points, dtype=np.float64)
    emit("Points", [("Points", points)])
    emit("PointData", sorted(mesh.point_data.items()))
    emit(
        "CellData",
        [
            (k, np.concatenate([np.asarray(a) for a in v]))
            for k, v in sorted(mesh.cell_data.items())
        ],
    )

    lines.append("</Piece>")
    lines.append("</StructuredGrid>")
    lines.append("</VTKFile>")
    with open(filename, "w") as f:
        f.write("\n".join(lines) + "\n")
