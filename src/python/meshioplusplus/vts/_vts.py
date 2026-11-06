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

import numpy as np

from .. import _provenance
from .._exceptions import ReadError, WriteError
from .._grid import _lattice_py, lattice_from_mesh
from .._mesh import Mesh
from .._region_field_data import field_arrays_for_write
from .._vtk_xml_read import field_arrays
from .._vtk_xml_read import load as _load_xml
from .._vtk_xml_read import read_pieces
from ..vti._vti import (
    _COMPRESSION_TO_ATTR,
    emit_data_arrays,
    extent_dims,
    read_extent_data,
)
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
    root, reader = _load_xml(filename, "StructuredGrid", "VTS")
    return read_pieces(root, reader, "StructuredGrid", "VTS", _read_piece)


def _read_piece(grid, piece, reader):
    dims, num_points, num_cells = extent_dims(grid, piece, "VTS")
    if num_points == 0:
        return Mesh(np.empty((0, 3)), [], field_data=field_arrays(piece, reader, "VTS"))

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

    point_data, cell_data = read_extent_data(
        piece, reader, num_points, num_cells, "VTS"
    )

    return Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_arrays(piece, reader, "VTS"),
    )


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

    def emit(section, items, field=False):
        emit_data_arrays(
            lines, section, items, binary, compression, header_type, field=field
        )

    emit(
        "FieldData", field_arrays_for_write(mesh, "VTS", numpy_to_vtu_type), field=True
    )
    lines.append(f'<Piece Extent="{ext}">')
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
