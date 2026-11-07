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

from .._exceptions import ReadError
from .._vtk_xml_read import load as _load_xml
from .._vtk_xml_read import read_pieces
from ..vti._vti import read_lattice_piece, write_lattice


def read(filename):
    root, reader = _load_xml(filename, "StructuredGrid", "VTS")
    return read_pieces(root, reader, "StructuredGrid", "VTS", _read_piece)


def _read_piece(grid, piece, reader):
    return read_lattice_piece(grid, piece, reader, "VTS", _read_points)


def _read_points(piece, reader, dims, num_points):
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
    return points


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    def geometry(lines, emit, dims, origin, spacing, compression):
        emit("Points", [("Points", np.asarray(mesh.points, dtype=np.float64))])

    write_lattice(
        filename,
        mesh,
        grid="StructuredGrid",
        prefix="VTS",
        not_lattice=(
            "StructuredGrid is a lattice, and this mesh is not one: it needs "
            "exactly one hexahedron block whose points tile an axis-aligned box "
            "with uniform spacing (the writer recovers WholeExtent from it, then "
            "writes the mesh's own points unchanged). A partial grid "
            "(voxelize's 'surface'/'inside' fill, or an octree) cannot be "
            "written as .vts either -- write it as .vtu, which stores the "
            "cells explicitly."
        ),
        binary=binary,
        compression=compression,
        header_type=header_type,
        geometry=geometry,
    )
