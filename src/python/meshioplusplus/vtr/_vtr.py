"""VTK XML RectilinearGrid, the pure-Python reference.

RectilinearGrid states the same `nx * ny * nz` hexahedron topology
`.vti`/`.vts` do, but its `<Coordinates>` are three independent, only
*monotonic* 1-D arrays -- point `(i, j, k)` sits at
`(xs[i], ys[j], zs[k])`, the tensor product. That is genuinely more general
than `.vti`'s uniform spacing: a graded mesh (finer near a wall, coarser far
from it) is a RectilinearGrid, never an ImageData.

The reader is fully general -- it never checks that the coordinates are
evenly spaced. The writer, unlike the reader, requires a UNIFORM dense
lattice (via `meshioplusplus._grid.lattice_from_mesh`, exactly as `.vti`'s
does): recovering three arbitrary per-axis arrays from an unstructured point
set needs a detector this module does not implement. A genuinely graded
mesh cannot be written as `.vtr` today -- a documented follow-up, not a
silent gap.
"""

from __future__ import annotations

import numpy as np

from .._exceptions import ReadError
from .._vtk_xml_read import load as _load_xml
from .._vtk_xml_read import read_pieces
from ..vti._vti import _encode_binary, read_lattice_piece, write_lattice
from ..vtu._vtu import numpy_to_vtu_type


def _read_axis(coordinates, name, expected, reader):
    arrays = coordinates.findall("DataArray")
    da = next((a for a in arrays if a.get("Name") == name), None)
    if da is None:
        if len(arrays) != 3:
            raise ReadError("VTR Coordinates must have three axis DataArrays")
        da = arrays["xyz".index(name[0])]
    arr = reader.read_data(da).reshape(-1)
    if arr.size != expected:
        raise ReadError(
            f"VTR {name} has {arr.size} entries, but WholeExtent needs {expected}"
        )
    return arr.astype(np.float64)


def read(filename):
    root, reader = _load_xml(filename, "RectilinearGrid", "VTR")
    return read_pieces(root, reader, "RectilinearGrid", "VTR", _read_piece)


def _read_piece(grid, piece, reader):
    return read_lattice_piece(grid, piece, reader, "VTR", _read_points)


def _read_points(piece, reader, dims, num_points):
    coords = piece.find("Coordinates")
    if coords is None:
        raise ReadError("VTR Piece has no Coordinates")
    xs = _read_axis(coords, "x_coordinates", dims[0] + 1, reader)
    ys = _read_axis(coords, "y_coordinates", dims[1] + 1, reader)
    zs = _read_axis(coords, "z_coordinates", dims[2] + 1, reader)
    gz, gy, gx = np.meshgrid(zs, ys, xs, indexing="ij")
    return np.stack([gx.reshape(-1), gy.reshape(-1), gz.reshape(-1)], axis=1)


def _write_coordinates(binary, header_type):
    def geometry(lines, emit, dims, origin, spacing, compression):
        lines.append("<Coordinates>")
        for k, name in enumerate(("x_coordinates", "y_coordinates", "z_coordinates")):
            axis = (
                origin[k] + np.arange(int(dims[k]) + 1, dtype=np.float64) * spacing[k]
            )
            vtu_type = numpy_to_vtu_type[np.dtype(np.float64)]
            head = f'<DataArray type="{vtu_type}" Name="{name}"'
            if binary:
                lines.append(head + ' format="binary">')
                lines.append(_encode_binary(axis, compression, header_type))
            else:
                lines.append(head + ' format="ascii">')
                lines.extend(f"{v:.11e}" for v in axis)
            lines.append("</DataArray>")
        lines.append("</Coordinates>")

    return geometry


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    write_lattice(
        filename,
        mesh,
        grid="RectilinearGrid",
        prefix="VTR",
        not_lattice=(
            "RectilinearGrid needs a UNIFORM dense lattice today: exactly one "
            "hexahedron block whose points tile an axis-aligned box with "
            "uniform per-axis spacing. A genuinely graded (non-uniform) mesh "
            "cannot be written as .vtr yet -- a documented follow-up, see "
            "doc/roadmap.md -- and a partial grid (voxelize's "
            "'surface'/'inside' fill, or an octree) cannot be written as "
            ".vtr either -- write it as .vtu, which stores the cells "
            "explicitly."
        ),
        binary=binary,
        compression=compression,
        header_type="UInt32" if header_type is None else header_type,
        geometry=_write_coordinates(
            binary, "UInt32" if header_type is None else header_type
        ),
    )
