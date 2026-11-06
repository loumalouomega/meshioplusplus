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

from .. import _provenance
from .._exceptions import ReadError, WriteError
from .._grid import lattice_from_mesh
from .._mesh import Mesh
from .._region_field_data import field_arrays_for_write
from .._vtk_xml_read import field_arrays
from .._vtk_xml_read import load as _load_xml
from .._vtk_xml_read import read_pieces
from ..vti._vti import (
    _COMPRESSION_TO_ATTR,
    _encode_binary,
    emit_data_arrays,
    extent_dims,
    read_extent_data,
)
from ..vtu._vtu import numpy_to_vtu_type


def _hex_conn(dims):
    nx, ny, nz = (int(d) for d in dims)
    px, py = nx + 1, ny + 1
    ci, cj, ck = np.meshgrid(
        np.arange(nx, dtype=np.int64),
        np.arange(ny, dtype=np.int64),
        np.arange(nz, dtype=np.int64),
        indexing="ij",
    )
    ci = ci.transpose(2, 1, 0).reshape(-1)
    cj = cj.transpose(2, 1, 0).reshape(-1)
    ck = ck.transpose(2, 1, 0).reshape(-1)
    base = (ck * py + cj) * px + ci
    top = base + px * py
    return np.stack(
        [
            base,
            base + 1,
            base + px + 1,
            base + px,
            top,
            top + 1,
            top + px + 1,
            top + px,
        ],
        axis=1,
    ).astype(np.int64)


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
    dims, num_points, num_cells = extent_dims(grid, piece, "VTR")
    if num_points == 0:
        return Mesh(np.empty((0, 3)), [], field_data=field_arrays(piece, reader, "VTR"))

    coords = piece.find("Coordinates")
    if coords is None:
        raise ReadError("VTR Piece has no Coordinates")
    xs = _read_axis(coords, "x_coordinates", dims[0] + 1, reader)
    ys = _read_axis(coords, "y_coordinates", dims[1] + 1, reader)
    zs = _read_axis(coords, "z_coordinates", dims[2] + 1, reader)

    gz, gy, gx = np.meshgrid(zs, ys, xs, indexing="ij")
    points = np.stack([gx.reshape(-1), gy.reshape(-1), gz.reshape(-1)], axis=1)
    cells = [] if num_cells == 0 else [("hexahedron", _hex_conn(dims))]

    point_data, cell_data = read_extent_data(
        piece, reader, num_points, num_cells, "VTR"
    )

    return Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_arrays(piece, reader, "VTR"),
    )


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    spec = lattice_from_mesh(mesh)
    if spec is None:
        raise WriteError(
            "RectilinearGrid needs a UNIFORM dense lattice today: exactly one "
            "hexahedron block whose points tile an axis-aligned box with "
            "uniform per-axis spacing. A genuinely graded (non-uniform) mesh "
            "cannot be written as .vtr yet -- a documented follow-up, see "
            "doc/roadmap.md -- and a partial grid (voxelize's "
            "'surface'/'inside' fill, or an octree) cannot be written as "
            ".vtr either -- write it as .vtu, which stores the cells "
            "explicitly."
        )
    dims, origin, spacing = spec
    if header_type is None:
        header_type = "UInt32"
    if compression not in (None, "zlib", "lzma", "lz4", "zstd"):
        raise WriteError(f"Unknown VTR compression '{compression}'")
    if not binary:
        compression = None

    ext = f"0 {int(dims[0])} 0 {int(dims[1])} 0 {int(dims[2])}"
    lines = ['<?xml version="1.0"?>']
    attrs = 'type="RectilinearGrid" version="0.1" byte_order="LittleEndian"'
    if compression:
        attrs += f' compressor="{_COMPRESSION_TO_ATTR[compression]}"'
    if header_type != "UInt32":
        attrs += f' header_type="{header_type}"'
    lines.append(f"<VTKFile {attrs}>")
    lines.append(_provenance.render_xml_comment(_provenance.SlotTier.BLOCK))
    lines.append(f'<RectilinearGrid WholeExtent="{ext}">')

    def emit(section, items, field=False):
        emit_data_arrays(
            lines, section, items, binary, compression, header_type, field=field
        )

    emit(
        "FieldData", field_arrays_for_write(mesh, "VTR", numpy_to_vtu_type), field=True
    )
    lines.append(f'<Piece Extent="{ext}">')
    lines.append("<Coordinates>")
    for k, name in enumerate(("x_coordinates", "y_coordinates", "z_coordinates")):
        axis = origin[k] + np.arange(int(dims[k]) + 1, dtype=np.float64) * spacing[k]
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

    emit("PointData", sorted(mesh.point_data.items()))
    emit(
        "CellData",
        [
            (k, np.concatenate([np.asarray(a) for a in v]))
            for k, v in sorted(mesh.cell_data.items())
        ],
    )

    lines.append("</Piece>")
    lines.append("</RectilinearGrid>")
    lines.append("</VTKFile>")
    with open(filename, "w") as f:
        f.write("\n".join(lines) + "\n")
