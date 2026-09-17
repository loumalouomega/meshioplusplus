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

import xml.etree.ElementTree as ET

import numpy as np

from .. import _provenance
from .._exceptions import ReadError, WriteError
from .._grid import lattice_from_mesh
from .._mesh import Mesh
from ..vti._vti import _COMPRESSION_TO_ATTR, _ArrayReader, _encode_binary, _parse_n
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
        [base, base + 1, base + px + 1, base + px, top, top + 1, top + px + 1, top + px],
        axis=1,
    ).astype(np.int64)


def _read_axis(coordinates, name, expected, reader):
    for da in coordinates.findall("DataArray"):
        if da.get("Name") != name:
            continue
        arr = reader.read_data(da)
        if arr.shape[0] != expected:
            raise ReadError(
                f"VTR {name} has {arr.shape[0]} entries, but WholeExtent needs {expected}"
            )
        return arr.astype(np.float64)
    raise ReadError(f"VTR Coordinates has no '{name}' DataArray")


def read(filename):
    tree = ET.parse(str(filename))
    root = tree.getroot()
    if root.tag != "VTKFile":
        raise ReadError("Expected tag 'VTKFile'")
    if root.get("type") != "RectilinearGrid":
        raise ReadError("Expected type RectilinearGrid")

    compression = root.get("compressor")
    if compression == "vtkLZMADataCompressor":
        raise ReadError("lzma-compressed VTR is not supported")
    header_type = root.get("header_type", "UInt32")
    byte_order = root.get("byte_order")

    appended = root.find("AppendedData")
    appended_data = None
    if appended is not None:
        encoding = appended.get("encoding", "base64")
        if encoding != "base64":
            raise ReadError(f"VTR appended data encoding '{encoding}' is not supported")
        text = appended.text or ""
        appended_data = text.strip().lstrip("_")

    grid = root.find("RectilinearGrid")
    if grid is None:
        raise ReadError("No RectilinearGrid found")
    whole = _parse_n(grid.get("WholeExtent"), 6, int)
    if whole is None:
        raise ReadError("RectilinearGrid has no readable WholeExtent")

    pieces = grid.findall("Piece")
    if not pieces:
        raise ReadError("No Piece found")
    if len(pieces) > 1:
        raise ReadError("multi-piece VTR is not supported")
    piece = pieces[0]
    piece_extent = _parse_n(piece.get("Extent"), 6, int)
    if piece_extent is not None and not np.array_equal(piece_extent, whole):
        raise ReadError(
            "VTR Piece Extent differs from WholeExtent; a partial piece "
            "is not supported"
        )

    dims = np.array([whole[2 * k + 1] - whole[2 * k] for k in range(3)], dtype=np.int64)
    if np.any(dims < 0):
        raise ReadError("VTR WholeExtent is inverted")
    num_points = int(np.prod(dims + 1))
    num_cells = int(np.prod(dims)) if np.all(dims > 0) else 0

    reader = _ArrayReader(header_type, byte_order, compression, appended_data)
    coords = piece.find("Coordinates")
    if coords is None:
        raise ReadError("VTR Piece has no Coordinates")
    xs = _read_axis(coords, "x_coordinates", dims[0] + 1, reader)
    ys = _read_axis(coords, "y_coordinates", dims[1] + 1, reader)
    zs = _read_axis(coords, "z_coordinates", dims[2] + 1, reader)

    gz, gy, gx = np.meshgrid(zs, ys, xs, indexing="ij")
    points = np.stack([gx.reshape(-1), gy.reshape(-1), gz.reshape(-1)], axis=1)
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
                    f"VTR {what} array '{name}' has {arr.shape[0]} rows, but the "
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
