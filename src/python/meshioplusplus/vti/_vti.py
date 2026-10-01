"""VTK XML ImageData, the pure-Python reference.

ImageData is a regular lattice whose geometry is three attributes -- ``Origin``,
``Spacing`` and ``WholeExtent`` -- rather than a point array. A :class:`Mesh` has
no implicit geometry, so this module expands the extent into explicit points and
one ``hexahedron`` block on read (through :mod:`meshioplusplus._grid`, the same
numbering :func:`meshioplusplus.grid` and :func:`meshioplusplus.voxelize`
produce) and derives the three attributes back from the geometry on write.

That means the writer **requires a lattice**: a mesh that is not one has no
extent to write, and a *partial* one (``voxelize``'s ``surface``/``inside``
fills, or an octree) cannot be one either, since ImageData has no way to express
a hole. Both raise :class:`WriteError` by name.

The array encoding -- base64, the 32 KiB block framing and the codecs -- is the
same VTK XML machinery :mod:`meshioplusplus.vtu` uses, and the decode path here
literally reuses ``VtuReader``'s methods rather than transcribing them.
"""

from __future__ import annotations

import base64

import numpy as np

from .. import _provenance
from .._exceptions import ReadError, WriteError
from .._grid import _lattice_py, lattice_from_mesh
from .._mesh import Mesh
from .._vtk_xml_read import load as _load_xml
from .._vtk_xml_read import field_arrays, piece_extent, read_pieces
from ..vtu._vtu import (
    _COMPRESSION_TO_ATTR,
    _chunk_it,
    _compressor_for,
    numpy_to_vtu_type,
    vtu_to_numpy_type,
)


def _parse_n(text, count, dtype):
    if text is None:
        return None
    parts = text.replace(",", " ").split()
    if len(parts) < count:
        return None
    return np.array([dtype(x) for x in parts[:count]])


def read(filename):
    root, reader = _load_xml(filename, "ImageData", "VTI")
    return read_pieces(root, reader, "ImageData", "VTI", _read_piece)


def _read_piece(grid, piece, reader):
    whole = piece_extent(grid, piece)
    origin = _parse_n(grid.get("Origin"), 3, float)
    if origin is None:
        origin = np.zeros(3)
    spacing = _parse_n(grid.get("Spacing"), 3, float)
    if spacing is None:
        spacing = np.ones(3)
    direction = _parse_n(grid.get("Direction"), 9, float)
    if direction is not None and not np.array_equal(direction, np.eye(3).reshape(-1)):
        raise ReadError("VTI with a non-identity Direction is not supported")

    dims = np.array([whole[2 * k + 1] - whole[2 * k] for k in range(3)], dtype=np.int64)
    if np.any(dims < 0):
        raise ReadError("VTI WholeExtent is inverted")
    # A point at extent index i sits at Origin + i * Spacing, so an extent that
    # does not start at zero translates the mesh's own lo corner.
    lo = np.array([origin[k] + whole[2 * k] * spacing[k] for k in range(3)])

    points, conn = _lattice_py(dims, lo, spacing)
    cells = [] if conn is None else [("hexahedron", conn)]
    num_points = points.shape[0]
    num_cells = 0 if conn is None else conn.shape[0]

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
                    f"VTI {what} array '{name}' has {arr.shape[0]} rows, but the "
                    f"extent has {expected} {what}s"
                )
            if what == "cell":
                if num_cells == 0:
                    continue
                sink[name] = [arr]
            else:
                sink[name] = arr

    return Mesh(points, cells, point_data=point_data, cell_data=cell_data,
                field_data=field_arrays(piece, reader, "VTI"))


def _encode_binary(data, compression, header_type):
    """One ``<DataArray>`` body as base64, in VTK's block framing.

    Transcribed from :mod:`meshioplusplus.vtu`'s two writer closures (whose form
    is tied to its streaming-etree writer and cannot be called from here). The
    cross-check that it agrees is not this docstring: ``tests/python/test_vti.py``
    reads every encoding back through the *C++* reader and vice versa.
    """
    hdr_dtype = vtu_to_numpy_type[header_type]
    data_bytes = data.tobytes()
    if not compression:
        header = np.array(len(data_bytes), dtype=hdr_dtype)
        return base64.b64encode(header.tobytes() + data_bytes).decode()

    max_block_size = 32768
    num_blocks = -int(-len(data_bytes) // max_block_size)
    last_block_size = len(data_bytes) - (num_blocks - 1) * max_block_size
    c = _compressor_for(_COMPRESSION_TO_ATTR[compression])
    blocks = [c.compress(b) for b in _chunk_it(data_bytes, max_block_size)]
    header = np.array(
        [num_blocks, max_block_size, last_block_size] + [len(b) for b in blocks],
        dtype=hdr_dtype,
    )
    return (
        base64.b64encode(header.tobytes()).decode()
        + base64.b64encode(b"".join(blocks)).decode()
    )


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    spec = lattice_from_mesh(mesh)
    if spec is None:
        raise WriteError(
            "ImageData is a regular lattice, and this mesh is not one: it needs "
            "exactly one hexahedron block whose points tile an axis-aligned box "
            "with uniform spacing. A partial grid (voxelize's 'surface'/'inside' "
            "fill, or an octree) cannot be written as .vti either -- write it as "
            ".vtu, which stores the cells explicitly."
        )
    dims, origin, spacing = spec
    if header_type is None:
        header_type = "UInt32"
    if compression not in (None, "zlib", "lzma", "lz4", "zstd"):
        raise WriteError(f"Unknown VTI compression '{compression}'")
    if not binary:
        compression = None

    ext = f"0 {int(dims[0])} 0 {int(dims[1])} 0 {int(dims[2])}"
    lines = ['<?xml version="1.0"?>']
    attrs = 'type="ImageData" version="0.1" byte_order="LittleEndian"'
    if compression:
        attrs += f' compressor="{_COMPRESSION_TO_ATTR[compression]}"'
    if header_type != "UInt32":
        attrs += f' header_type="{header_type}"'
    lines.append(f"<VTKFile {attrs}>")
    lines.append(_provenance.render_xml_comment(_provenance.SlotTier.BLOCK))

    # %.17g, matching the C++ writer exactly: the stream default of six
    # significant digits would lose ~10 digits of a real origin, placing the grid
    # 1e-7 off its own points with nothing downstream to flag it.
    def num(v):
        return f"{float(v):.17g}"

    lines.append(
        f'<ImageData WholeExtent="{ext}" '
        f'Origin="{num(origin[0])} {num(origin[1])} {num(origin[2])}" '
        f'Spacing="{num(spacing[0])} {num(spacing[1])} {num(spacing[2])}">'
    )
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

    emit("PointData", sorted(mesh.point_data.items()))
    # A lattice has exactly one block, so cell_data's per-block list has one
    # entry; concatenating anyway keeps a mis-shaped input from being written
    # silently truncated.
    emit(
        "CellData",
        [
            (k, np.concatenate([np.asarray(a) for a in v]))
            for k, v in sorted(mesh.cell_data.items())
        ],
    )

    lines.append("</Piece>")
    lines.append("</ImageData>")
    lines.append("</VTKFile>")
    with open(filename, "w") as f:
        f.write("\n".join(lines) + "\n")
