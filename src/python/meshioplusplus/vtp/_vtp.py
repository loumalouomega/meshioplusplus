"""
I/O for VTK XML PolyData (.vtp) files, c.f.
<https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html>

The same VTK-XML container as VTU with a <PolyData> grid holding
<Verts>/<Lines>/<Polys>/<Strips> connectivity+offsets sections. Only surface
cells are representable: ``vertex`` (Verts), ``line`` (Lines) and
``triangle``/``quad``/``polygon`` (Polys). Cell data follows VTK's canonical
PolyData cell order — Verts, Lines, Polys, Strips — in both directions.
Triangle strips, poly-vertex/poly-line rows and multiple pieces are not
supported.
"""

import base64
import lzma
import zlib
from xml.etree import ElementTree as ET

import numpy as np

from .. import _provenance
from .._exceptions import ReadError, WriteError
from .._mesh import Mesh
from .._vtk_common import vtk_cells_from_data
from ..vtu._vtu import numpy_to_vtu_type, vtu_to_numpy_type

_SECTION_TAGS = ("Verts", "Lines", "Polys", "Strips")


def _read_binary_data(data, dtype, byte_order, header_type, compression):
    header_dtype = vtu_to_numpy_type[header_type]
    if byte_order is not None:
        bo = "<" if byte_order == "LittleEndian" else ">"
        header_dtype = header_dtype.newbyteorder(bo)
        dtype = dtype.newbyteorder(bo)

    if compression is None:
        byte_string = base64.b64decode(data)
        num_header_bytes = np.dtype(header_dtype).itemsize
        total_num_bytes = np.frombuffer(byte_string[:num_header_bytes], header_dtype)[0]
        # The header may have been base64-encoded separately (padding).
        if len(byte_string) == num_header_bytes:
            header_len = len(base64.b64encode(byte_string))
            byte_string = base64.b64decode(data[header_len:])
        else:
            byte_string = byte_string[num_header_bytes:]
        return np.frombuffer(byte_string[:total_num_bytes], dtype=dtype)

    # compressed: header = [num_blocks, max_block_size, last_block_size, sizes...]
    num_bytes_per_item = np.dtype(header_dtype).itemsize
    num_chars = -(-num_bytes_per_item // 3) * 4  # base64 chars for the first item
    byte_string = base64.b64decode(data[:num_chars])[:num_bytes_per_item]
    num_blocks = int(np.frombuffer(byte_string, header_dtype)[0])

    num_header_bytes = num_bytes_per_item * (3 + num_blocks)
    num_header_chars = -(-num_header_bytes // 3) * 4
    header = np.frombuffer(base64.b64decode(data[:num_header_chars]), header_dtype)
    block_sizes = header[3:]

    byte_array = base64.b64decode(data[num_header_chars:])
    byte_offsets = np.concatenate([[0], np.cumsum(block_sizes)]).astype(np.int64)

    c = {"vtkLZMADataCompressor": lzma, "vtkZLibDataCompressor": zlib}[compression]
    return np.concatenate(
        [
            np.frombuffer(
                c.decompress(byte_array[byte_offsets[k] : byte_offsets[k + 1]]),
                dtype=dtype,
            )
            for k in range(num_blocks)
        ]
    )


def _read_data_array(elem, byte_order, header_type, compression):
    fmt = elem.get("format", "ascii")
    dtype = vtu_to_numpy_type[elem.get("type")]
    num_components = int(elem.get("NumberOfComponents", 0))

    if fmt == "ascii":
        data = np.array((elem.text or "").split(), dtype=dtype)
    elif fmt == "binary":
        data = _read_binary_data(
            (elem.text or "").strip(), dtype, byte_order, header_type, compression
        )
    else:
        raise ReadError(f"VTP '{fmt}' data is not supported")

    if num_components > 1:
        data = data.reshape(-1, num_components)
    return data


def read(filename):
    tree = ET.parse(filename)
    root = tree.getroot()
    if root.tag != "VTKFile":
        raise ReadError("Expected tag 'VTKFile'")
    if root.get("type") != "PolyData":
        raise ReadError("Expected type PolyData")

    byte_order = root.get("byte_order")
    compression = root.get("compressor")
    if compression not in (None, "vtkZLibDataCompressor", "vtkLZMADataCompressor"):
        raise ReadError(f"Unknown VTP compressor '{compression}'")
    header_type = root.get("header_type", "UInt32")

    if root.find("AppendedData") is not None:
        raise ReadError("appended VTP data is not supported")

    grid = root.find("PolyData")
    if grid is None:
        raise ReadError("No PolyData found")
    pieces = grid.findall("Piece")
    if len(pieces) != 1:
        raise ReadError("Only single-piece PolyData is supported")
    piece = pieces[0]

    def read_data(elem):
        return _read_data_array(elem, byte_order, header_type, compression)

    points = None
    point_data = {}
    cell_data_raw = {}
    sections = {}
    for child in piece:
        if child.tag == "Points":
            points = read_data(child.find("DataArray"))
            if points.ndim == 1:
                points = points.reshape(-1, 3)
        elif child.tag == "PointData":
            for da in child.findall("DataArray"):
                point_data[da.get("Name")] = read_data(da)
        elif child.tag == "CellData":
            for da in child.findall("DataArray"):
                cell_data_raw[da.get("Name")] = read_data(da)
        elif child.tag in _SECTION_TAGS:
            arrays = {}
            for da in child.findall("DataArray"):
                arrays[da.get("Name")] = read_data(da)
            sections[child.tag] = (
                np.asarray(arrays.get("connectivity", []), dtype=np.int64),
                np.asarray(arrays.get("offsets", []), dtype=np.int64),
            )

    if "Strips" in sections and sections["Strips"][1].size > 0:
        raise ReadError("triangle-strip VTP cells are not supported")

    # Concatenate sections in VTK's canonical PolyData cell order (Verts,
    # Lines, Polys), synthesizing a VTK type id per row so the shared VTK
    # reconstruction can build the blocks and split cell_data.
    conn_parts = []
    offset_parts = []
    type_parts = []
    conn_base = 0
    for tag, kind in (("Verts", 0), ("Lines", 1), ("Polys", 2)):
        if tag not in sections:
            continue
        conn, offsets = sections[tag]
        if offsets.size == 0:
            continue
        sizes = np.diff(np.concatenate([[0], offsets]))
        if kind == 0:
            if not np.all(sizes == 1):
                raise ReadError("poly-vertex VTP cells are not supported")
            types = np.full(sizes.shape, 1, dtype=np.int64)  # VTK_VERTEX
        elif kind == 1:
            if not np.all(sizes == 2):
                raise ReadError("poly-line VTP cells are not supported")
            types = np.full(sizes.shape, 3, dtype=np.int64)  # VTK_LINE
        else:
            types = np.where(sizes == 3, 5, np.where(sizes == 4, 9, 7))
        conn_parts.append(conn)
        offset_parts.append(offsets + conn_base)
        type_parts.append(types)
        conn_base += conn.size

    if conn_parts:
        connectivity = np.concatenate(conn_parts)
        offsets = np.concatenate(offset_parts)
        types = np.concatenate(type_parts)
        cells, cell_data = vtk_cells_from_data(
            connectivity, offsets, types, cell_data_raw
        )
    else:
        cells, cell_data = [], {}

    if points is None:
        points = np.empty((0, 3))
    return Mesh(points, cells, point_data=point_data, cell_data=cell_data)


def _chunk_it(array, n):
    k = 0
    while len(array[k * n : (k + 1) * n]) > 0:
        yield array[k * n : (k + 1) * n]
        k += 1


def _classify_block(cell_block):
    cell_type = cell_block.type
    if cell_type == "vertex":
        return 0
    if cell_type == "line":
        return 1
    if cell_type in ("triangle", "quad", "polygon"):
        return 2
    raise WriteError(f"VTP: PolyData cannot hold '{cell_type}' cells")


def _block_rows(cell_block):
    data = cell_block.data
    if isinstance(data, np.ndarray) and data.ndim == 2:
        return [np.asarray(row, dtype=np.int64) for row in data]
    return [np.asarray(row, dtype=np.int64) for row in data]


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    if header_type is None:
        header_type = "UInt32"
    if compression not in (None, "zlib", "lzma"):
        raise WriteError(f"Unknown compression '{compression}'")

    # Classify blocks and build the VTK canonical order (Verts, Lines, Polys)
    # as a stable partition of the mesh's block order.
    kinds = [_classify_block(cb) for cb in mesh.cells]
    block_order = [
        bi for want in (0, 1, 2) for bi, kind in enumerate(kinds) if kind == want
    ]
    section_rows = {0: [], 1: [], 2: []}
    for bi in block_order:
        rows = _block_rows(mesh.cells[bi])
        if kinds[bi] == 0 and any(len(r) != 1 for r in rows):
            raise WriteError("VTP: vertex cells must have exactly one node")
        section_rows[kinds[bi]].extend(rows)

    points = np.asarray(mesh.points)
    if points.shape[1] < 3:
        points = np.column_stack(
            [points, np.zeros((points.shape[0], 3 - points.shape[1]), points.dtype)]
        )

    def data_array_str(name, data, ncomp):
        vtu_type = numpy_to_vtu_type[data.dtype.newbyteorder("=")]
        out = [f'<DataArray type="{vtu_type}" Name="{name}"']
        if ncomp > 0:
            out.append(f' NumberOfComponents="{ncomp}"')
        out.append(f' format="{"binary" if binary else "ascii"}">\n')
        flat = np.ascontiguousarray(data).astype(
            data.dtype.newbyteorder("="), copy=False
        )
        if binary:
            data_bytes = flat.tobytes()
            if compression:
                max_block_size = 32768
                num_blocks = -(-len(data_bytes) // max_block_size)
                c = {"lzma": lzma, "zlib": zlib}[compression]
                blocks = [c.compress(b) for b in _chunk_it(data_bytes, max_block_size)]
                last_block_size = len(data_bytes) - (num_blocks - 1) * max_block_size
                header = np.array(
                    [num_blocks, max_block_size, last_block_size]
                    + [len(b) for b in blocks],
                    dtype=vtu_to_numpy_type[header_type],
                )
                out.append(base64.b64encode(header.tobytes()).decode())
                out.append(base64.b64encode(b"".join(blocks)).decode())
            else:
                header = np.array(len(data_bytes), dtype=vtu_to_numpy_type[header_type])
                out.append(base64.b64encode(header.tobytes() + data_bytes).decode())
            out.append("\n")
        else:
            fmt = "{:.11e}\n" if vtu_type.startswith("Float") else "{:d}\n"
            out.append("".join(fmt.format(v) for v in flat.reshape(-1)))
        out.append("</DataArray>\n")
        return "".join(out)

    def section_str(tag, rows):
        if not rows:
            return ""
        connectivity = np.concatenate(rows).astype(np.int64)
        offsets = np.cumsum([len(r) for r in rows]).astype(np.int64)
        return (
            f"<{tag}>\n"
            + data_array_str("connectivity", connectivity, 0)
            + data_array_str("offsets", offsets, 0)
            + f"</{tag}>\n"
        )

    out = []
    out.append('<?xml version="1.0"?>\n')
    out.append('<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian"')
    if binary and compression:
        compressor = {
            "zlib": "vtkZLibDataCompressor",
            "lzma": "vtkLZMADataCompressor",
        }[compression]
        out.append(f' compressor="{compressor}"')
    if header_type != "UInt32":
        out.append(f' header_type="{header_type}"')
    out.append(">\n")
    out.append(_provenance.render_xml_comment(_provenance.SlotTier.BLOCK) + "\n")
    out.append("<PolyData>\n")
    out.append(
        f'<Piece NumberOfPoints="{points.shape[0]}"'
        f' NumberOfVerts="{len(section_rows[0])}"'
        f' NumberOfLines="{len(section_rows[1])}"'
        f' NumberOfStrips="0"'
        f' NumberOfPolys="{len(section_rows[2])}">\n'
    )
    out.append("<Points>\n")
    out.append(data_array_str("Points", points, 3))
    out.append("</Points>\n")
    out.append(section_str("Verts", section_rows[0]))
    out.append(section_str("Lines", section_rows[1]))
    out.append(section_str("Polys", section_rows[2]))

    if mesh.point_data:
        out.append("<PointData>\n")
        for name in sorted(mesh.point_data):
            data = np.asarray(mesh.point_data[name])
            ncomp = data.shape[1] if data.ndim == 2 else 0
            out.append(data_array_str(name, data, ncomp))
        out.append("</PointData>\n")

    if mesh.cell_data:
        # Cell data follows the reordered (Verts, Lines, Polys) block order.
        out.append("<CellData>\n")
        for name in sorted(mesh.cell_data):
            blocks = mesh.cell_data[name]
            data = np.concatenate([np.asarray(blocks[bi]) for bi in block_order])
            ncomp = data.shape[1] if data.ndim == 2 else 0
            out.append(data_array_str(name, data, ncomp))
        out.append("</CellData>\n")

    out.append("</Piece>\n</PolyData>\n</VTKFile>\n")
    payload = "".join(out).encode()

    if hasattr(filename, "write"):
        try:
            filename.write(payload)
        except TypeError:
            filename.write(payload.decode())
    else:
        with open(filename, "wb") as f:
            f.write(payload)
