"""
I/O for VTK XML PolyData (.vtp) files, c.f.
<https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html>

The same VTK-XML container as VTU with a <PolyData> grid holding
<Verts>/<Lines>/<Polys>/<Strips> connectivity+offsets sections. Only surface
cells are representable: ``vertex`` (Verts), ``line`` (Lines) and
``triangle``/``quad``/``polygon`` (Polys). Cell data follows VTK's canonical
PolyData cell order — Verts, Lines, Polys, Strips — in both directions.
Multiple pieces concatenate in document order without welding. Triangle
strips and poly-vertex/poly-line rows are not supported.
"""

import base64
import lzma
import zlib

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._mesh import Mesh
from .._region_field_data import (
    FILE_INDEX_KEY,
    file_to_global_from,
    is_region_field_name,
    regions_from_field_arrays,
    regions_to_field_arrays,
)
from .._vtk_common import vtk_cells_from_data
from ..vtu._vtu import numpy_to_vtu_type, vtu_to_numpy_type

_SECTION_TAGS = ("Verts", "Lines", "Polys", "Strips")


def read(filename):
    from .._vtk_xml_read import load, read_pieces

    root, reader = load(filename, "PolyData", "VTP")
    return read_pieces(root, reader, "PolyData", "VTP", _read_piece)


def _read_piece(grid, piece, reader):
    def read_data(elem):
        data = reader.read_data(elem)
        # PolyData historically exposes scalar arrays as one-dimensional.
        if elem.get("NumberOfComponents", "0") == "1":
            data = data.reshape(-1)
        return data

    points = None
    point_data = {}
    field_data = {}
    cell_data_raw = {}
    sections = {}

    # <FieldData> belongs to the dataset: VTK writes it on the grid, before the
    # <Piece>, and also accepts it inside one (the piece's overriding the grid's).
    # A non-numeric array (`type="String"`) has no numpy dtype here: skipped with
    # a warning rather than failing a file whose field data used to be ignored.
    for holder in (piece,):
        for fd in holder.findall("FieldData"):
            for da in fd.findall("DataArray"):
                if da.get("type") not in vtu_to_numpy_type:
                    warn(
                        f"VTP: skipping <FieldData> array '{da.get('Name')}' of "
                        f"type '{da.get('type')}' (only numeric arrays are read)"
                    )
                    continue
                field_data[da.get("Name")] = read_data(da)

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

    def count(name):
        text = piece.get(name, "")
        if not text or any(ch not in "0123456789" for ch in text):
            raise ReadError(f"VTP: invalid {name}")
        value = int(text)
        if value > np.iinfo(np.int64).max:
            raise ReadError(f"VTP: invalid {name}")
        return value

    num_points = count("NumberOfPoints")
    for tag in _SECTION_TAGS:
        name = "NumberOf" + tag
        actual = sections[tag][1].size if tag in sections else 0
        if piece.get(name) is not None and count(name) != actual:
            raise ReadError("VTP: cell count differs from section offsets")
    if points is not None and len(points) != num_points:
        raise ReadError("VTP Points length differs from NumberOfPoints")
    if any(len(arr) != num_points for arr in point_data.values()):
        raise ReadError("VTP PointData length differs from NumberOfPoints")

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
        if offsets.size and (
            np.any(np.diff(np.concatenate([[0], offsets])) <= 0)
            or offsets[-1] != conn.size
        ):
            raise ReadError("VTP: offsets do not span connectivity")
        if conn.size and (np.min(conn) < 0 or np.max(conn) >= num_points):
            raise ReadError("VTP: point index out of range")
        if offsets.size == 0:
            if conn.size:
                raise ReadError("VTP: offsets do not span connectivity")
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
            if np.any(sizes < 3):
                raise ReadError("VTP: polygon has fewer than three points")
            types = np.where(sizes == 3, 5, np.where(sizes == 4, 9, 7))
        conn_parts.append(conn)
        offset_parts.append(offsets + conn_base)
        type_parts.append(types)
        conn_base += conn.size

    num_cells = sum(len(types) for types in type_parts)
    if any(len(arr) != num_cells for arr in cell_data_raw.values()):
        raise ReadError("VTP CellData length differs from cell count")
    if points is None and num_points:
        raise ReadError("VTP Piece has no Points")

    file_to_global = None
    if conn_parts:
        connectivity = np.concatenate(conn_parts)
        offsets = np.concatenate(offset_parts)
        types = np.concatenate(type_parts)
        # The hidden file index follows the cells (see _region_field_data).
        cell_data_raw[FILE_INDEX_KEY] = np.arange(len(types), dtype=np.int64)
        cells, cell_data = vtk_cells_from_data(
            connectivity, offsets, types, cell_data_raw
        )
        file_to_global = file_to_global_from(cell_data)
    else:
        cells, cell_data = [], {}

    if points is None:
        points = np.empty((0, 3))
    regions, field_data = regions_from_field_arrays(
        field_data,
        len(points),
        file_to_global if file_to_global is not None else np.empty(0, np.int64),
        "vtp",
    )
    mesh = Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
    )
    if regions:
        mesh.regions = regions
    mesh._vtk_file_to_global = (
        file_to_global if file_to_global is not None else np.empty(0, np.int64)
    )
    return mesh


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

    def data_array_str(name, data, ncomp, ntuples=None):
        if name == "vtkGhostType" and data.dtype != np.uint8:
            # VTK's reserved ghost-flag name: always UInt8 on disk (see the VTU writer).
            data = data.astype(np.uint8)
        vtu_type = numpy_to_vtu_type[data.dtype.newbyteorder("=")]
        out = [f'<DataArray type="{vtu_type}" Name="{name}"']
        if ntuples is not None:
            # <FieldData> arrays carry an explicit tuple count (VTK requires it).
            out.append(f' NumberOfTuples="{ntuples}"')
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
    # Field data is dataset-global: on the grid, before the <Piece>, where VTK's
    # own writers put it. A value that is not a numeric array has no VTK type and
    # is skipped with a warning; the caller's mapping is never modified.
    field_arrays = {}
    for name in sorted(mesh.field_data):
        if is_region_field_name(name):
            warn(
                f"VTP: field_data '{name}' uses the region naming convention; not written"
            )
            continue
        try:
            arr = np.asarray(mesh.field_data[name])
            arr = arr.astype(arr.dtype.newbyteorder("="), copy=False)
        except Exception:
            arr = None
        if arr is None or arr.dtype not in numpy_to_vtu_type:
            warn(f"VTP: field_data '{name}' is not a numeric array; not written")
            continue
        field_arrays[name] = arr.reshape(arr.shape[0], -1) if arr.ndim > 2 else arr
    # Named regions, their cells numbered in the file's (Verts, Lines, Polys) order.
    if getattr(mesh, "regions", None):
        sizes = [len(cb.data) for cb in mesh.cells]
        bases = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
        global_to_file = np.empty(int(bases[-1]), dtype=np.int64)
        next_file = 0
        for bi in block_order:
            n = sizes[bi]
            global_to_file[bases[bi] : bases[bi] + n] = np.arange(
                next_file, next_file + n
            )
            next_file += n
        field_arrays.update(regions_to_field_arrays(mesh, global_to_file))
    if field_arrays:
        out.append("<FieldData>\n")
        for name, arr in field_arrays.items():
            ncomp = arr.shape[1] if arr.ndim == 2 else 0
            out.append(
                data_array_str(name, arr, ncomp, arr.shape[0] if arr.ndim else 1)
            )
        out.append("</FieldData>\n")
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
