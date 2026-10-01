"""Shared, private VTK XML framing for the PolyData and structured readers."""

import math

import numpy as np

from ._common import warn
from ._exceptions import ReadError
from ._mesh import Mesh
from ._region_field_data import is_region_field_name, regions_from_field_arrays
from ._regions import Region
from .vtu._vtu import VtuReader, _compressor_for, _load_root


class ArrayReader:
    """Reuse VTU's decoder rather than maintaining another binary framing path."""

    _ordered = VtuReader._ordered
    read_appended = VtuReader.read_appended
    read_uncompressed_binary = VtuReader.read_uncompressed_binary
    read_compressed_binary = VtuReader.read_compressed_binary
    read_data = VtuReader.read_data

    def __init__(
        self,
        header_type,
        byte_order,
        compression,
        appended_data=None,
        raw_appended=None,
    ):
        if header_type not in ("UInt32", "UInt64"):
            raise ReadError(f"Unknown VTK XML header type '{header_type}'")
        if byte_order not in (None, "LittleEndian", "BigEndian"):
            raise ReadError(f"Unknown VTK XML byte order '{byte_order}'")
        if compression is not None:
            _compressor_for(compression)
        self.header_type = header_type
        self.byte_order = byte_order
        self.compression = compression
        self.appended_data = appended_data
        self.raw_appended = raw_appended


def load(filename, dataset_type, format_name):
    root, raw = _load_root(filename)
    if root.tag != "VTKFile":
        raise ReadError("Expected tag 'VTKFile'")
    if root.get("type") != dataset_type:
        raise ReadError(f"Expected type {dataset_type}")
    compression = root.get("compressor")
    # Preserve the structured readers' intentional lzma parity contract.
    # PolyData still has a Python-only lzma path, as it did before.
    if compression == "vtkLZMADataCompressor" and dataset_type != "PolyData":
        raise ReadError(f"lzma-compressed {format_name} is not supported")
    text = None
    appended = root.find("AppendedData")
    if appended is not None and raw is None:
        encoding = appended.get("encoding", "base64")
        if encoding != "base64":
            raise ReadError(f"Unknown {format_name} AppendedData encoding '{encoding}'")
        text = "".join((appended.text or "").split())
        if not text.startswith("_"):
            raise ReadError(f"{format_name}: AppendedData does not start with '_'")
        text = text[1:]
    reader = ArrayReader(
        root.get("header_type", "UInt32"),
        root.get("byte_order"),
        compression,
        text,
        raw,
    )
    return root, reader


def piece_extent(grid, piece):
    """Validate a piece's own extent, without assuming it fills WholeExtent."""

    def parse(text):
        try:
            values = [int(x) for x in text.split()]
        except (AttributeError, ValueError) as exc:
            raise ReadError("VTK: malformed extent") from exc
        if len(values) != 6:
            raise ReadError("VTK: malformed extent")
        dims = [values[2 * k + 1] - values[2 * k] + 1 for k in range(3)]
        if min(dims) < 1 or math.prod(dims) > np.iinfo(np.int64).max:
            raise ReadError("VTK: inverted or overflowing extent")
        return np.array(values, dtype=np.int64)

    whole = parse(grid.get("WholeExtent"))
    text = piece.get("Extent")
    if text is None and len(grid.findall("Piece")) > 1:
        raise ReadError("VTK: missing Piece Extent")
    extent = whole if text is None else parse(text)
    if np.any(extent[::2] < whole[::2]) or np.any(extent[1::2] > whole[1::2]):
        raise ReadError("VTK: Piece Extent is outside WholeExtent")
    return extent


def field_arrays(holder, reader, fmt):
    from .vtu._vtu import vtu_to_numpy_type

    out = {}
    for section in holder.findall("FieldData"):
        for da in section.findall("DataArray"):
            if da.get("type") not in vtu_to_numpy_type:
                warn(f"{fmt}: skipping non-numeric field array '{da.get('Name')}'")
                continue
            out[da.get("Name")] = reader.read_data(da)
    return out


def read_pieces(root, reader, dataset_type, fmt, read_piece):
    """Piece-major concatenation; never weld boundaries or discard ghost cells."""
    grid = root.find(dataset_type)
    if grid is None:
        raise ReadError(f"No {dataset_type} found")
    pieces = grid.findall("Piece")
    if not pieces:
        raise ReadError("No Piece found")
    meshes = [read_piece(grid, piece, reader) for piece in pieces]
    if dataset_type != "PolyData":
        for mesh in meshes:
            count = sum(len(cb.data) for cb in mesh.cells)
            mesh.regions, mesh.field_data = regions_from_field_arrays(
                mesh.field_data, len(mesh.points), np.arange(count), fmt
            )
    points = [mesh.points for mesh in meshes]
    if any(p.dtype != points[0].dtype or p.shape[1:] != points[0].shape[1:] for p in points):
        raise ReadError("VTK: pieces disagree on point dtype or components")
    point_data, cell_raw = {}, {}
    for attr, sink in (("point_data", point_data), ("cell_data", cell_raw)):
        names = sorted(set().union(*(getattr(mesh, attr) for mesh in meshes)))
        for name in names:
            parts = []
            for mesh in meshes:
                data = getattr(mesh, attr)
                if name not in data:
                    break
                if attr == "cell_data":
                    parts.extend(data[name])
                else:
                    parts.append(data[name])
            else:
                if parts and all(p.dtype == parts[0].dtype and p.shape[1:] == parts[0].shape[1:] for p in parts):
                    sink[name] = np.concatenate(parts)
                    continue
            warn(f"{fmt}: data '{name}' is missing from, or differs between, pieces; dropped")

    cells, file_maps, local_regions = [], [], {}
    point_base = cell_base = 0
    fields = field_arrays(grid, reader, fmt)
    passthrough = {}
    for mesh in meshes:
        count = sum(len(cb.data) for cb in mesh.cells)
        file_maps.append(getattr(mesh, "_vtk_file_to_global", np.arange(count)) + cell_base)
        for cb in mesh.cells:
            data = cb.data + point_base
            if cells and cells[-1][0] == cb.type and cells[-1][1].shape[1:] == data.shape[1:]:
                cells[-1] = (cb.type, np.concatenate([cells[-1][1], data]))
            else:
                cells.append((cb.type, data))
        for name, arr in mesh.field_data.items():
            if is_region_field_name(name):
                passthrough[name] = arr
            else:
                fields[name] = arr
        for region in mesh.regions:
            entries = region.entries.copy()
            if region.kind == "point":
                entries += point_base
            elif region.kind == "cell":
                entries += cell_base
            else:
                entries[:, 0] += cell_base
            key = region.kind, region.name
            if key in local_regions:
                entries = np.concatenate([local_regions[key].entries, entries])
            local_regions[key] = Region(region.name, region.kind, entries, region.dim, region.tag)
        point_base += len(mesh.points)
        cell_base += count
    regions, fields = regions_from_field_arrays(fields, point_base, np.concatenate(file_maps), fmt)
    fields.update(passthrough)
    all_regions = {(r.kind, r.name): r for r in regions}
    all_regions.update(local_regions)
    cell_data = {}
    boundaries = np.cumsum([0] + [len(data) for _, data in cells])
    for name, arr in cell_raw.items():
        cell_data[name] = [arr[a:b].copy() for a, b in zip(boundaries[:-1], boundaries[1:])]
    return Mesh(np.concatenate(points), cells, point_data=point_data, cell_data=cell_data,
                field_data=fields, regions=list(all_regions.values()))
