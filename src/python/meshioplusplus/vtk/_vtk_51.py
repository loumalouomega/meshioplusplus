import functools

import numpy as np

from .._common import info, join_strings, replace_space
from .._files import open_file
from .._vtk_common import meshio_to_vtk_order, vtk_cells_from_data
from . import _legacy_common as _legacy

# VTK 5.1 data types
vtk_to_numpy_dtype_name = {
    "float": "float32",
    "double": "float64",
    "int": "int32",
    "char": "int8",
    "short": "int16",
    "long": "int64",
    "unsigned_char": "uint8",
    "unsigned_short": "uint16",
    "unsigned_int": "uint32",
    "unsigned_long": "uint64",
    "vtktypeint8": "int8",
    "vtktypeint16": "int16",
    "vtktypeint32": "int32",
    "vtktypeint64": "int64",
    "vtktypeuint8": "uint8",
    "vtktypeuint16": "uint16",
    "vtktypeuint32": "uint32",
    "vtktypeuint64": "uint64",
}

_read_points = functools.partial(
    _legacy._read_points, dtype_names=vtk_to_numpy_dtype_name
)
numpy_to_vtk_dtype = {v: k for k, v in vtk_to_numpy_dtype_name.items()}


def read(filename):
    with open_file(filename, "rb") as f:
        out = read_buffer(f)
    return out


def read_buffer(f):
    return _legacy.read_buffer(
        f, _read_section, _check_mesh, _cells_from_info, vtk_to_numpy_dtype_name
    )


def _cells_from_info(info):
    return vtk_cells_from_data(
        info.connectivity, info.offsets, info.types, info.cell_data_raw
    )


def _read_section(f, info):
    if _legacy.read_common_section(f, info, _read_points, ", "):
        return
    if info.section == "CELLS":
        info.active = "CELLS"
        try:
            line = f.readline().decode()
        except UnicodeDecodeError:
            line = ""

        assert line.startswith("OFFSETS")
        # vtk DataFile Version 5.1 - appearing in Paraview 5.8.1 outputs
        # No specification found for this file format.
        # See the question on ParaView Discourse Forum:
        # https://discourse.paraview.org/t/specification-of-vtk-datafile-version-5-1/5127
        info.num_offsets = int(info.split[1])
        info.num_items = int(info.split[2])
        dtype = np.dtype(vtk_to_numpy_dtype_name[line.split()[1]])
        offsets = _legacy.read_int_data(
            f, info.is_ascii, info.num_offsets, dtype, "Expected newline"
        )

        line = f.readline().decode()
        assert line.startswith("CONNECTIVITY")
        dtype = np.dtype(vtk_to_numpy_dtype_name[line.split()[1]])
        connectivity = _legacy.read_int_data(
            f, info.is_ascii, info.num_items, dtype, "Expected newline"
        )
        info.connectivity = connectivity
        assert offsets[0] == 0
        assert offsets[-1] == len(connectivity)
        info.offsets = offsets[1:]


def _check_mesh(info):
    _legacy.check_mesh(info)
    if info.dataset["type"] != "UNSTRUCTURED_GRID":
        # Structured generators return the counted legacy layout, not 5.1 CSR.
        width = {1: 1, 3: 2, 9: 4, 12: 8}[int(info.types[0])]
        info.connectivity = info.connectivity.reshape(-1, width + 1)[:, 1:].reshape(-1)
        info.offsets = np.arange(1, len(info.types) + 1, dtype=np.int64) * width


def write(filename, mesh, binary=True):
    points = _legacy.pad_2d_to_3d(mesh, binary)

    if mesh.point_sets:
        info(
            "VTK format cannot write point_sets. Converting them to point_data...",
            highlight=False,
        )
        key, _ = join_strings(list(mesh.point_sets.keys()))
        key, _ = replace_space(key)
        mesh.point_sets_to_data(key)

    if mesh.cell_sets:
        info(
            "VTK format cannot write cell_sets. Converting them to cell_data...",
            highlight=False,
        )
        key, _ = join_strings(list(mesh.cell_sets.keys()))
        key, _ = replace_space(key)
        mesh.cell_sets_to_data(key)

    with open_file(filename, "wb") as f:
        _legacy.write_unstructured(
            f, "5.1", mesh, points, binary, _write_cells, numpy_to_vtk_dtype
        )


def _write_cells(f, cells, binary):
    total_num_cells = sum(len(c.data) for c in cells)
    total_num_idx = sum(c.data.size for c in cells)
    f.write(f"CELLS {total_num_cells + 1} {total_num_idx}\n".encode())

    # offsets
    offsets = [[0]]
    k = 0
    for cell_block in cells:
        m, n = cell_block.data.shape
        offsets.append(np.arange(k + n, k + (m + 1) * n, n))
        k = offsets[-1][-1]
    offsets = np.concatenate(offsets)

    if binary:
        f.write(b"OFFSETS vtktypeint64\n")
        # force big-endian and int64
        offsets.astype(">i8").tofile(f, sep="")
        f.write(b"\n")

        f.write(b"CONNECTIVITY vtktypeint64\n")
        for cell_block in cells:
            d = cell_block.data
            cell_idx = meshio_to_vtk_order(cell_block.type)
            if cell_idx is not None:
                d = d[:, cell_idx]
            # force big-endian and int64
            d.astype(">i8").tofile(f, sep="")
        f.write(b"\n")
    else:
        # ascii
        f.write(b"OFFSETS vtktypeint64\n")
        offsets.tofile(f, sep="\n")
        f.write(b"\n")

        f.write(b"CONNECTIVITY vtktypeint64\n")
        for cell_block in cells:
            d = cell_block.data
            cell_idx = meshio_to_vtk_order(cell_block.type)
            if cell_idx is not None:
                d = d[:, cell_idx]
            d.tofile(f, sep="\n")
            f.write(b"\n")

    _legacy.write_cell_types(f, cells, binary)


def _write_field_data(f, data, binary):
    _legacy.write_field_data(f, data, binary, numpy_to_vtk_dtype)
