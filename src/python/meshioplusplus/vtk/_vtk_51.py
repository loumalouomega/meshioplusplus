import functools

import numpy as np

from .. import _provenance
from .._common import info, join_strings, replace_space, warn
from .._exceptions import ReadError
from .._files import open_file
from .._vtk_common import meshio_to_vtk_order, meshio_to_vtk_type, vtk_cells_from_data
from . import _legacy_common as _legacy
from ._legacy_common import (
    _generate_cells,
    _generate_points,
    _read_cell_types,
    _skip_meta,
    vtk_dataset_types,
)

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
    if info.section == "METADATA":
        _skip_meta(f)

    elif info.section == "DATASET":
        info.active = "DATASET"
        info.dataset["type"] = info.split[1].upper()
        if info.dataset["type"] not in vtk_dataset_types:
            legal = ", ".join(vtk_dataset_types)
            raise ReadError(
                f"Only VTK '{legal}' supported (not {info.dataset['type']})."
            )

    elif info.section == "POINTS":
        info.active = "POINTS"
        info.num_points = int(info.split[1])
        data_type = info.split[2].lower()
        info.points = _read_points(f, data_type, info.is_ascii, info.num_points)

    elif info.section == "CELLS":
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
        offsets = _read_int_data(f, info.is_ascii, info.num_offsets, dtype)

        line = f.readline().decode()
        assert line.startswith("CONNECTIVITY")
        dtype = np.dtype(vtk_to_numpy_dtype_name[line.split()[1]])
        connectivity = _read_int_data(f, info.is_ascii, info.num_items, dtype)
        info.connectivity = connectivity
        assert offsets[0] == 0
        assert offsets[-1] == len(connectivity)
        info.offsets = offsets[1:]

    elif info.section == "CELL_TYPES":
        info.active = "CELL_TYPES"
        info.num_items = int(info.split[1])
        info.types = _read_cell_types(f, info.is_ascii, info.num_items)

    elif info.section == "POINT_DATA":
        info.active = "POINT_DATA"
        info.num_items = int(info.split[1])

    elif info.section == "CELL_DATA":
        info.active = "CELL_DATA"
        info.num_items = int(info.split[1])

    elif info.section == "LOOKUP_TABLE":
        info.num_items = int(info.split[2])
        np.fromfile(f, count=info.num_items * 4, sep=" ", dtype=float)
        # rgba = data.reshape((info.num_items, 4))

    elif info.section == "COLOR_SCALARS":
        nValues = int(info.split[2])
        # re-use num_items from active POINT/CELL_DATA
        num_items = info.num_items
        dtype = np.ubyte
        if info.is_ascii:
            dtype = float
        np.fromfile(f, count=num_items * nValues, dtype=dtype)


def _check_mesh(info):
    if info.dataset["type"] == "UNSTRUCTURED_GRID":
        if info.connectivity is None:
            raise ReadError("Required section CELLS not found.")
        if info.types is None:
            raise ReadError("Required section CELL_TYPES not found.")
    elif info.dataset["type"] == "STRUCTURED_POINTS":
        dim = info.dataset["DIMENSIONS"]
        ori = info.dataset["ORIGIN"]
        spa = (
            info.dataset["SPACING"]
            if "SPACING" in info.dataset
            else info.dataset["ASPECT_RATIO"]
        )
        axis = [
            np.linspace(ori[i], ori[i] + (dim[i] - 1.0) * spa[i], dim[i])
            for i in range(3)
        ]
        info.points = _generate_points(axis)
        info.connectivity, info.types = _generate_cells(dim=info.dataset["DIMENSIONS"])
    elif info.dataset["type"] == "RECTILINEAR_GRID":
        axis = [
            info.dataset["X_COORDINATES"],
            info.dataset["Y_COORDINATES"],
            info.dataset["Z_COORDINATES"],
        ]
        info.points = _generate_points(axis)
        info.connectivity, info.types = _generate_cells(dim=info.dataset["DIMENSIONS"])
    elif info.dataset["type"] == "STRUCTURED_GRID":
        info.connectivity, info.types = _generate_cells(dim=info.dataset["DIMENSIONS"])

    if info.dataset["type"] != "UNSTRUCTURED_GRID":
        # Structured generators return the counted legacy layout, not 5.1 CSR.
        width = {1: 1, 3: 2, 9: 4, 12: 8}[int(info.types[0])]
        info.connectivity = info.connectivity.reshape(-1, width + 1)[:, 1:].reshape(-1)
        info.offsets = np.arange(1, len(info.types) + 1, dtype=np.int64) * width


def _read_int_data(f, is_ascii, num_items, dtype):
    if is_ascii:
        c = np.fromfile(f, count=num_items, sep=" ", dtype=dtype)
    else:
        dtype = dtype.newbyteorder(">")
        c = np.fromfile(f, count=num_items, dtype=dtype)
        line = f.readline().decode()
        if line != "\n":
            raise ReadError("Expected newline")
    return c


def _pad(array):
    return np.pad(array, ((0, 0), (0, 1)), "constant")


def write(filename, mesh, binary=True):
    if mesh.points.shape[1] == 2:
        warn(
            "VTK requires 3D points, but 2D points given. "
            "Appending 0 third component."
        )
        points = _pad(mesh.points)
    else:
        points = mesh.points

    if mesh.point_data:
        for name, values in mesh.point_data.items():
            if len(values.shape) == 2 and values.shape[1] == 2:
                warn(
                    "VTK requires 3D vectors, but 2D vectors given. "
                    f"Appending 0 third component to {name}."
                )
                mesh.point_data[name] = _pad(values)

    for name, data in mesh.cell_data.items():
        for k, values in enumerate(data):
            if len(values.shape) == 2 and values.shape[1] == 2:
                warn(
                    "VTK requires 3D vectors, but 2D vectors given. "
                    f"Appending 0 third component to {name}."
                )
                data[k] = _pad(data[k])

    if not binary:
        warn("VTK ASCII files are only meant for debugging.")

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
        f.write(b"# vtk DataFile Version 5.1\n")
        f.write(f"{_provenance.lines(_provenance.SlotTier.SINGLE_LINE)[0]}\n".encode())
        f.write(("BINARY\n" if binary else "ASCII\n").encode())
        f.write(b"DATASET UNSTRUCTURED_GRID\n")

        # write points and cells
        _write_points(f, points, binary)
        _write_cells(f, mesh.cells, binary)

        # write point data
        if mesh.point_data:
            num_points = mesh.points.shape[0]
            f.write(f"POINT_DATA {num_points}\n".encode())
            _write_field_data(f, mesh.point_data, binary)

        # write cell data
        if mesh.cell_data:
            total_num_cells = sum(len(c.data) for c in mesh.cells)
            f.write(f"CELL_DATA {total_num_cells}\n".encode())
            _write_field_data(f, mesh.cell_data, binary)


def _write_points(f, points, binary):
    dtype = numpy_to_vtk_dtype[points.dtype.name]
    f.write(f"POINTS {len(points)} {dtype}\n".encode())

    if binary:
        # Binary data must be big endian, see
        # <https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        # if points.dtype.byteorder == "<" or (
        #     points.dtype.byteorder == "=" and sys.byteorder == "little"
        # ):
        #     logging.warn("Converting to new byte order")
        points.astype(points.dtype.newbyteorder(">")).tofile(f, sep="")
    else:
        # ascii
        points.tofile(f, sep=" ")
    f.write(b"\n")


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

    # write cell types
    f.write(f"CELL_TYPES {total_num_cells}\n".encode())
    if binary:
        for c in cells:
            vtk_type = meshio_to_vtk_type[c.type]
            np.full(len(c.data), vtk_type, dtype=np.dtype(">i4")).tofile(f, sep="")
        f.write(b"\n")
    else:
        # ascii
        for c in cells:
            vtk_type = meshio_to_vtk_type[c.type]
            np.full(len(c.data), vtk_type).tofile(f, sep="\n")
            f.write(b"\n")


def _write_field_data(f, data, binary):
    _legacy.write_field_data(f, data, binary, numpy_to_vtk_dtype)
