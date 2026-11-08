"""Readers the legacy VTK 4.2 and 5.1 layouts share.

The two version modules differ in their dtype names (5.1 adds the
``vtktype*`` spellings), their cell sections and their writers; what they read
alike lives here. A reader that maps a VTK type name takes the version's own
table as ``dtype_names``; the ``FIELD`` writer takes its inverse.
"""

from functools import reduce

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._mesh import Mesh
from .._vtk_common import Info, meshio_to_vtk_type

# supported vtk dataset types
vtk_dataset_types = [
    "UNSTRUCTURED_GRID",
    "STRUCTURED_POINTS",
    "STRUCTURED_GRID",
    "RECTILINEAR_GRID",
]
# additional infos per dataset type
vtk_dataset_infos = {
    "UNSTRUCTURED_GRID": [],
    "STRUCTURED_POINTS": [
        "DIMENSIONS",
        "ORIGIN",
        "SPACING",
        "ASPECT_RATIO",  # alternative for SPACING in version 1.0 and 2.0
    ],
    "STRUCTURED_GRID": ["DIMENSIONS"],
    "RECTILINEAR_GRID": [
        "DIMENSIONS",
        "X_COORDINATES",
        "Y_COORDINATES",
        "Z_COORDINATES",
    ],
}

# all main sections in vtk
vtk_sections = [
    "METADATA",
    "DATASET",
    "POINTS",
    "CELLS",
    "CELL_TYPES",
    "POINT_DATA",
    "CELL_DATA",
    "LOOKUP_TABLE",
    "COLOR_SCALARS",
]


def read_buffer(f, read_section, check_mesh, cells_from_info, dtype_names):
    """One legacy file after its version line.

    ``read_section(f, info)`` and ``check_mesh(info)`` are the version's own
    (their messages differ); ``cells_from_info(info)`` turns the parsed cell
    sections into ``(cells, cell_data)``."""
    # initialize output data
    info = Info()

    # skip title comment
    f.readline()

    data_type = f.readline().decode().strip().upper()
    if data_type not in ["ASCII", "BINARY"]:
        raise ReadError(f"Unknown VTK data type '{data_type}'.")
    info.is_ascii = data_type == "ASCII"

    while True:
        line = f.readline().decode()
        if not line:
            # EOF
            break

        line = line.strip()
        if len(line) == 0:
            continue

        info.split = line.split()
        info.section = info.split[0].upper()

        if info.section in vtk_sections:
            read_section(f, info)
        else:
            read_subsection(f, info, dtype_names)

    check_mesh(info)

    cells, cell_data = cells_from_info(info)

    return Mesh(
        info.points,
        cells,
        point_data=info.point_data,
        cell_data=cell_data,
        field_data=info.field_data,
    )


def read_subsection(f, info, dtype_names):
    if info.active == "POINT_DATA":
        d = info.point_data
    elif info.active == "CELL_DATA":
        d = info.cell_data_raw
    elif info.active == "DATASET":
        d = info.dataset
    else:
        d = info.field_data

    if info.section in vtk_dataset_infos[info.dataset["type"]]:
        if info.section[1:] == "_COORDINATES":
            info.num_points = int(info.split[1])
            data_type = info.split[2].lower()
            d[info.section] = _read_coords(
                f, data_type, info.is_ascii, info.num_points, dtype_names
            )
        else:
            if info.section == "DIMENSIONS":
                d[info.section] = list(map(int, info.split[1:]))
            else:
                d[info.section] = list(map(float, info.split[1:]))
            if len(d[info.section]) != 3:
                raise ReadError(
                    f"Wrong number of info in section '{info.section}'. "
                    f"Need 3, got {len(d[info.section])}."
                )
    elif info.section == "SCALARS":
        d.update(
            _read_scalar_field(
                f, info.num_items, info.split, info.is_ascii, dtype_names
            )
        )
    elif info.section == "VECTORS":
        d.update(
            _read_field(f, info.num_items, info.split, [3], info.is_ascii, dtype_names)
        )
    elif info.section == "TENSORS":
        d.update(
            _read_field(
                f, info.num_items, info.split, [3, 3], info.is_ascii, dtype_names
            )
        )
    elif info.section == "FIELD":
        d.update(_read_fields(f, int(info.split[2]), info.is_ascii, dtype_names))
    else:
        raise ReadError(f"Unknown section '{info.section}'.")


def write_field_data(f, data, binary, numpy_to_vtk_dtype):
    """One ``FIELD FieldData`` section; ``numpy_to_vtk_dtype`` is the version's
    numpy-name -> VTK type-name table."""
    f.write((f"FIELD FieldData {len(data)}\n").encode())
    for name, values in data.items():
        if isinstance(values, list):
            values = np.concatenate(values)
        if len(values.shape) == 1:
            num_tuples = values.shape[0]
            num_components = 1
        else:
            num_tuples = values.shape[0]
            num_components = values.shape[1]

        if " " in name:
            raise WriteError(f"VTK doesn't support spaces in field names ('{name}').")

        vtk_dtype = numpy_to_vtk_dtype[values.dtype.name]
        f.write(f"{name} {num_components} {num_tuples} {vtk_dtype}\n".encode())
        if binary:
            values.astype(values.dtype.newbyteorder(">")).tofile(f, sep="")
        else:
            # ascii
            values.tofile(f, sep=" ")
            # np.savetxt(f, points)
        f.write(b"\n")


def _generate_cells(dim):
    ele_dim = [d - 1 for d in dim if d > 1]
    # TODO use math.prod when requiring Python 3.8+? this would save the int conversion
    # <https://github.com/microsoft/pyright/issues/1226>
    ele_no = int(np.prod(ele_dim))
    spatial_dim = len(ele_dim)

    if spatial_dim == 0:
        return np.array([1, 0], dtype=np.int64), np.array([1], dtype=np.int64)
    if spatial_dim == 1:
        # cells are lines in 1D
        cells = np.empty((ele_no, 3), dtype=int)
        cells[:, 0] = 2
        cells[:, 1] = np.arange(ele_no, dtype=int)
        cells[:, 2] = cells[:, 1] + 1
        cell_types = np.full(ele_no, 3, dtype=int)
    elif spatial_dim == 2:
        # cells are quad in 2D
        cells = np.empty((ele_no, 5), dtype=int)
        cells[:, 0] = 4
        cells[:, 1] = np.arange(0, ele_no, dtype=int)
        cells[:, 1] += np.arange(0, ele_no, dtype=int) // ele_dim[0]
        cells[:, 2] = cells[:, 1] + 1
        cells[:, 3] = cells[:, 1] + 2 + ele_dim[0]
        cells[:, 4] = cells[:, 3] - 1
        cell_types = np.full(ele_no, 9, dtype=int)
    else:
        # cells are hex in 3D
        cells = np.empty((ele_no, 9), dtype=int)
        cells[:, 0] = 8
        cells[:, 1] = np.arange(ele_no)
        cells[:, 1] += (ele_dim[0] + ele_dim[1] + 1) * (
            np.arange(ele_no) // (ele_dim[0] * ele_dim[1])
        )
        cells[:, 1] += (np.arange(ele_no) % (ele_dim[0] * ele_dim[1])) // ele_dim[0]
        cells[:, 2] = cells[:, 1] + 1
        cells[:, 3] = cells[:, 1] + 2 + ele_dim[0]
        cells[:, 4] = cells[:, 3] - 1
        cells[:, 5] = cells[:, 1] + (1 + ele_dim[0]) * (1 + ele_dim[1])
        cells[:, 6] = cells[:, 5] + 1
        cells[:, 7] = cells[:, 5] + 2 + ele_dim[0]
        cells[:, 8] = cells[:, 7] - 1
        cell_types = np.full(ele_no, 12, dtype=int)

    return cells.reshape(-1), cell_types


def _generate_points(axis):
    x_dim = len(axis[0])
    y_dim = len(axis[1])
    z_dim = len(axis[2])
    pnt_no = x_dim * y_dim * z_dim
    x_id, y_id, z_id = np.mgrid[0:x_dim, 0:y_dim, 0:z_dim]
    points = np.empty((pnt_no, 3), dtype=axis[0].dtype)
    # VTK sorts points and cells in Fortran order
    points[:, 0] = axis[0][x_id.reshape(-1, order="F")]
    points[:, 1] = axis[1][y_id.reshape(-1, order="F")]
    points[:, 2] = axis[2][z_id.reshape(-1, order="F")]
    return points


def _read_coords(f, data_type, is_ascii, num_points, dtype_names):
    dtype = np.dtype(dtype_names[data_type])
    if is_ascii:
        coords = np.fromfile(f, count=num_points, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        coords = np.fromfile(f, count=num_points, dtype=dtype)
        line = f.readline().decode()
        if line != "\n":
            raise ReadError()
    return coords


def _read_points(f, data_type, is_ascii, num_points, dtype_names):
    dtype = np.dtype(dtype_names[data_type])
    if is_ascii:
        points = np.fromfile(f, count=num_points * 3, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        points = np.fromfile(f, count=num_points * 3, dtype=dtype)
        line = f.readline().decode()
        if line != "\n":
            raise ReadError()
    return points.reshape((num_points, 3))


def _read_cell_types(f, is_ascii, num_items):
    if is_ascii:
        ct = np.fromfile(f, count=int(num_items), sep=" ", dtype=int)
    else:
        # binary
        ct = np.fromfile(f, count=int(num_items), dtype=">i4")
        line = f.readline().decode()
        # Sometimes, there's no newline at the end
        if line.strip() != "":
            raise ReadError()
    return ct


def _read_scalar_field(f, num_data, split, is_ascii, dtype_names):
    data_name = split[1]
    data_type = split[2].lower()
    try:
        num_comp = int(split[3])
    except IndexError:
        num_comp = 1

    # The standard says:
    # > The parameter numComp must range between (1,4) inclusive; [...]
    if not (0 < num_comp < 5):
        raise ReadError("The parameter numComp must range between (1,4) inclusive")

    dtype = np.dtype(dtype_names[data_type])
    lt, _ = f.readline().decode().split()
    if lt.upper() != "LOOKUP_TABLE":
        raise ReadError()

    if is_ascii:
        data = np.fromfile(f, count=num_data * num_comp, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        data = np.fromfile(f, count=num_data * num_comp, dtype=dtype)
        line = f.readline().decode()
        if line != "\n":
            raise ReadError()

    data = data.reshape(-1, num_comp)
    return {data_name: data}


def _read_field(f, num_data, split, shape, is_ascii, dtype_names):
    data_name = split[1]
    data_type = split[2].lower()

    dtype = np.dtype(dtype_names[data_type])
    # prod()
    # <https://stackoverflow.com/q/2104782/353337>
    k = reduce((lambda x, y: x * y), shape)

    if is_ascii:
        data = np.fromfile(f, count=k * num_data, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        data = np.fromfile(f, count=k * num_data, dtype=dtype)
        line = f.readline().decode()
        if line != "\n":
            raise ReadError()

    data = data.reshape(-1, *shape)
    return {data_name: data}


def _read_fields(f, num_fields, is_ascii, dtype_names):
    data = {}
    for _ in range(num_fields):
        line = f.readline().decode().split()
        if line[0] == "METADATA":
            _skip_meta(f)
            name, shape0, shape1, data_type = f.readline().decode().split()
        else:
            name, shape0, shape1, data_type = line

        shape0 = int(shape0)
        shape1 = int(shape1)
        dtype = np.dtype(dtype_names[data_type.lower()])

        if is_ascii:
            dat = np.fromfile(f, count=shape0 * shape1, sep=" ", dtype=dtype)
        else:
            # Binary data is big endian, see
            # <https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
            dtype = dtype.newbyteorder(">")
            dat = np.fromfile(f, count=shape0 * shape1, dtype=dtype)
            line = f.readline().decode()
            if line != "\n":
                raise ReadError()

        if shape0 != 1:
            dat = dat.reshape((shape1, shape0))

        data[name] = dat

    return data


def _skip_meta(f):
    # skip possible metadata
    # https://vtk.org/doc/nightly/html/IOLegacyInformationFormat.html
    while True:
        line = f.readline().decode().strip()
        if not line:
            # end of metadata is a blank line
            break


# -- sections and writer pieces both versions share -------------------------------


def read_common_section(f, info, read_points, dataset_sep):
    """Read every section but ``CELLS``, whose layout differs by version.

    Returns ``False`` for a section left to the caller. ``read_points`` is the
    version's point reader (its dtype table bound); ``dataset_sep`` joins the
    supported dataset names in the unsupported-dataset message, which the two
    versions word apart (``"', '"`` in 4.2, ``", "`` in 5.1).
    """
    if info.section == "METADATA":
        _skip_meta(f)

    elif info.section == "DATASET":
        info.active = "DATASET"
        info.dataset["type"] = info.split[1].upper()
        if info.dataset["type"] not in vtk_dataset_types:
            legal = dataset_sep.join(vtk_dataset_types)
            raise ReadError(
                f"Only VTK '{legal}' supported (not {info.dataset['type']})."
            )

    elif info.section == "POINTS":
        info.active = "POINTS"
        info.num_points = int(info.split[1])
        data_type = info.split[2].lower()
        info.points = read_points(f, data_type, info.is_ascii, info.num_points)

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

    else:
        return False
    return True


def check_mesh(info):
    """Check an unstructured grid's sections, or build a structured one's cells.

    Structured datasets get points (where the file has none) and the counted
    legacy cell layout, ``(n, p0, p1, ...)`` per cell.
    """
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


def read_int_data(f, is_ascii, num_items, dtype=np.dtype("int32"), message=None):
    """``num_items`` integers, big-endian when binary (then a newline follows).

    A binary block not followed by its newline raises ``ReadError(message)``
    (a bare ``ReadError()`` when ``message`` is ``None``).
    """
    if is_ascii:
        c = np.fromfile(f, count=num_items, sep=" ", dtype=dtype)
    else:
        dtype = dtype.newbyteorder(">")
        c = np.fromfile(f, count=num_items, dtype=dtype)
        line = f.readline().decode()
        if line != "\n":
            raise ReadError() if message is None else ReadError(message)
    return c


def _pad_3d(array):
    return np.pad(array, ((0, 0), (0, 1)), "constant")


def pad_2d_to_3d(mesh, binary):
    """The points to write, with 2-D points and vectors given a zero third component.

    2-D ``point_data`` and ``cell_data`` vectors are padded in place on ``mesh``,
    each with a warning, and an ASCII write warns that it is for debugging.
    """
    if mesh.points.shape[1] == 2:
        warn(
            "VTK requires 3D points, but 2D points given. "
            "Appending 0 third component."
        )
        points = _pad_3d(mesh.points)
    else:
        points = mesh.points

    if mesh.point_data:
        for name, values in mesh.point_data.items():
            if len(values.shape) == 2 and values.shape[1] == 2:
                warn(
                    "VTK requires 3D vectors, but 2D vectors given. "
                    f"Appending 0 third component to {name}."
                )
                mesh.point_data[name] = _pad_3d(values)

    for name, data in mesh.cell_data.items():
        for k, values in enumerate(data):
            if len(values.shape) == 2 and values.shape[1] == 2:
                warn(
                    "VTK requires 3D vectors, but 2D vectors given. "
                    f"Appending 0 third component to {name}."
                )
                data[k] = _pad_3d(data[k])

    if not binary:
        warn("VTK ASCII files are only meant for debugging.")
    return points


def write_points(f, points, binary, numpy_to_vtk_dtype):
    """The ``POINTS`` section; binary data is big-endian, as the format requires."""
    dtype = numpy_to_vtk_dtype[points.dtype.name]
    f.write(f"POINTS {len(points)} {dtype}\n".encode())
    if binary:
        # <https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>
        points.astype(points.dtype.newbyteorder(">")).tofile(f, sep="")
    else:
        points.tofile(f, sep=" ")
    f.write(b"\n")


def write_cell_types(f, cells, binary):
    """The ``CELL_TYPES`` section, one VTK type per cell (big-endian int32)."""
    total_num_cells = sum(len(c.data) for c in cells)
    f.write(f"CELL_TYPES {total_num_cells}\n".encode())
    if binary:
        for c in cells:
            vtk_type = meshio_to_vtk_type[c.type]
            np.full(len(c.data), vtk_type, dtype=np.dtype(">i4")).tofile(f, sep="")
        f.write(b"\n")
    else:
        for c in cells:
            vtk_type = meshio_to_vtk_type[c.type]
            np.full(len(c.data), vtk_type).tofile(f, sep="\n")
            f.write(b"\n")


def write_unstructured(
    f, version, mesh, points, binary, write_cells, numpy_to_vtk_dtype
):
    """One ``UNSTRUCTURED_GRID`` file: header, points, cells and both data sections.

    ``write_cells(f, cells, binary)`` is the version's own ``CELLS`` writer.
    """
    f.write(f"# vtk DataFile Version {version}\n".encode())
    f.write(f"{_provenance.lines(_provenance.SlotTier.SINGLE_LINE)[0]}\n".encode())
    f.write(("BINARY\n" if binary else "ASCII\n").encode())
    f.write(b"DATASET UNSTRUCTURED_GRID\n")

    write_points(f, points, binary, numpy_to_vtk_dtype)
    write_cells(f, mesh.cells, binary)

    if mesh.point_data:
        num_points = mesh.points.shape[0]
        f.write(f"POINT_DATA {num_points}\n".encode())
        write_field_data(f, mesh.point_data, binary, numpy_to_vtk_dtype)

    if mesh.cell_data:
        total_num_cells = sum(len(c.data) for c in mesh.cells)
        f.write(f"CELL_DATA {total_num_cells}\n".encode())
        write_field_data(f, mesh.cell_data, binary, numpy_to_vtk_dtype)
