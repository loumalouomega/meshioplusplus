"""Readers the legacy VTK 4.2 and 5.1 layouts share.

The two version modules differ in their dtype names (5.1 adds the
``vtktype*`` spellings), their cell sections and their writers; what they read
alike lives here. A reader that maps a VTK type name takes the version's own
table as ``dtype_names``.
"""

from functools import reduce

import numpy as np

from .._exceptions import ReadError


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
