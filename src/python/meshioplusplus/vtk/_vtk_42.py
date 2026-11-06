"""
I/O for VTK <https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf>.
"""

import functools

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError
from .._files import open_file
from .._mesh import CellBlock
from .._vtk_common import (
    HEX_ORDER_TYPES,
    hex_swap_edges,
    meshio_to_vtk_order,
    meshio_to_vtk_type,
    vtk_to_meshio_order,
    vtk_to_meshio_type,
)
from . import _legacy_common as _legacy
from ._legacy_common import (
    _generate_cells,
    _generate_points,
    _read_cell_types,
    _skip_meta,
    vtk_dataset_types,
)

vtk_type_to_numnodes = np.array(
    [
        0,  # empty
        1,  # vertex
        -1,  # poly_vertex
        2,  # line
        -1,  # poly_line
        3,  # triangle
        -1,  # triangle_strip
        -1,  # polygon
        -1,  # pixel
        4,  # quad
        4,  # tetra
        -1,  # voxel
        8,  # hexahedron
        6,  # wedge
        5,  # pyramid
        10,  # penta_prism
        12,  # hexa_prism
        -1,
        -1,
        -1,
        -1,
        3,  # line3
        6,  # triangle6
        8,  # quad8
        10,  # tetra10
        20,  # hexahedron20
        15,  # wedge15
        13,  # pyramid13
        9,  # quad9
        27,  # hexahedron27
        6,  # quad6
        12,  # wedge12
        18,  # wedge18
        24,  # hexahedron24
        7,  # triangle7
        4,  # line4
    ]
)


# These are all VTK data types.
vtk_to_numpy_dtype_name = {
    "bit": "bool",
    "unsigned_char": "uint8",
    "char": "int8",
    "unsigned_short": "uint16",
    "short": "int16",
    "unsigned_int": "uint32",
    "int": "int32",
    "unsigned_long": "uint64",
    "long": "int64",
    "float": "float32",
    "double": "float64",
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
    return translate_cells(info.connectivity, info.types, info.cell_data_raw)


def _read_section(f, info):
    if info.section == "METADATA":
        _skip_meta(f)

    elif info.section == "DATASET":
        info.active = "DATASET"
        info.dataset["type"] = info.split[1].upper()
        if info.dataset["type"] not in vtk_dataset_types:
            raise ReadError(
                "Only VTK '{}' supported (not {}).".format(
                    "', '".join(vtk_dataset_types), info.dataset["type"]
                )
            )

    elif info.section == "POINTS":
        info.active = "POINTS"
        info.num_points = int(info.split[1])
        data_type = info.split[2].lower()
        info.points = _read_points(f, data_type, info.is_ascii, info.num_points)

    elif info.section == "CELLS":
        info.active = "CELLS"
        info.num_items = int(info.split[2])
        info.connectivity = _read_int_data(f, info.is_ascii, info.num_items)

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


def _read_int_data(f, is_ascii, num_items, dtype=np.dtype("int32")):
    if is_ascii:
        c = np.fromfile(f, count=num_items, sep=" ", dtype=dtype)
    else:
        dtype = dtype.newbyteorder(">")
        c = np.fromfile(f, count=num_items, dtype=dtype)
        line = f.readline().decode()
        if line != "\n":
            raise ReadError()
    return c


def translate_cells(connectivity, types, cell_data_raw):
    # https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
    # Translate it into the cells array.
    # `data` is a one-dimensional vector with
    # (num_points0, p0, p1, ... ,pk, numpoints1, p10, p11, ..., p1k, ...
    # or a tuple with (offsets, connectivity)
    has_polygon = np.any(types == meshio_to_vtk_type["polygon"])

    cells = []
    cell_data = {}
    if has_polygon:
        numnodes = np.empty(len(types), dtype=int)
        # If some polygons are in the VTK file, loop over the cells
        numcells = len(types)
        offsets = np.empty(len(types), dtype=int)
        offsets[0] = 0
        for idx in range(numcells - 1):
            numnodes[idx] = connectivity[offsets[idx]]
            offsets[idx + 1] = offsets[idx] + numnodes[idx] + 1

        idx = numcells - 1
        numnodes[idx] = connectivity[offsets[idx]]
        if not np.all(numnodes == connectivity[offsets]):
            raise ReadError()

        # TODO: cell_data
        for idx, vtk_cell_type in enumerate(types):
            start = offsets[idx] + 1

            new_order = vtk_to_meshio_order(vtk_cell_type, offsets.dtype)
            if new_order is None:
                new_order = np.arange(numnodes[idx], dtype=offsets.dtype)

            cell_idx = start + new_order

            cell = connectivity[cell_idx]

            cell_type = vtk_to_meshio_type[vtk_cell_type]

            if (
                len(cells) > 0
                and cells[-1][0] == cell_type
                # the following check if needed for polygons; key can be equal, but
                # still needs to go into a new block
                and len(cell) == len(cells[-1][1][-1])
            ):
                cells[-1][1].append(cell)
            else:
                # open up a new cell block
                cells.append((cell_type, [cell]))
    else:
        # Infer offsets from the cell types. This is much faster than manually going
        # through the data array. Slight disadvantage: This doesn't work for cells with
        # a custom number of points.
        if np.any(types >= len(vtk_type_to_numnodes)):
            illegal_types = ", ".join(
                str(item) for item in set(types[types >= len(vtk_type_to_numnodes)])
            )
            raise ReadError(
                f"File contains cells that meshio cannot handle (types {illegal_types})"
            )

        numnodes = vtk_type_to_numnodes[types]
        if not np.all(numnodes > 0):
            raise ReadError("File contains cells that meshio cannot handle.")

        offsets = np.cumsum(numnodes + 1) - (numnodes + 1)

        if not np.all(numnodes == connectivity[offsets]):
            raise ReadError()
        idx0 = 1

        b = np.concatenate(
            [[0], np.where(types[:-1] != types[1:])[0] + 1, [len(types)]]
        )
        for start, end in zip(b[:-1], b[1:]):
            if start == end:
                continue
            meshio_type = vtk_to_meshio_type[types[start]]
            n = numnodes[start]
            new_order = vtk_to_meshio_order(types[start], dtype=offsets.dtype)
            if new_order is None:
                cell_idx = idx0 + np.arange(0, n, dtype=offsets.dtype)
            else:
                cell_idx = idx0 + new_order
            indices = np.add.outer(offsets[start:end], cell_idx)
            cells.append((meshio_type, connectivity[indices]))
            for name, d in cell_data_raw.items():
                if name not in cell_data:
                    cell_data[name] = []
                cell_data[name].append(d[start:end])

    return cells, cell_data


def write(filename, mesh, binary=True):
    def pad(array):
        return np.pad(array, ((0, 0), (0, 1)), "constant")

    if mesh.points.shape[1] == 2:
        warn(
            "VTK requires 3D points, but 2D points given. "
            + "Appending 0 third component."
        )
        points = pad(mesh.points)
    else:
        points = mesh.points

    if mesh.point_data:
        for name, values in mesh.point_data.items():
            if len(values.shape) == 2 and values.shape[1] == 2:
                warn(
                    "VTK requires 3D vectors, but 2D vectors given. "
                    + f"Appending 0 third component to {name}."
                )
                mesh.point_data[name] = pad(values)

    for name, data in mesh.cell_data.items():
        for k, values in enumerate(data):
            if len(values.shape) == 2 and values.shape[1] == 2:
                warn(
                    "VTK requires 3D vectors, but 2D vectors given. "
                    + f"Appending 0 third component to {name}."
                )
                data[k] = pad(data[k])

    if not binary:
        warn("VTK ASCII files are only meant for debugging.")

    with open_file(filename, "wb") as f:
        f.write(b"# vtk DataFile Version 4.2\n")
        f.write(f"{_provenance.lines(_provenance.SlotTier.SINGLE_LINE)[0]}\n".encode())
        f.write(("BINARY\n" if binary else "ASCII\n").encode())
        f.write(b"DATASET UNSTRUCTURED_GRID\n")

        # write points and cells
        _write_points(f, points, binary)
        _write_cells(f, mesh.cells, binary, mesh.cell_data.get("HigherOrderDegrees"))

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
    f.write(
        "POINTS {} {}\n".format(
            len(points), numpy_to_vtk_dtype[points.dtype.name]
        ).encode()
    )

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


def _old_hex_order(cells, degrees):
    """The blocks with their higher-order hexahedra in the order VTK reads
    from a 4.2 file: edges 10 and 11 swapped (``hex_swap_edges``)."""
    out, warned = [], False
    for k, c in enumerate(cells):
        if meshio_to_vtk_type.get(c.type) not in HEX_ORDER_TYPES:
            out.append(c)
            continue
        data = np.array(c.data, dtype=np.int64)
        deg = None
        if degrees is not None and len(degrees) == len(cells):
            deg = np.asarray(degrees[k]).reshape(-1)
            deg = deg.reshape(-1, 3) if deg.size == 3 * len(data) else None
        for r in range(len(data)):
            if (
                not hex_swap_edges(data[r], None if deg is None else deg[r])
                and not warned
            ):
                warn(
                    f"VTK 4.2 writer: a {c.type} of {data.shape[1]} nodes has no degrees "
                    "that explain its size; VTK will read its edges 10 and 11 swapped "
                    "(write 5.1)"
                )
                warned = True
        out.append(CellBlock(c.type, data))
    return out


def _write_cells(f, cells, binary, degrees=None):
    # VTK reads the higher-order hexahedra of a file older than 5.1 with edges
    # 10 and 11 swapped, so they are written so.
    cells = _old_hex_order(cells, degrees)
    total_num_cells = sum(len(c.data) for c in cells)
    total_num_idx = sum(c.data.size for c in cells)
    # For each cell, the number of nodes is stored
    total_num_idx += total_num_cells
    f.write(f"CELLS {total_num_cells} {total_num_idx}\n".encode())
    if binary:
        for c in cells:
            n = c.data.shape[1]
            d = c.data

            new_order = meshio_to_vtk_order(c.type)
            if new_order is not None:
                d = d[:, new_order]

            dtype = np.dtype(">i4")
            # One must force endianness here:
            # <https://github.com/numpy/numpy/issues/15088>
            np.column_stack(
                [np.full(d.shape[0], n, dtype=dtype), d.astype(dtype)],
            ).astype(dtype).tofile(f, sep="")
        f.write(b"\n")
    else:
        # ascii
        for c in cells:
            n = c.data.shape[1]

            d = c.data
            new_order = meshio_to_vtk_order(c.type)
            if new_order is not None:
                d = d[:, new_order]

            # prepend a column with the value n
            np.column_stack(
                [np.full(c.data.shape[0], n, dtype=c.data.dtype), d]
            ).tofile(f, sep="\n")
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
