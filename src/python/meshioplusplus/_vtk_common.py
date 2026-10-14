import numpy as np

from ._common import num_nodes_per_cell, warn
from ._exceptions import ReadError
from ._mesh import CellBlock

# https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
vtk_to_meshio_type = {
    0: "empty",
    1: "vertex",
    # 2: 'poly_vertex',
    3: "line",
    # 4: 'poly_line',
    5: "triangle",
    # 6: 'triangle_strip',
    7: "polygon",
    8: "pixel",
    9: "quad",
    10: "tetra",
    # 11: 'voxel',
    12: "hexahedron",
    13: "wedge",
    14: "pyramid",
    15: "penta_prism",
    16: "hexa_prism",
    21: "line3",
    22: "triangle6",
    23: "quad8",
    24: "tetra10",
    25: "hexahedron20",
    26: "wedge15",
    27: "pyramid13",
    28: "quad9",
    29: "hexahedron27",
    30: "quad6",
    31: "wedge12",
    32: "wedge18",
    33: "hexahedron24",
    34: "triangle7",
    35: "line4",
    42: "polyhedron",
    #
    # 60: VTK_HIGHER_ORDER_EDGE,
    # 61: VTK_HIGHER_ORDER_TRIANGLE,
    # 62: VTK_HIGHER_ORDER_QUAD,
    # 63: VTK_HIGHER_ORDER_POLYGON,
    # 64: VTK_HIGHER_ORDER_TETRAHEDRON,
    # 65: VTK_HIGHER_ORDER_WEDGE,
    # 66: VTK_HIGHER_ORDER_PYRAMID,
    # 67: VTK_HIGHER_ORDER_HEXAHEDRON,
    # Arbitrary order Lagrange elements
    68: "VTK_LAGRANGE_CURVE",
    69: "VTK_LAGRANGE_TRIANGLE",
    70: "VTK_LAGRANGE_QUADRILATERAL",
    71: "VTK_LAGRANGE_TETRAHEDRON",
    72: "VTK_LAGRANGE_HEXAHEDRON",
    73: "VTK_LAGRANGE_WEDGE",
    74: "VTK_LAGRANGE_PYRAMID",
    # Arbitrary order Bezier elements
    75: "VTK_BEZIER_CURVE",
    76: "VTK_BEZIER_TRIANGLE",
    77: "VTK_BEZIER_QUADRILATERAL",
    78: "VTK_BEZIER_TETRAHEDRON",
    79: "VTK_BEZIER_HEXAHEDRON",
    80: "VTK_BEZIER_WEDGE",
    81: "VTK_BEZIER_PYRAMID",
}
meshio_to_vtk_type = {v: k for k, v in vtk_to_meshio_type.items()}


def vtk_to_meshio_order(vtk_type, dtype=int):
    # meshio uses the same node ordering as VTK for most cell types. However, for the
    # linear wedge, the ordering of the gmsh Prism [1] is adopted since this is found in
    # most codes (Abaqus, Ansys, Nastran,...). In the vtkWedge [2], the normal of the
    # (0,1,2) triangle points outwards, while in gmsh this normal points inwards.
    # [1] http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    # [2] https://vtk.org/doc/nightly/html/classvtkWedge.html
    if vtk_type == 13:
        return np.array([0, 2, 1, 3, 5, 4], dtype=dtype)
    return None


def meshio_to_vtk_order(meshio_type, dtype=int):
    if meshio_type == "wedge":
        return np.array([0, 2, 1, 3, 5, 4], dtype=dtype)
    return None


def vtk_cells_from_data(connectivity, offsets, types, cell_data_raw):
    # Translate it into the cells array.
    # `connectivity` is a one-dimensional vector with
    # (p00, p01, ... ,p0k, p10, p11, ..., p1k, ...
    # `offsets` is a pointer array that points to the first position of p0, p1, etc.
    if len(offsets) != len(types):
        raise ReadError(f"len(offsets) != len(types) ({len(offsets)} != {len(types)})")

    # A grid with no cells (an empty piece of a partitioned file) has no blocks.
    if len(types) == 0:
        return [], {}

    # identify cell blocks
    breaks = np.where(types[:-1] != types[1:])[0] + 1
    # all cells with indices between start[k] and end[k] have the same type
    start_end = list(
        zip(
            np.concatenate([[0], breaks]),
            np.concatenate([breaks, [len(types)]]),
        )
    )

    cells = []
    cell_data = {}

    for start, end in start_end:
        try:
            meshio_type = vtk_to_meshio_type[types[start]]
        except KeyError:
            warn(
                f"File contains cells that meshio cannot handle (type {types[start]})."
            )
            continue

        # cells with varying number of points
        special_cells = [
            "polygon",
            "VTK_LAGRANGE_CURVE",
            "VTK_LAGRANGE_TRIANGLE",
            "VTK_LAGRANGE_QUADRILATERAL",
            "VTK_LAGRANGE_TETRAHEDRON",
            "VTK_LAGRANGE_HEXAHEDRON",
            "VTK_LAGRANGE_WEDGE",
            "VTK_LAGRANGE_PYRAMID",
            "VTK_BEZIER_CURVE",
            "VTK_BEZIER_TRIANGLE",
            "VTK_BEZIER_QUADRILATERAL",
            "VTK_BEZIER_TETRAHEDRON",
            "VTK_BEZIER_HEXAHEDRON",
            "VTK_BEZIER_WEDGE",
            "VTK_BEZIER_PYRAMID",
        ]
        if meshio_type in special_cells:
            # Polygons have unknown and varying number of nodes per cell.

            # Index where the previous block of cells stopped. Needed to know the number
            # of nodes for the first cell in the block.
            first_node = 0 if start == 0 else offsets[start - 1]

            # Start off the cell-node relation for each cell in this block
            start_cn = np.hstack((first_node, offsets[start:end]))
            # Find the size of each cell except the last
            sizes = np.diff(start_cn)

            # find where the cell blocks start and end
            b = np.diff(sizes)
            c = np.concatenate([[0], np.where(b != 0)[0] + 1, [len(sizes)]])

            # Loop over all cell sizes, find all cells with this size, and assign
            # connectivity
            for cell_block_start, cell_block_end in zip(c, c[1:]):
                items = np.arange(cell_block_start, cell_block_end)
                sz = sizes[cell_block_start]

                new_order = vtk_to_meshio_order(types[start], dtype=offsets.dtype)
                if new_order is None:
                    new_order = np.arange(sz, dtype=offsets.dtype)
                new_order -= sz

                indices = np.add.outer(start_cn[items + 1], new_order)
                cells.append(CellBlock(meshio_type, connectivity[indices]))

                # Store cell data for this set of cells
                for name, d in cell_data_raw.items():
                    if name not in cell_data:
                        cell_data[name] = []
                    cell_data[name].append(d[start + items])
        else:
            # Non-polygonal cell. Same number of nodes per cell makes everything easier.
            n = num_nodes_per_cell[meshio_type]

            new_order = vtk_to_meshio_order(types[start], dtype=offsets.dtype)
            if new_order is None:
                new_order = np.arange(n, dtype=offsets.dtype)
            new_order -= n

            indices = np.add.outer(offsets[start:end], new_order)
            cells.append(CellBlock(meshio_type, connectivity[indices]))
            for name, d in cell_data_raw.items():
                if name not in cell_data:
                    cell_data[name] = []
                cell_data[name].append(d[start:end])

    return cells, cell_data


class Info:
    """Info Container for the VTK reader."""

    def __init__(self):
        self.points = None
        self.field_data = {}
        self.cell_data_raw = {}
        self.point_data = {}
        self.dataset = {}
        self.connectivity = None
        self.offsets = None
        self.types = None
        self.active = None
        self.is_ascii = False
        self.split = []
        self.num_items = 0
        # One of the problem in reading VTK files are POINT_DATA and CELL_DATA fields.
        # They can contain a number of SCALARS+LOOKUP_TABLE tables, without giving and
        # indication of how many there are. Hence, SCALARS must be treated like a
        # first-class section.  To associate it with POINT/CELL_DATA, we store the
        # `active` section in this variable.
        self.section = None


# VTK 9.1 swapped the nodes of edges 10 and 11 of its higher-order hexahedra
# (both run along the third direction) and marks files in the new order by
# version: XML 2.2 and later, legacy 5.1. meshio++ holds the new order, so a
# reader converts an older file and a writer of an older version converts on
# the way out (the twin of the core's detail/vtk_hex_order.hpp). The swap is
# its own inverse.
HEX_ORDER_TYPES = (72, 79)  # VTK_LAGRANGE_HEXAHEDRON, VTK_BEZIER_HEXAHEDRON


def hex_swap_edges(nodes, degrees=None):
    """Swap edges 10 and 11 of one higher-order hexahedron in place; False when
    the degrees (given, or from the node count) do not explain its size."""
    n = len(nodes)
    if degrees is None:
        q = int(round(n ** (1.0 / 3.0)))
        p = (q - 1, q - 1, q - 1)
    else:
        p = tuple(int(d) for d in degrees)
    if min(p) < 1 or (p[0] + 1) * (p[1] + 1) * (p[2] + 1) != n:
        return False
    m = p[2] - 1
    s = 8 + 4 * (p[0] - 1) + 4 * (p[1] - 1) + 2 * m
    a = nodes[s : s + m].copy()
    nodes[s : s + m] = nodes[s + m : s + 2 * m]
    nodes[s + m : s + 2 * m] = a
    return True


def hex_swap_stream(conn, offsets, types, degrees, what):
    """Swap edges 10 and 11 of every higher-order hexahedron of a VTK cell
    stream (``offsets`` are end offsets into ``conn``, changed in place);
    ``degrees`` is the ``HigherOrderDegrees`` cell array or None."""
    types = np.asarray(types).reshape(-1)
    if degrees is not None:
        degrees = np.asarray(degrees).reshape(-1)
        degrees = degrees.reshape(-1, 3) if degrees.size == 3 * len(types) else None
    warned = False
    for c in np.flatnonzero(np.isin(types, HEX_ORDER_TYPES)):
        start = 0 if c == 0 else int(offsets[c - 1])
        end = int(offsets[c])
        if end > start and not hex_swap_edges(
            conn[start:end], None if degrees is None else degrees[c]
        ):
            if not warned:
                warn(
                    f"{what}: a higher-order hexahedron of {end - start} nodes in a file "
                    "older than VTK 9.1 has no degrees that explain its size; its edge "
                    "order is left as written"
                )
                warned = True
