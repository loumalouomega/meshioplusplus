"""Node-ordering permutations between file formats and meshio++'s (VTK) cell
layouts, keyed by ``(format, cell type)``. Twin of ``detail/node_order.cpp``.

Every entry holds both directions as *gather* tables:

- ``to_meshio[k]``: the file slot meshio++ node ``k`` comes from on read,
  ``meshio[k] = file[to_meshio[k]]``;
- ``from_meshio[j]``: the meshio++ node file slot ``j`` comes from on write,
  ``file[j] = meshio[from_meshio[j]]``.

A type with no entry uses the identity. Where each table was pinned is listed
in ``doc/node_ordering.md``.
"""

from __future__ import annotations

from typing import NamedTuple

__all__ = ["NodeOrder", "node_order", "node_order_keys", "to_meshio", "from_meshio"]


class NodeOrder(NamedTuple):
    to_meshio: tuple
    from_meshio: tuple


def _inverse(table):
    inv = [0] * len(table)
    for i, j in enumerate(table):
        inv[j] = i
    return tuple(inv)


_TO, _FROM = "to_meshio", "from_meshio"

# Each table is written in the direction its source gives it (the comments in
# detail/node_order.cpp say where each one comes from); the other is derived.
_SOURCES = [
    ("med", "tetra", _TO, [0, 1, 3, 2]),
    ("med", "pyramid", _TO, [0, 3, 2, 1, 4]),
    ("med", "wedge", _TO, [3, 4, 5, 0, 1, 2]),
    ("med", "hexahedron", _TO, [4, 5, 6, 7, 0, 1, 2, 3]),
    ("med", "tetra10", _TO, [0, 1, 3, 2, 4, 8, 7, 6, 5, 9]),
    ("med", "pyramid13", _TO, [0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10]),
    ("med", "wedge15", _TO, [3, 4, 5, 0, 1, 2, 9, 10, 11, 6, 7, 8, 12, 13, 14]),
    (
        "med",
        "wedge18",
        _TO,
        [3, 4, 5, 0, 1, 2, 9, 10, 11, 6, 7, 8] + list(range(12, 18)),
    ),
    (
        "med",
        "hexahedron20",
        _TO,
        [4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11] + [16, 17, 18, 19],
    ),
    (
        "med",
        "hexahedron27",
        _TO,
        [4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11]
        + [16, 17, 18, 19, 24, 22, 21, 23, 25, 20, 26],
    ),
    ("code_aster", "wedge15", _TO, list(range(9)) + [12, 13, 14, 9, 10, 11]),
    (
        "code_aster",
        "wedge18",
        _TO,
        list(range(9)) + [12, 13, 14, 9, 10, 11, 15, 16, 17],
    ),
    (
        "code_aster",
        "hexahedron20",
        _TO,
        list(range(12)) + [16, 17, 18, 19, 12, 13, 14, 15],
    ),
    (
        "code_aster",
        "hexahedron27",
        _TO,
        list(range(12))
        + [16, 17, 18, 19, 12, 13, 14, 15]
        + [24, 22, 21, 23, 20, 25, 26],
    ),
    ("frd", "hexahedron20", _TO, list(range(12)) + [16, 17, 18, 19, 12, 13, 14, 15]),
    ("frd", "wedge15", _TO, list(range(9)) + [12, 13, 14, 9, 10, 11]),
    ("frd", "line3", _TO, [0, 2, 1]),
    # MSC Patran neutral file: hex20/wedge15 put the vertical mid-edges before
    # the top ring (Patran Reference Manual, Element Library).
    ("patran", "hexahedron20", _TO, list(range(12)) + [16, 17, 18, 19, 12, 13, 14, 15]),
    ("patran", "wedge15", _TO, list(range(9)) + [12, 13, 14, 9, 10, 11]),
    # FEBio .feb/.xplt: hex27's mid-height face centres run y-, x+, y+, x-
    # (FECore FEHex27); every other FEBio type is in meshio++'s order.
    ("febio", "hexahedron27", _TO, list(range(20)) + [23, 21, 20, 22, 24, 25, 26]),
    # Elmer mesh directory: 820/827 put the vertical mid-edge nodes before the top
    # ring, and the 827 mid-height face centres run y-, x+, y+, x- (elements.def;
    # ElmerSolver's VTU writer applies the same permutation).
    ("elmer", "hexahedron20", _TO, list(range(12)) + [16, 17, 18, 19, 12, 13, 14, 15]),
    (
        "elmer",
        "hexahedron27",
        _TO,
        list(range(12)) + [16, 17, 18, 19, 12, 13, 14, 15, 23, 21, 20, 22, 24, 25, 26],
    ),
    # COMSOL .mphtxt/.mphbin: corners in tensor order (x fastest), then every
    # other node of the quadratic lattice in lexicographic (z, y, x) order;
    # checked against real COMSOL files and AWS Palace's COMSOL->gmsh tables.
    ("mphtxt", "quad", _TO, [0, 1, 3, 2]),
    ("mphtxt", "hexahedron", _TO, [0, 1, 3, 2, 4, 5, 7, 6]),
    ("mphtxt", "pyramid", _TO, [0, 1, 3, 2, 4]),
    ("mphtxt", "triangle6", _TO, [0, 1, 2, 3, 5, 4]),
    ("mphtxt", "quad9", _TO, [0, 1, 3, 2, 4, 7, 8, 5, 6]),
    ("mphtxt", "tetra10", _TO, [0, 1, 2, 3, 4, 6, 5, 7, 8, 9]),
    (
        "mphtxt",
        "hexahedron27",
        _TO,
        [
            0,
            1,
            3,
            2,
            4,
            5,
            7,
            6,
            8,
            11,
            12,
            9,
            22,
            25,
            26,
            23,
            13,
            15,
            21,
            19,
            16,
            18,
            14,
            20,
            10,
            24,
            17,
        ],
    ),
    (
        "mphtxt",
        "wedge18",
        _TO,
        [0, 1, 2, 3, 4, 5, 6, 8, 7, 15, 17, 16, 9, 11, 14, 10, 13, 12],
    ),
    ("mphtxt", "pyramid14", _TO, [0, 1, 3, 2, 4, 5, 8, 9, 6, 10, 11, 13, 12, 7]),
    ("unv", "line3", _FROM, [0, 2, 1]),
    ("unv", "triangle6", _FROM, [0, 3, 1, 4, 2, 5]),
    ("unv", "quad8", _FROM, [0, 4, 1, 5, 2, 6, 3, 7]),
    ("unv", "quad9", _FROM, [0, 4, 1, 5, 2, 6, 3, 7, 8]),
    ("unv", "tetra10", _FROM, [0, 4, 1, 5, 2, 6, 7, 8, 9, 3]),
    ("unv", "pyramid13", _FROM, [0, 5, 1, 6, 2, 7, 3, 8, 9, 10, 11, 12, 4]),
    ("unv", "wedge15", _FROM, [0, 6, 1, 7, 2, 8, 12, 13, 14, 3, 9, 4, 10, 5, 11]),
    (
        "unv",
        "hexahedron20",
        _FROM,
        [0, 8, 1, 9, 2, 10, 3, 11, 16, 17, 18, 19] + [4, 12, 5, 13, 6, 14, 7, 15],
    ),
]


def _build():
    out = {}
    for fmt, cell_type, direction, table in _SOURCES:
        table = tuple(table)
        if direction == _TO:
            out[(fmt, cell_type)] = NodeOrder(table, _inverse(table))
        else:
            out[(fmt, cell_type)] = NodeOrder(_inverse(table), table)
    return out


_ORDERS = _build()


def node_order(fmt, cell_type):
    """The :class:`NodeOrder` for ``(fmt, cell_type)``, or ``None`` for the
    identity."""
    return _ORDERS.get((fmt, cell_type))


def node_order_keys():
    """Every ``(format, cell type)`` key, sorted."""
    return sorted(_ORDERS)


def to_meshio(fmt, cell_type, data):
    """Reorder a ``(n, k)`` connectivity array from the file's node order to
    meshio++'s; unchanged when the table is the identity or of another width."""
    order = _ORDERS.get((fmt, cell_type))
    if order is None or data.shape[1] != len(order.to_meshio):
        return data
    return data[:, list(order.to_meshio)]


def from_meshio(fmt, cell_type, data):
    """Reorder a ``(n, k)`` connectivity array from meshio++'s node order to the
    file's; the inverse of :func:`to_meshio`."""
    order = _ORDERS.get((fmt, cell_type))
    if order is None or data.shape[1] != len(order.from_meshio):
        return data
    return data[:, list(order.from_meshio)]
