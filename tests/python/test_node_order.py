"""The node-ordering registry (``_node_order.py``, twin of ``detail/node_order.cpp``).

Every table is a permutation whose two directions are inverses. The tables of
gmsh, CGNS, GiD, Kratos and Exodus are pinned geometrically: an element built
in the file's order from that format's own documentation (mid-edge nodes on
edge midpoints, face and body centres on centroids) must read as a valid
meshio++ element. For Code_Aster the tables are pinned against Code_Aster's own
readers:

- ``bibfor/prepost/inigms.F90`` (gmsh -> Aster, used by ``PRE_GMSH``) composed
  with meshio++'s gmsh tables must give the ``code_aster`` table exactly;
- ``bibfor/third_party_interf/lrmtyp.F90`` (MED <-> Aster, used by
  ``LIRE_MAILLAGE(FORMAT='MED')``) composed with the ``med`` tables must give
  the same element, up to a symmetry of the reference cell.

The two Code_Aster tables below are transcribed, 0-based, from those sources.
"""

import numpy as np
import pytest

from meshioplusplus import _node_order
from meshioplusplus.gmsh.common import _gmsh_to_meshio_order

# ---------------------------------------------------------------------------
# meshio++ (VTK) reference elements: corners, then the mid-edge nodes in the
# order of ``edges``, then the centres of ``faces``, then the body centre.
# ---------------------------------------------------------------------------

_TET = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
_TET_E = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
# fmt: off
_HEX = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
_HEX_E = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
          (0, 4), (1, 5), (2, 6), (3, 7)]
_HEX_F = [(0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4), (3, 7, 6, 2), (0, 3, 2, 1), (4, 5, 6, 7)]
# fmt: on
_WED = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]]
_WED_E = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]
_WED_F = [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
_PYR = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 1]]
_PYR_E = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)]
_TRI = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
_TRI_E = [(0, 1), (1, 2), (2, 0)]
_QUA = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
_QUA_E = [(0, 1), (1, 2), (2, 3), (3, 0)]

# cell type -> (corners, edges, faces, has body centre)
_SHAPES = {
    "line": ([[0, 0, 0], [1, 0, 0]], [], [], False),
    "line3": ([[0, 0, 0], [1, 0, 0]], [(0, 1)], [], False),
    "triangle6": (_TRI, _TRI_E, [], False),
    "quad8": (_QUA, _QUA_E, [], False),
    "quad9": (_QUA, _QUA_E, [(0, 1, 2, 3)], False),
    "tetra": (_TET, [], [], False),
    "tetra10": (_TET, _TET_E, [], False),
    "pyramid": (_PYR, [], [], False),
    "pyramid13": (_PYR, _PYR_E, [], False),
    "wedge": (_WED, [], [], False),
    "wedge15": (_WED, _WED_E, [], False),
    "wedge18": (_WED, _WED_E, _WED_F, False),
    "hexahedron": (_HEX, [], [], False),
    "hexahedron20": (_HEX, _HEX_E, [], False),
    "hexahedron27": (_HEX, _HEX_E, _HEX_F, True),
    "quad": (_QUA, [], [], False),
    "triangle": (_TRI, [], [], False),
    "pyramid14": (_PYR, _PYR_E, [(0, 1, 2, 3)], False),
}


def _reference(cell_type):
    corners, edges, faces, body = _SHAPES[cell_type]
    c = np.array(corners, dtype=float)
    parts = [c]
    parts += [(c[a] + c[b])[None] / 2 for a, b in edges]
    parts += [c[list(f)].mean(0)[None] for f in faces]
    if body:
        parts.append(c.mean(0)[None])
    return np.vstack(parts)


def _is_valid(x, cell_type):
    """Mid-edge nodes on midpoints, centres on centroids, and a positive volume
    (a positive area normal along +z for 2-D cells) for ``x`` in meshio order."""
    corners, edges, faces, body = _SHAPES[cell_type]
    nc = len(corners)
    c = x[:nc]
    k = nc
    for a, b in edges:
        if not np.allclose(x[k], (c[a] + c[b]) / 2):
            return False
        k += 1
    for f in faces:
        if not np.allclose(x[k], c[list(f)].mean(0)):
            return False
        k += 1
    if body and not np.allclose(x[k], c.mean(0)):
        return False
    if nc == 2:
        return True
    if nc in (3, 4) and cell_type not in ("tetra", "tetra10"):
        return np.cross(c[1] - c[0], c[2] - c[0])[2] > 0
    tet = {4: (0, 1, 2, 3), 5: (0, 1, 3, 4), 6: (0, 1, 2, 3), 8: (0, 1, 3, 4)}[nc]
    m = np.array([c[tet[i]] - c[tet[0]] for i in (1, 2, 3)])
    return np.linalg.det(m) > 0


_ALL = _node_order.node_order_keys()


@pytest.mark.parametrize("fmt, cell_type", _ALL, ids=[f"{f}-{t}" for f, t in _ALL])
def test_every_table_is_a_permutation_and_its_inverse(fmt, cell_type):
    order = _node_order.node_order(fmt, cell_type)
    n = len(_reference(cell_type))
    assert sorted(order.to_meshio) == list(range(n))
    assert [order.to_meshio[j] for j in order.from_meshio] == list(range(n))


# ---------------------------------------------------------------------------
# File-order reference elements, transcribed from each format's own
# documentation rather than from the tables: every file slot past the corners
# is the centroid of the listed corners (meshio++ corner numbering; the corners
# of these formats are meshio++'s). Built in file order and read through
# `to_meshio`, each must come out a valid meshio++ element.
# ---------------------------------------------------------------------------

_HEX_CORNERS_ONLY = [(k,) for k in range(8)]
_WED_CORNERS_ONLY = [(k,) for k in range(6)]
_PYR_CORNERS_ONLY = [(k,) for k in range(5)]
# Bottom ring, verticals, top ring: CGNS SIDS HEXA_20, SEACAS Ioss Hex20 and
# Kratos's Hexahedra3D20 (PointsLocalCoordinates, GenerateEdges).
_HEX20_VERTICALS_FIRST = _HEX_CORNERS_ONLY + [
    (0, 1), (1, 2), (2, 3), (3, 0),
    (0, 4), (1, 5), (2, 6), (3, 7),
    (4, 5), (5, 6), (6, 7), (7, 4),
]  # fmt: skip
_WED15_VERTICALS_FIRST = _WED_CORNERS_ONLY + [
    (0, 1), (1, 2), (2, 0),
    (0, 3), (1, 4), (2, 5),
    (3, 4), (4, 5), (5, 3),
]  # fmt: skip
_Z_MINUS, _Z_PLUS = (0, 1, 2, 3), (4, 5, 6, 7)
_X_MINUS, _X_PLUS = (0, 3, 7, 4), (1, 2, 6, 5)
_Y_MINUS, _Y_PLUS = (0, 1, 5, 4), (2, 3, 7, 6)
_BODY = tuple(range(8))
# gmsh reference manual, "Node ordering".
_GMSH_HEX20 = _HEX_CORNERS_ONLY + [
    (0, 1), (0, 3), (0, 4), (1, 2), (1, 5), (2, 3),
    (2, 6), (3, 7), (4, 5), (4, 7), (5, 6), (6, 7),
]  # fmt: skip
_GMSH_WED15 = _WED_CORNERS_ONLY + [
    (0, 1), (0, 2), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (3, 5), (4, 5),
]  # fmt: skip
_GMSH_PYR13 = _PYR_CORNERS_ONLY + [
    (0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (2, 4), (3, 4),
]  # fmt: skip

_FILE_ORDER = {
    ("gmsh", "tetra10"): [
        (0,),
        (1,),
        (2,),
        (3,),
        (0, 1),
        (1, 2),
        (0, 2),
        (0, 3),
        (2, 3),
        (1, 3),
    ],
    ("gmsh", "hexahedron20"): _GMSH_HEX20,
    ("gmsh", "hexahedron27"): _GMSH_HEX20
    + [_Z_MINUS, _Y_MINUS, _X_MINUS, _X_PLUS, _Y_PLUS, _Z_PLUS, _BODY],
    ("gmsh", "wedge15"): _GMSH_WED15,
    ("gmsh", "wedge18"): _GMSH_WED15 + [(0, 1, 4, 3), (0, 2, 5, 3), (1, 2, 5, 4)],
    ("gmsh", "pyramid13"): _GMSH_PYR13,
    ("gmsh", "pyramid14"): _GMSH_PYR13 + [(0, 1, 2, 3)],
    # CGNS SIDS, "Unstructured Grid Element Numbering Conventions".
    ("cgns", "hexahedron20"): _HEX20_VERTICALS_FIRST,
    ("cgns", "hexahedron27"): _HEX20_VERTICALS_FIRST
    + [_Z_MINUS, _Y_MINUS, _X_PLUS, _Y_PLUS, _X_MINUS, _Z_PLUS, _BODY],
    ("cgns", "wedge15"): _WED15_VERTICALS_FIRST,
    ("cgns", "wedge18"): _WED15_VERTICALS_FIRST
    + [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    # Kratos's geometry classes (hexahedra_3d_27.h, prism_3d_15.h), which the
    # GiD writer emits unchanged for these two types and `.mdpa` stores.
    ("gid", "hexahedron27"): _HEX20_VERTICALS_FIRST
    + [_Z_MINUS, _Y_MINUS, _X_PLUS, _Y_PLUS, _X_MINUS, _Z_PLUS, _BODY],
    ("gid", "wedge15"): _WED15_VERTICALS_FIRST,
    ("mdpa", "hexahedron20"): _HEX20_VERTICALS_FIRST,
    ("mdpa", "hexahedron27"): _HEX20_VERTICALS_FIRST
    + [_Z_MINUS, _Y_MINUS, _X_PLUS, _Y_PLUS, _X_MINUS, _Z_PLUS, _BODY],
    ("mdpa", "wedge15"): _WED15_VERTICALS_FIRST,
    # SEACAS Ioss (Hex20, Hex27, Wedge15): Hex27 puts the body centre first.
    ("exodus", "hexahedron20"): _HEX20_VERTICALS_FIRST,
    ("exodus", "hexahedron27"): _HEX20_VERTICALS_FIRST
    + [_BODY, _Z_MINUS, _Z_PLUS, _X_MINUS, _X_PLUS, _Y_MINUS, _Y_PLUS],
    ("exodus", "wedge15"): _WED15_VERTICALS_FIRST,
}


@pytest.mark.parametrize(
    "fmt, cell_type",
    sorted(_FILE_ORDER),
    ids=[f"{f}-{t}" for f, t in sorted(_FILE_ORDER)],
)
def test_a_file_order_element_reads_as_a_valid_element(fmt, cell_type):
    corners = np.array(_SHAPES[cell_type][0], dtype=float)
    in_file = np.array(
        [corners[list(s)].mean(0) for s in _FILE_ORDER[(fmt, cell_type)]]
    )
    assert len(in_file) == len(_reference(cell_type))
    back = _node_order.to_meshio(fmt, cell_type, np.arange(len(in_file))[None])[0]
    assert _is_valid(in_file[back], cell_type)


def test_every_moved_format_has_a_file_order_reference():
    # The tables moved out of the gmsh, CGNS, GiD, Kratos and Exodus readers
    # are each pinned by a transcription of their own format's documentation.
    moved = {"gmsh", "cgns", "gid", "mdpa", "exodus"}
    keys = {key for key in _ALL if key[0] in moved}
    assert keys == set(_FILE_ORDER)


def test_matches_the_cpp_twin():
    from meshioplusplus import _core

    table = _core.node_order_table()
    python = {
        key: (
            list(_node_order.node_order(*key).to_meshio),
            list(_node_order.node_order(*key).from_meshio),
        )
        for key in _ALL
    }
    assert {key: (list(a), list(b)) for key, (a, b) in table.items()} == python


# Code_Aster's gmsh reader: aster[j] = gmsh[_INIGMS[type][j]] (inigms.F90's
# nuconn, 0-based; gmsh's pyramid13 and wedge15/18 and hexahedron20/27).
_INIGMS = {
    "tetra10": [0, 1, 2, 3, 4, 5, 6, 7, 9, 8],
    "hexahedron20": list(range(9)) + [11, 13, 9, 10, 12, 14, 15, 16, 18, 19, 17],
    "hexahedron27": list(range(9))
    + [11, 13, 9, 10, 12, 14, 15, 16, 18, 19, 17]
    + [20, 21, 23, 24, 22, 25, 26],
    "wedge15": list(range(7)) + [9, 7, 8, 10, 11, 12, 14, 13],
    "wedge18": list(range(7)) + [9, 7, 8, 10, 11, 12, 14, 13, 15, 17, 16],
    "pyramid13": list(range(6)) + [8, 10, 6, 7, 9, 11, 12],
}


@pytest.mark.parametrize("cell_type", sorted(_INIGMS))
def test_code_aster_table_is_code_asters_gmsh_reader(cell_type):
    gmsh_row = np.arange(len(_INIGMS[cell_type]))[None]
    aster_row = gmsh_row[:, _INIGMS[cell_type]]
    via_gmsh = _gmsh_to_meshio_order(cell_type, gmsh_row)
    via_aster = _node_order.to_meshio("code_aster", cell_type, aster_row)
    np.testing.assert_array_equal(via_aster, via_gmsh)


# Code_Aster's MED reader: med[k] = aster[_LRMTYP[type][k]] (lrmtyp.F90's
# nuanom, 0-based); the types it leaves alone are absent.
_LRMTYP = {
    # fmt: off
    "tetra": [0, 2, 1, 3],
    "tetra10": [0, 2, 1, 3, 6, 5, 4, 7, 9, 8],
    "wedge": [0, 2, 1, 3, 5, 4],
    "wedge15": [0, 2, 1, 3, 5, 4, 8, 7, 6, 14, 13, 12, 9, 11, 10],
    "wedge18": [0, 2, 1, 3, 5, 4, 8, 7, 6, 14, 13, 12, 9, 11, 10, 17, 16, 15],
    "pyramid": [0, 3, 2, 1, 4],
    "pyramid13": [0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10],
    "hexahedron": [0, 3, 2, 1, 4, 7, 6, 5],
    "hexahedron20": [0, 3, 2, 1, 4, 7, 6, 5, 11, 10, 9, 8,
                     19, 18, 17, 16, 12, 15, 14, 13],
    "hexahedron27": [0, 3, 2, 1, 4, 7, 6, 5, 11, 10, 9, 8,
                     19, 18, 17, 16, 12, 15, 14, 13, 20, 24, 23, 22, 21, 25, 26],
    # fmt: on
}


@pytest.mark.parametrize("cell_type", sorted(_LRMTYP))
def test_code_aster_and_med_tables_agree_with_code_asters_med_reader(cell_type):
    # A .mail element, converted to MED the way Code_Aster does, must read back
    # through the MED table as a valid element: the same cell, possibly with its
    # corners relabelled by a symmetry of the reference element.
    ref = _reference(cell_type)
    aster = ref[list(_aster_from_meshio(cell_type))]
    med = aster[_LRMTYP[cell_type]]
    back = _node_order.to_meshio("med", cell_type, np.arange(len(ref))[None])[0]
    assert _is_valid(med[back], cell_type)


def _aster_from_meshio(cell_type):
    order = _node_order.node_order("code_aster", cell_type)
    return order.from_meshio if order else range(len(_reference(cell_type)))
