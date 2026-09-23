"""The node-ordering registry (``_node_order.py``, twin of ``detail/node_order.cpp``).

Every table is checked geometrically. A meshio++ reference element (mid-edge
nodes on edge midpoints, face and body centres on centroids) is written into the
file's order and read back through another table; the result must again be a
valid element. For Code_Aster the tables are also pinned against Code_Aster's
own readers:

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
def test_every_table_round_trips_a_reference_element(fmt, cell_type):
    order = _node_order.node_order(fmt, cell_type)
    ref = _reference(cell_type)
    assert sorted(order.to_meshio) == list(range(len(ref)))
    in_file = ref[list(order.from_meshio)]
    back = in_file[list(order.to_meshio)]
    np.testing.assert_array_equal(back, ref)
    assert _is_valid(back, cell_type)


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
