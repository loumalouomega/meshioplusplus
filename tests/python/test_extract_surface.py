import numpy as np
import pytest

import meshioplusplus
from meshioplusplus._surface import _extract_surface_py

from . import helpers


def _count(mesh, cell_type):
    return sum(len(b.data) for b in mesh.cells if b.type == cell_type)


def test_single_hex_six_quads():
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )
    mesh = meshioplusplus.Mesh(
        pts, [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))]
    )
    s = meshioplusplus.extract_surface(mesh)
    assert _count(s, "quad") == 6
    assert len(s.points) == 8
    np.testing.assert_allclose(np.sort(s.points, axis=0), np.sort(pts, axis=0))


def test_two_hex_stack_drops_internal_face():
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
            [0, 0, 2],
            [1, 0, 2],
            [1, 1, 2],
            [0, 1, 2],
        ],
        dtype=float,
    )
    conn = np.array([[0, 1, 2, 3, 4, 5, 6, 7], [4, 5, 6, 7, 8, 9, 10, 11]])
    mesh = meshioplusplus.Mesh(pts, [("hexahedron", conn)])
    s = meshioplusplus.extract_surface(mesh, record_parent_ids=True)
    assert _count(s, "quad") == 10  # 12 - 2 shared
    parents = s.cell_data["surface:parent_cell"][0].ravel()
    assert set(parents.tolist()) <= {0, 1}


def test_tetra_four_triangles():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    mesh = meshioplusplus.Mesh(pts, [("tetra", np.array([[0, 1, 2, 3]]))])
    s = meshioplusplus.extract_surface(mesh)
    assert _count(s, "triangle") == 4
    assert len(s.points) == 4


def test_2d_patch_boundary_edges():
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))])
    s = meshioplusplus.extract_surface(mesh)
    assert _count(s, "line") == 4


def test_no_supported_cells_raises():
    pts = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    mesh = meshioplusplus.Mesh(pts, [("line", np.array([[0, 1]]))])
    with pytest.raises(ValueError):
        meshioplusplus.extract_surface(mesh)


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tet_mesh,
        helpers.hex_mesh,
        helpers.wedge_mesh,
        helpers.pyramid_mesh,
        helpers.tet10_mesh,
        helpers.hex20_mesh,
        helpers.tri_mesh,
        helpers.quad_mesh,
        helpers.tri_quad_mesh,
    ],
)
@pytest.mark.parametrize("parents", [False, True])
def test_cpp_matches_python(mesh, parents):
    from meshioplusplus import _core

    a = _core.extract_surface(mesh, parents)
    b = _extract_surface_py(mesh, parents)
    assert [bl.type for bl in a.cells] == [bl.type for bl in b.cells]
    for ca, cb in zip(a.cells, b.cells):
        np.testing.assert_array_equal(ca.data, cb.data)
    np.testing.assert_allclose(a.points, b.points)
    if parents:
        for pa, pb in zip(
            a.cell_data["surface:parent_cell"], b.cell_data["surface:parent_cell"]
        ):
            np.testing.assert_array_equal(pa, pb)
