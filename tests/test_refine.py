"""Uniform refinement (``meshioplusplus.refine``)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import refine


def _two_triangles():
    """Two triangles sharing edge (0, 2): 4 points, 5 unique edges."""
    return meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )


def _unit_cube():
    return meshioplusplus.Mesh(
        np.array(
            [
                [0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
            ]
        ),
        [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))],
    )


def _unit_tet():
    return meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )


def _unit_quad():
    return meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0]]),
        [("quad", np.array([[0, 1, 2, 3]]))],
    )


def _signed_tet_volume(pts, row):
    a, b, c, d = (pts[i] for i in row)
    return np.dot(np.cross(b - a, c - a), d - a) / 6.0


# --- 2D --------------------------------------------------------------------- #


def test_triangle_splits_into_four_with_shared_mid_edge_nodes():
    out = refine(_two_triangles())
    assert len(out.cells) == 1
    assert out.cells[0].type == "triangle"
    assert len(out.cells[0].data) == 8
    # 4 original + 5 unique edges -- NOT 4 + 2*3, which is what an undeduped
    # per-cell midpoint would give.
    assert len(out.points) == 9


def test_triangle_mid_nodes_sit_at_edge_midpoints():
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [2, 0, 0], [0, 2, 0]]),
        [("triangle", np.array([[0, 1, 2]]))],
    )
    out = refine(mesh)
    assert len(out.points) == 6
    # New nodes are appended in edge order (0,1), (1,2), (2,0).
    assert np.allclose(out.points[3], [1, 0, 0])
    assert np.allclose(out.points[4], [1, 1, 0])
    assert np.allclose(out.points[5], [0, 1, 0])


def test_linear_point_data_is_interpolated_exactly():
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [2, 0, 0], [0, 2, 0]]),
        [("triangle", np.array([[0, 1, 2]]))],
    )
    # f(x, y) = 3x + 5y + 1
    mesh.point_data["f"] = np.array([1.0, 7.0, 11.0])
    out = refine(mesh)
    expected = 3 * out.points[:, 0] + 5 * out.points[:, 1] + 1
    assert np.allclose(out.point_data["f"], expected)


def test_quad_splits_into_four_with_centre_at_corner_mean():
    out = refine(_unit_quad())
    assert out.cells[0].type == "quad"
    assert len(out.cells[0].data) == 4
    # 4 corners + 4 edge mids + 1 face centre.
    assert len(out.points) == 9
    assert np.allclose(out.points[8], [1, 1, 0])


# --- 3D --------------------------------------------------------------------- #


def test_tetra_splits_into_eight_conserving_volume_and_orientation():
    mesh = _unit_tet()
    before = meshioplusplus.compute_stats(mesh)
    out = refine(mesh)
    after = meshioplusplus.compute_stats(out)

    assert out.cells[0].type == "tetra"
    assert len(out.cells[0].data) == 8
    assert after["signed_volume"] == pytest.approx(before["signed_volume"])
    assert after["num_inverted"] == 0

    pts = np.asarray(out.points)
    for row in out.cells[0].data:
        assert _signed_tet_volume(pts, row) > 0


def test_hex_splits_into_eight_with_27_nodes():
    mesh = _unit_cube()
    before = meshioplusplus.compute_stats(mesh)
    out = refine(mesh)
    after = meshioplusplus.compute_stats(out)

    assert out.cells[0].type == "hexahedron"
    assert len(out.cells[0].data) == 8
    # 8 corners + 12 edge mids + 6 face centres + 1 body.
    assert len(out.points) == 27
    assert np.allclose(out.points[26], [0.5, 0.5, 0.5])
    assert after["signed_volume"] == pytest.approx(before["signed_volume"])
    assert after["num_inverted"] == 0


def test_wedge_splits_into_eight_with_18_nodes():
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]]),
        [("wedge", np.array([[0, 1, 2, 3, 4, 5]]))],
    )
    before = meshioplusplus.compute_stats(mesh)
    out = refine(mesh)
    after = meshioplusplus.compute_stats(out)

    assert out.cells[0].type == "wedge"
    assert len(out.cells[0].data) == 8
    # 6 corners + 9 edge mids + 3 quad-face centres, no body node.
    assert len(out.points) == 18
    assert after["signed_volume"] == pytest.approx(before["signed_volume"])
    assert after["num_inverted"] == 0


def test_adjacent_hexes_refine_conformingly():
    """Face-centre nodes must be shared, or every interior face becomes boundary."""
    mesh = meshioplusplus.Mesh(
        np.array(
            [
                [0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
                [0, 0, 2], [1, 0, 2], [1, 1, 2], [0, 1, 2],
            ]
        ),
        [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7], [4, 5, 6, 7, 8, 9, 10, 11]]))],
    )
    out = refine(mesh)
    # A conforming 2x2x4 lattice of hexahedra has 3x3x5 nodes.
    assert len(out.points) == 45

    before = sum(len(cb.data) for cb in meshioplusplus.extract_surface(mesh).cells)
    after = sum(len(cb.data) for cb in meshioplusplus.extract_surface(out).cells)
    # Every boundary facet splits into 4; no interior facet may show up.
    assert after == before * 4


# --- levels, data, maps ------------------------------------------------------ #


def test_levels_two_equals_refining_twice():
    mesh = _two_triangles()
    direct = refine(mesh, levels=2)
    twice = refine(refine(mesh))
    assert np.array_equal(direct.cells[0].data, twice.cells[0].data)
    assert np.allclose(direct.points, twice.points)


def test_zero_levels_is_unchanged():
    mesh = _two_triangles()
    out = refine(mesh, levels=0)
    assert np.array_equal(out.cells[0].data, mesh.cells[0].data)
    assert np.allclose(out.points, mesh.points)


def test_cell_data_is_replicated_to_children():
    mesh = _two_triangles()
    mesh.cell_data["tag"] = [np.array([10, 20])]
    out = refine(mesh)
    assert list(out.cell_data["tag"][0]) == [10, 10, 10, 10, 20, 20, 20, 20]


def test_record_parent_ids_names_the_original_ancestor():
    out = refine(_two_triangles(), levels=2, record_parent_ids=True)
    parents = np.asarray(out.cell_data["refine:parent_cell"][0]).reshape(-1)
    assert len(parents) == 32  # 2 parents x 16 grandchildren
    assert set(parents[:16]) == {0}
    assert set(parents[16:]) == {1}


def test_cell_sets_expand_to_children():
    mesh = _two_triangles()
    mesh.cell_sets = {"left": [np.array([0])]}
    out = refine(mesh)
    assert list(out.cell_sets["left"][0]) == [0, 1, 2, 3]


def test_mixed_type_mesh_refines_each_block_by_its_own_template():
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [2, 0, 0], [2, 1, 0], [1, 1, 0], [0, 1, 0]]),
        [
            ("triangle", np.array([[0, 1, 4], [0, 4, 5]])),
            ("quad", np.array([[1, 2, 3, 4]])),
        ],
    )
    out = refine(mesh)
    assert [cb.type for cb in out.cells] == ["triangle", "quad"]
    assert len(out.cells[0].data) == 8
    assert len(out.cells[1].data) == 4


# --- rejections, parity, I/O ------------------------------------------------- #


@pytest.mark.parametrize(
    "cell_type, conn",
    [
        ("tetra10", np.arange(10).reshape(1, 10)),
        ("triangle6", np.arange(6).reshape(1, 6)),
        ("pyramid", np.arange(5).reshape(1, 5)),
    ],
)
def test_unsupported_cell_types_raise(cell_type, conn):
    mesh = meshioplusplus.Mesh(np.zeros((16, 3)), [(cell_type, conn)])
    with pytest.raises(ValueError):
        refine(mesh)


def test_ragged_block_raises():
    mesh = meshioplusplus.Mesh(
        np.zeros((5, 3)), [("polygon", [[0, 1, 2, 3], [3, 2, 4]])]
    )
    with pytest.raises(ValueError):
        refine(mesh)


@pytest.mark.parametrize(
    "factory", [_two_triangles, _unit_quad, _unit_tet, _unit_cube]
)
def test_cpp_matches_python(factory):
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._refine import _refine_py

    mesh = factory()
    got = core.refine(mesh, 2, False)["mesh"]
    ref, _, _ = _refine_py(mesh, 2, False)

    assert [cb.type for cb in got.cells] == [cb.type for cb in ref.cells]
    for a, b in zip(got.cells, ref.cells):
        assert np.array_equal(a.data, b.data)
    assert np.allclose(got.points, ref.points)


def test_roundtrip_write_read(tmp_path):
    out = refine(_unit_cube())
    path = tmp_path / "refined.vtu"
    meshioplusplus.write(path, out)
    back = meshioplusplus.read(path)
    assert len(back.points) == len(out.points)
    assert np.array_equal(back.cells[0].data, out.cells[0].data)
