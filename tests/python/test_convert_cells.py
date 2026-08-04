"""Tests for the convert_cells operation (linearize / simplexify / elevate)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import convert_cells


def _unit_cube():
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
        float,
    )
    return meshioplusplus.Mesh(
        pts, [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))]
    )


def _tetra10():
    # 4 corners of a unit tet plus the 6 exact edge midpoints, in VTK order:
    # 4=mid(0,1), 5=mid(1,2), 6=mid(0,2), 7=mid(0,3), 8=mid(1,3), 9=mid(2,3).
    corners = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    edges = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
    mids = np.array([0.5 * (corners[a] + corners[b]) for a, b in edges])
    pts = np.vstack([corners, mids])
    return meshioplusplus.Mesh(
        pts, [("tetra10", np.array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]))]
    )


def _quad_grid():
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0]], float
    )
    return meshioplusplus.Mesh(pts, [("quad", np.array([[0, 1, 4, 3], [1, 2, 5, 4]]))])


def _signed_tet_volume(pts, row):
    a, b, c, d = (pts[i] for i in row)
    return float(np.dot(np.cross(b - a, c - a), d - a)) / 6.0


# --------------------------------------------------------------------------- #
# linearize                                                                    #
# --------------------------------------------------------------------------- #
def test_linearize_tetra10_keeps_corners_and_prunes_mid_nodes():
    mesh = _tetra10()
    out = convert_cells(mesh, mode="linearize")

    assert len(out.cells) == 1
    assert out.cells[0].type == "tetra"
    assert len(out.cells[0].data) == 1
    # The 6 mid-edge nodes are unreferenced and pruned.
    assert len(out.points) == 4
    assert out.cells[0].data.tolist() == [[0, 1, 2, 3]]
    assert np.allclose(out.points, mesh.points[:4])


def test_linearize_carries_point_data_through_the_prune():
    mesh = _tetra10()
    mesh.point_data["T"] = np.arange(10, dtype=float)
    out = convert_cells(mesh, mode="linearize")
    assert np.allclose(out.point_data["T"], [0.0, 1.0, 2.0, 3.0])


def test_linearize_leaves_linear_cells_alone():
    mesh = _unit_cube()
    out = convert_cells(mesh, mode="linearize")
    assert out.cells[0].type == "hexahedron"
    assert len(out.points) == 8
    assert np.array_equal(out.cells[0].data, mesh.cells[0].data)


# --------------------------------------------------------------------------- #
# simplexify                                                                   #
# --------------------------------------------------------------------------- #
def test_simplexify_hexahedron_gives_six_oriented_tetra():
    mesh = _unit_cube()
    out = convert_cells(mesh, mode="simplexify")

    assert len(out.cells) == 1
    assert out.cells[0].type == "tetra"
    assert len(out.cells[0].data) == 6
    # The decomposition reuses the parent's own corner nodes.
    assert len(out.points) == 8

    vols = [_signed_tet_volume(out.points, row) for row in out.cells[0].data]
    assert all(v > 0 for v in vols), "simplexify emitted an inverted tetrahedron"
    assert np.isclose(sum(vols), 1.0)


def test_simplexify_quad_gives_two_triangles():
    mesh = _quad_grid()
    out = convert_cells(mesh, mode="simplexify")
    assert out.cells[0].type == "triangle"
    assert len(out.cells[0].data) == 4  # 2 quads -> 2 triangles each
    assert len(out.points) == len(mesh.points)


def test_simplexify_ngon_fans_into_n_minus_two_triangles():
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [1.5, 1, 0], [0.5, 1.5, 0], [-0.5, 1, 0]], float
    )
    mesh = meshioplusplus.Mesh(pts, [("polygon", np.array([[0, 1, 2, 3, 4]]))])
    out = convert_cells(mesh, mode="simplexify")
    assert out.cells[0].type == "triangle"
    assert len(out.cells[0].data) == 3


def test_simplexify_replicates_cell_data_one_array_per_block():
    mesh = _unit_cube()
    mesh.cell_data["mat"] = [np.array([7], dtype=np.int64)]
    out = convert_cells(mesh, mode="simplexify", record_parent_ids=True)

    # Exactly one cell_data array per output block, sized to the output cells.
    assert len(out.cell_data["mat"]) == len(out.cells)
    assert out.cell_data["mat"][0].tolist() == [7] * 6
    assert out.cell_data["convert:parent_cell"][0].tolist() == [0] * 6


def test_simplexify_polyhedron_fans_into_tetrahedra():
    """Since v9.17.0 a polyhedron decomposes instead of raising."""
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    mesh = meshioplusplus.Mesh(
        pts, [("polyhedron4", [[[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]]])]
    )
    out = convert_cells(mesh, mode="simplexify")
    assert len(out.cells) == 1
    assert out.cells[0].type == "tetra"
    # 4 faces x 3 edges = 12 children; 1 cell centroid + 4 face centroids added.
    assert len(out.cells[0].data) == 12
    assert len(out.points) == 4 + 5


# --------------------------------------------------------------------------- #
# elevate                                                                      #
# --------------------------------------------------------------------------- #
def test_elevate_triangle_puts_mid_nodes_at_edge_midpoints():
    pts = np.array([[0, 0, 0], [2, 0, 0], [0, 2, 0]], float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2]]))])
    mesh.point_data["T"] = np.array([0.0, 10.0, 20.0])

    out = convert_cells(mesh, mode="elevate")

    assert out.cells[0].type == "triangle6"
    assert len(out.points) == 6
    assert out.cells[0].data.tolist() == [[0, 1, 2, 3, 4, 5]]
    # node 3 = mid(0,1), node 4 = mid(1,2), node 5 = mid(2,0)
    assert np.allclose(out.points[3], [1, 0, 0])
    assert np.allclose(out.points[4], [1, 1, 0])
    assert np.allclose(out.points[5], [0, 1, 0])
    # point_data on a new node is the mean of its edge endpoints.
    assert np.allclose(out.point_data["T"][3:], [5.0, 15.0, 10.0])


def test_elevate_shares_mid_nodes_between_adjacent_cells():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2], [1, 3, 2]]))])
    out = convert_cells(mesh, mode="elevate")
    # 5 unique edges over the two triangles -> 5 new nodes.
    assert len(out.points) == 4 + 5
    # Both triangles are in one block; the shared edge (1,2) is cell 0's edge 1
    # and cell 1's edge 2.
    assert out.cells[0].data[0][4] == out.cells[0].data[1][5]


def test_elevate_full_lagrange_target_raises():
    pts = np.zeros((9, 3))
    mesh = meshioplusplus.Mesh(
        pts, [("quad9", np.array([[0, 1, 2, 3, 4, 5, 6, 7, 8]]))]
    )
    with pytest.raises(Exception):
        convert_cells(mesh, mode="elevate")


# --------------------------------------------------------------------------- #
# round-trip identities                                                        #
# --------------------------------------------------------------------------- #
def test_elevate_then_linearize_returns_the_original():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    mesh = meshioplusplus.Mesh(pts, [("tetra", np.array([[0, 1, 2, 3]]))])

    back = convert_cells(convert_cells(mesh, mode="elevate"), mode="linearize")

    assert back.cells[0].type == "tetra"
    assert np.array_equal(back.cells[0].data, mesh.cells[0].data)
    assert np.allclose(back.points, mesh.points)


def test_linearize_then_elevate_restores_the_topology():
    mesh = _tetra10()
    up = convert_cells(convert_cells(mesh, mode="linearize"), mode="elevate")

    assert up.cells[0].type == "tetra10"
    assert len(up.cells[0].data) == len(mesh.cells[0].data)
    assert len(up.points) == len(mesh.points)
    # The fixture's mid nodes are the exact edge midpoints, so the geometry
    # comes back too (mid nodes may be renumbered, hence the sorted compare).
    assert np.allclose(
        np.sort(up.points, axis=0), np.sort(np.asarray(mesh.points), axis=0)
    )


# --------------------------------------------------------------------------- #
# C++/numpy parity, sets, I/O, validation                                      #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mode", ["linearize", "simplexify", "elevate"])
def test_cpp_matches_python(mode):
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._convert_cells import _convert_cells_py

    mesh = _unit_cube()
    got = core.convert_cells(mesh, mode, False)["mesh"]
    ref, _, _ = _convert_cells_py(mesh, mode, False)

    assert [cb.type for cb in got.cells] == [cb.type for cb in ref.cells]
    for a, b in zip(got.cells, ref.cells):
        assert np.array_equal(a.data, b.data)
    assert np.allclose(got.points, ref.points)


def test_cell_sets_expand_to_children_under_simplexify():
    mesh = _quad_grid()
    mesh.cell_sets = {"left": [np.array([0])]}
    out = convert_cells(mesh, mode="simplexify")
    # Quad 0 becomes triangles 0 and 1.
    assert sorted(out.cell_sets["left"][0].tolist()) == [0, 1]


def test_roundtrip_write_read(tmp_path):
    out = convert_cells(_unit_cube(), mode="simplexify")
    p = tmp_path / "cc.vtu"
    meshioplusplus.write(p, out)
    back = meshioplusplus.read(p)
    assert len(back.points) == len(out.points)
    assert back.cells[0].type == "tetra"
    assert len(back.cells[0].data) == 6


def test_unknown_mode_raises():
    with pytest.raises(ValueError):
        convert_cells(_unit_cube(), mode="nope")
