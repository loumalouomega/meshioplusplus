"""Tests for polyhedral coarsening (the ``agglomerate`` operation)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import agglomerate

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(
    _core is None, reason="agglomerate has no pure-Python fallback"
)


def _two_hexes():
    # Two unit hexahedra sharing exactly one face (x=1 plane), spanning
    # x in [0,1] and [1,2].
    pts = np.array(
        [
            [0, 0, 0],
            [0, 1, 0],
            [0, 1, 1],
            [0, 0, 1],
            [1, 0, 0],
            [1, 1, 0],
            [1, 1, 1],
            [1, 0, 1],
            [2, 0, 0],
            [2, 1, 0],
            [2, 1, 1],
            [2, 0, 1],
        ],
        float,
    )
    conn = np.array([[0, 1, 2, 3, 4, 5, 6, 7], [4, 5, 6, 7, 8, 9, 10, 11]])
    return meshioplusplus.Mesh(pts, [("hexahedron", conn)])


@needs_core
def test_two_adjacent_hexes_merge_into_one_polyhedron_conserving_volume():
    mesh = _two_hexes()
    out = agglomerate(mesh, target_group_size=2)
    assert len(out.cells) == 1
    assert out.cells[0].type == "polyhedron"
    assert len(out.cells[0].data) == 1

    before = meshioplusplus.compute_stats(mesh)["signed_volume"]
    after = meshioplusplus.compute_stats(out)["signed_volume"]
    assert before == 2.0
    assert after == 2.0


@needs_core
def test_identity_grouping_round_trips_a_mesh_with_a_boundary_block_first():
    hexes = _two_hexes()
    pts = hexes.points
    quad = meshioplusplus.Mesh(pts, [("quad", np.array([[0, 1, 2, 3]]))])
    mesh = meshioplusplus.merge([quad, hexes])
    mesh.cell_data["material"] = [np.array([100]), np.array([7, 8])]

    out = agglomerate(mesh, target_group_size=1)
    assert len(out.cells) == 2
    assert out.cells[0].type == "quad"
    assert len(out.cells[0].data) == 1
    assert out.cells[1].type == "polyhedron"
    assert len(out.cells[1].data) == 2

    assert meshioplusplus.compute_stats(out)["signed_volume"] == pytest.approx(2.0)

    assert list(out.cell_data["material"][0]) == [100]
    assert list(out.cell_data["material"][1]) == [7, 8]


@needs_core
def test_zero_target_group_size_raises_value_error():
    with pytest.raises(ValueError):
        agglomerate(_two_hexes(), target_group_size=0)


@needs_core
def test_a_non_manifold_face_is_refused_by_name():
    # Three tetrahedra sharing one triangular face.
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1], [1, 1, 1]],
        float,
    )
    conn = np.array([[0, 1, 2, 3], [0, 2, 1, 4], [0, 1, 2, 5]])
    mesh = meshioplusplus.Mesh(pts, [("tetra", conn)])
    with pytest.raises(ValueError):
        agglomerate(mesh)


@needs_core
def test_non_volume_blocks_pass_through_unchanged_when_no_volume_cells_exist():
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    out = agglomerate(mesh)
    assert len(out.cells) == 1
    assert out.cells[0].type == "triangle"
    assert len(out.cells[0].data) == len(mesh.cells[0].data)


@needs_core
def test_point_and_cell_regions_survive_but_side_regions_are_dropped_empty():
    mesh = _two_hexes()
    mesh.point_sets = {"corner": np.array([0], dtype=np.int64)}
    mesh.cell_sets = {"both": [np.array([0, 1], dtype=np.int64)]}

    out = agglomerate(mesh, target_group_size=2)
    assert list(out.point_sets["corner"]) == [0]
    assert list(out.cell_sets["both"][0]) == [0]


def test_agglomerate_is_exported():
    assert "agglomerate" in meshioplusplus.__all__
    assert meshioplusplus.agglomerate is agglomerate
