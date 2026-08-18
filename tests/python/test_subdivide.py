"""Tests for polyhedral refinement (the ``subdivide`` operation)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import subdivide

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(
    _core is None, reason="subdivide has no pure-Python fallback"
)


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


def _l_prism():
    # A non-convex L-shaped prism: [0,2]x[0,2] extruded to z=1 with the
    # [1,2]x[1,2] corner removed -- cross-section area 3, so volume 3. A
    # convex fixture cannot distinguish subdivide's corner-average apex from
    # `compute_stats`'/`poly_measure`'s volume centroid.
    pts = np.array(
        [
            [0, 0, 0],
            [2, 0, 0],
            [2, 1, 0],
            [1, 1, 0],
            [1, 2, 0],
            [0, 2, 0],
            [0, 0, 1],
            [2, 0, 1],
            [2, 1, 1],
            [1, 1, 1],
            [1, 2, 1],
            [0, 2, 1],
        ],
        float,
    )
    faces = [
        [0, 5, 4, 3, 2, 1],
        [6, 7, 8, 9, 10, 11],
        [0, 1, 7, 6],
        [1, 2, 8, 7],
        [2, 3, 9, 8],
        [3, 4, 10, 9],
        [4, 5, 11, 10],
        [5, 0, 6, 11],
    ]
    return meshioplusplus.Mesh(pts, [("polyhedron12", [faces])])


@needs_core
def test_hexahedron_subdivides_into_six_pyramidal_children():
    mesh = _unit_cube()
    out = subdivide(mesh)
    assert len(out.cells) == 1
    cb = out.cells[0]
    assert cb.type == "polyhedron"
    assert len(cb.data) == 6
    assert len(out.points) == len(mesh.points) + 1


@needs_core
def test_non_convex_polyhedron_conserves_volume():
    mesh = _l_prism()
    before = meshioplusplus.compute_stats(mesh)["signed_volume"]
    assert before == pytest.approx(3.0, abs=1e-12)

    out = subdivide(mesh)
    after = meshioplusplus.compute_stats(out)["signed_volume"]
    assert after == pytest.approx(before, abs=1e-9)


@needs_core
def test_record_parent_ids_attaches_the_array():
    mesh = _unit_cube()
    out = subdivide(mesh, record_parent_ids=True)
    assert "subdivide:parent_cell" in out.cell_data
    parents = out.cell_data["subdivide:parent_cell"][0]
    assert len(parents) == 6
    assert np.all(np.asarray(parents) == 0)


@needs_core
def test_non_three_d_blocks_pass_through_unchanged():
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    out = subdivide(mesh)
    assert len(out.cells) == 1
    assert out.cells[0].type == "triangle"
    assert len(out.cells[0].data) == 2
    assert len(out.points) == len(mesh.points)


@needs_core
def test_an_open_face_set_raises_value_error():
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
    # The top face {4,5,6,7} is deliberately omitted, so the faces do not
    # close.
    faces = [
        [0, 3, 2, 1],
        [0, 1, 5, 4],
        [2, 3, 7, 6],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
    ]
    mesh = meshioplusplus.Mesh(pts, [("polyhedron8", [faces])])
    with pytest.raises(ValueError):
        subdivide(mesh)


@needs_core
def test_point_and_cell_regions_survive_but_side_regions_are_dropped():
    mesh = _unit_cube()
    mesh.point_sets = {"corner": np.array([0], dtype=np.int64)}
    mesh.cell_sets = {"all": [np.array([0], dtype=np.int64)]}

    out = subdivide(mesh)
    assert list(out.point_sets["corner"]) == [0]
    assert sorted(int(c) for c in out.cell_sets["all"][0]) == [0, 1, 2, 3, 4, 5]


def test_subdivide_is_exported():
    assert "subdivide" in meshioplusplus.__all__
    assert meshioplusplus.subdivide is subdivide
