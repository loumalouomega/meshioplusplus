"""Tests for surface remeshing (the ``remesh`` operation)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import remesh

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="remesh has no pure-Python fallback")


def _icosahedron():
    t = (1.0 + np.sqrt(5.0)) / 2.0
    pts = np.array(
        [
            [-1, t, 0],
            [1, t, 0],
            [-1, -t, 0],
            [1, -t, 0],
            [0, -1, t],
            [0, 1, t],
            [0, -1, -t],
            [0, 1, -t],
            [t, 0, -1],
            [t, 0, 1],
            [-t, 0, -1],
            [-t, 0, 1],
        ],
        dtype=float,
    )
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    faces = np.array(
        [
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ],
        dtype=np.int64,
    )
    return meshioplusplus.Mesh(pts, [("triangle", faces)])


@needs_core
def test_produces_the_requested_number_of_clusters():
    out, report = remesh(_icosahedron(), 150, return_report=True)
    assert out.points.shape[0] == 150
    assert report["num_clusters"] == 150
    assert out.cells[0].type == "triangle"
    assert out.cells[0].data.shape[0] > 0


@needs_core
def test_output_is_watertight():
    out = remesh(_icosahedron(), 150)
    conn = out.cells[0].data
    edges = {}
    for tri in conn:
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            key = (a, b) if a < b else (b, a)
            edges[key] = edges.get(key, 0) + 1
    assert set(edges.values()) == {2}


@needs_core
def test_quadric_metric_is_accepted_and_deterministic():
    a, ra = remesh(_icosahedron(), 80, metric="quadric", return_report=True)
    b, rb = remesh(_icosahedron(), 80, metric="quadric", return_report=True)
    assert a.points.shape == b.points.shape
    assert ra["num_iterations"] == rb["num_iterations"]
    np.testing.assert_array_equal(a.points, b.points)


@needs_core
def test_unknown_metric_raises():
    with pytest.raises(Exception):
        remesh(_icosahedron(), 50, metric="bogus")


@needs_core
def test_too_few_clusters_raises():
    with pytest.raises(Exception):
        remesh(_icosahedron(), 3)


@needs_core
def test_rejects_volume_cells_by_name():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    conn = np.array([[0, 1, 2, 3]], dtype=np.int64)
    mesh = meshioplusplus.Mesh(pts, [("tetra", conn)])
    with pytest.raises(Exception, match="extract_surface"):
        remesh(mesh, 10)


@needs_core
def test_point_data_is_dropped_field_data_is_carried():
    mesh = _icosahedron()
    mesh.point_data["temperature"] = np.full(mesh.points.shape[0], 3.0)
    mesh.field_data["solver"] = np.array([1.5])
    out = remesh(mesh, 40)
    assert "temperature" not in out.point_data
    assert "solver" in out.field_data


@needs_core
def test_subdivide_zero_disables_refinement():
    out, report = remesh(_icosahedron(), 10, subdivide=0, return_report=True)
    assert report["subdivide_applied"] == 0

