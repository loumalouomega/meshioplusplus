"""Tests for volumetric retetrahedralization (the ``remesh_volume`` operation)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import remesh_volume

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="remesh_volume has no pure-Python fallback")


def _octahedron():
    pts = np.array(
        [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float
    )
    faces = np.array(
        [
            [0, 2, 4],
            [2, 1, 4],
            [1, 3, 4],
            [3, 0, 4],
            [2, 0, 5],
            [1, 2, 5],
            [3, 1, 5],
            [0, 3, 5],
        ],
        dtype=np.int64,
    )
    return meshioplusplus.Mesh(pts, [("triangle", faces)])


def _box_volume(half=1.0):
    pts = np.array(
        [
            [-half, -half, -half],
            [half, -half, -half],
            [half, half, -half],
            [-half, half, -half],
            [-half, -half, half],
            [half, -half, half],
            [half, half, half],
            [-half, half, half],
        ],
        dtype=float,
    )
    tets = np.array(
        [[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]],
        dtype=np.int64,
    )
    return meshioplusplus.Mesh(pts, [("tetra", tets)])


@needs_core
def test_produces_a_single_tetra_block():
    out, report = remesh_volume(
        _octahedron(), cell_size=0.4, watertight_check="off", return_report=True
    )
    assert len(out.cells) == 1
    assert out.cells[0].type == "tetra"
    assert out.cells[0].data.shape[0] == report["num_tets"]
    assert report["num_tets"] > 0
    assert report["num_vertices_warped"] >= 0
    assert report["num_tets_rejected"] >= 0
    assert report["num_non_manifold_edges"] >= 0
    assert report["input_quality"]["watertight"]


@needs_core
def test_accepts_a_volume_mesh_by_extracting_its_boundary():
    out = remesh_volume(_box_volume(), cell_size=0.4, watertight_check="off")
    assert out.cells[0].data.shape[0] > 0


@needs_core
def test_resolution_and_cell_size_are_both_accepted():
    a = remesh_volume(_octahedron(), cell_size=0.4, watertight_check="off")
    b = remesh_volume(_octahedron(), resolution=(6, 6, 6), watertight_check="off")
    assert a.cells[0].data.shape[0] > 0
    assert b.cells[0].data.shape[0] > 0


@needs_core
def test_neither_resolution_nor_cell_size_raises():
    with pytest.raises(Exception):
        remesh_volume(_octahedron())


@needs_core
def test_both_resolution_and_cell_size_raises():
    with pytest.raises(Exception):
        remesh_volume(_octahedron(), resolution=(6, 6, 6), cell_size=0.4)


@needs_core
def test_negative_warp_fraction_raises():
    with pytest.raises(Exception):
        remesh_volume(_octahedron(), cell_size=0.4, warp_fraction=-0.1)


@needs_core
def test_oversized_output_is_refused_by_name():
    with pytest.raises(Exception, match="max_tets"):
        remesh_volume(_octahedron(), cell_size=0.02, max_tets=10, watertight_check="off")


@needs_core
def test_zero_warp_fraction_is_exactly_watertight():
    # mWarpFraction = 0 has no reuse-a-warped-position step, so the output's
    # own boundary is mathematically watertight -- doc/remesh_volume.md's
    # measured, honest oracle.
    out, report = remesh_volume(
        _octahedron(), cell_size=0.4, warp_fraction=0.0, watertight_check="off", return_report=True
    )
    assert report["num_non_manifold_edges"] == 0
    assert report["num_vertices_warped"] == 0


@needs_core
def test_point_data_is_dropped_field_data_is_carried():
    mesh = _octahedron()
    mesh.point_data["temperature"] = np.full(mesh.points.shape[0], 3.0)
    mesh.field_data["solver"] = np.array([1.5])
    out = remesh_volume(mesh, cell_size=0.4, watertight_check="off")
    assert "temperature" not in out.point_data
    assert "solver" in out.field_data


@needs_core
def test_result_is_stable_across_repeated_runs():
    mesh = _octahedron()
    a = remesh_volume(mesh, cell_size=0.3, watertight_check="off")
    b = remesh_volume(mesh, cell_size=0.3, watertight_check="off")
    np.testing.assert_array_equal(a.points, b.points)
    np.testing.assert_array_equal(a.cells[0].data, b.cells[0].data)
