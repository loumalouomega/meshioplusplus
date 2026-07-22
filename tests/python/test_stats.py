"""Tests for the geometric statistics operation."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import compute_stats
from meshioplusplus._stats import _compute_stats_py


def _cube_hex():
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


def _cube_tets():
    # unit cube split into 6 tetrahedra (Kuhn triangulation)
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
    tets = np.array(
        [
            [0, 1, 2, 6],
            [0, 2, 3, 6],
            [0, 3, 7, 6],
            [0, 7, 4, 6],
            [0, 4, 5, 6],
            [0, 5, 1, 6],
        ]
    )
    return meshioplusplus.Mesh(pts, [("tetra", tets)])


def test_bbox_and_centroid():
    mesh = _cube_hex()
    s = compute_stats(mesh)
    assert s["bbox_min"] == (0.0, 0.0, 0.0)
    assert s["bbox_max"] == (1.0, 1.0, 1.0)
    assert s["extent"] == (1.0, 1.0, 1.0)
    assert np.allclose(s["centroid"], (0.5, 0.5, 0.5))
    assert s["num_points"] == 8
    assert s["cell_type_counts"] == {"hexahedron": 1}


def test_unit_cube_volume_hex():
    s = compute_stats(_cube_hex())
    assert abs(s["signed_volume"] - 1.0) < 1e-9
    assert abs(s["unsigned_volume"] - 1.0) < 1e-9
    assert abs(s["total_area"] - 6.0) < 1e-9  # cube surface area
    assert s["num_inverted"] == 0


def test_unit_cube_volume_tets():
    s = compute_stats(_cube_tets())
    assert abs(s["signed_volume"] - 1.0) < 1e-9
    assert abs(s["unsigned_volume"] - 1.0) < 1e-9
    assert s["num_inverted"] == 0


def test_triangle_area():
    pts = np.array([[0, 0, 0], [2, 0, 0], [0, 3, 0]], float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2]]))])
    s = compute_stats(mesh)
    assert abs(s["total_area"] - 3.0) < 1e-9  # 0.5 * 2 * 3


def test_inverted_element_counted():
    # a hexahedron with swapped top/bottom is inverted (negative volume)
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
    conn = np.array([[4, 5, 6, 7, 0, 1, 2, 3]])  # inverted winding
    mesh = meshioplusplus.Mesh(pts, [("hexahedron", conn)])
    s = compute_stats(mesh)
    assert s["signed_volume"] < 0
    assert s["num_inverted"] == 1


def test_cpp_matches_python():
    pytest.importorskip("meshioplusplus._core")
    from meshioplusplus import _core

    for mesh in (_cube_hex(), _cube_tets()):
        got = _core.compute_stats(mesh)
        ref = _compute_stats_py(mesh)
        assert abs(got["total_area"] - ref["total_area"]) < 1e-9
        assert abs(got["signed_volume"] - ref["signed_volume"]) < 1e-9
        assert abs(got["unsigned_volume"] - ref["unsigned_volume"]) < 1e-9
        assert got["num_inverted"] == ref["num_inverted"]
        assert np.allclose(got["centroid"], ref["centroid"])
