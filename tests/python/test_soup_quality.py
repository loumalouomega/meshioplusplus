"""Tests for ``_sdf.soup_quality``, the watertightness report ``compute_sdf``'s
check and ``soup_quality`` (``_normals``) share."""

import numpy as np

from meshioplusplus._normals import _quality
from meshioplusplus._sdf import soup_quality


def _tet_faces():
    return np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]])


def test_a_closed_consistent_surface_is_watertight():
    r = soup_quality(_tet_faces(), 0)
    assert r == {
        "boundary_edges": 0,
        "non_manifold_edges": 0,
        "inconsistent_pairs": 0,
        "degenerate_triangles": 0,
        "watertight": True,
    }


def test_each_defect_is_counted():
    open_surface = _tet_faces()[:3]
    assert soup_quality(open_surface, 0)["boundary_edges"] == 3
    flipped = _tet_faces().copy()
    flipped[0] = flipped[0][::-1]
    assert soup_quality(flipped, 0)["inconsistent_pairs"] == 3
    fin = np.vstack([_tet_faces(), [[0, 1, 4]]])
    r = soup_quality(fin, 0)
    assert r["non_manifold_edges"] == 1 and not r["watertight"]
    # The caller's degenerate count decides the verdict too.
    assert not soup_quality(_tet_faces(), 1)["watertight"]
    assert soup_quality(np.zeros((0, 3), dtype=np.int64), 0)["watertight"]


def test_normals_quality_counts_degenerates_itself():
    corners = np.array(
        [
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            [[0, 0, 0], [1, 0, 0], [2, 0, 0]],  # collinear: zero area
        ],
        dtype=float,
    )
    verts = np.array([[0, 1, 2], [0, 1, 3]])
    r = _quality(verts, corners)
    assert r["degenerate_triangles"] == 1
    assert r["boundary_edges"] == soup_quality(verts, 1)["boundary_edges"]
