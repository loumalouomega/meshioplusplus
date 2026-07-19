"""Tests for the mesh cleanup operation (weld / prune / de-dup)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import clean
from meshioplusplus._clean import _clean_py


def test_weld_reduces_point_count():
    # node 3 coincides with node 0
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 0]], float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2], [3, 1, 2]]))])
    out, rep = clean(
        mesh,
        weld=True,
        atol=1e-9,
        remove_orphans=False,
        drop_degenerate=False,
        drop_duplicate_cells=False,
        return_report=True,
    )
    assert len(out.points) == 3
    assert rep["points_welded"] == 1
    # connectivity stays in range
    for cb in out.cells:
        assert cb.data.max() < len(out.points)


def test_remove_orphans():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [9, 9, 9]], float)  # node 3 orphan
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2]]))])
    out, rep = clean(mesh, remove_orphans=True, return_report=True)
    assert len(out.points) == 3
    assert rep["points_removed_orphan"] == 1


def test_drop_degenerate_and_duplicate():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
    mesh = meshioplusplus.Mesh(
        pts,
        [
            (
                "triangle",
                np.array(
                    [
                        [0, 1, 2],  # keep
                        [0, 1, 2],  # exact duplicate -> dropped
                        [0, 0, 1],  # repeated node -> degenerate
                    ]
                ),
            )
        ],
    )
    out, rep = clean(
        mesh, drop_degenerate=True, drop_duplicate_cells=True, return_report=True
    )
    assert rep["cells_dropped_duplicate"] == 1
    assert rep["cells_dropped_degenerate"] == 1
    assert len(out.cells[0].data) == 1
    for cb in out.cells:
        assert cb.data.max() < len(out.points)


def test_default_no_weld():
    # coincident nodes are NOT welded with the default options
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 0]], float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2], [3, 1, 2]]))])
    out, rep = clean(mesh, return_report=True)
    assert rep["points_welded"] == 0
    # both triangles kept (not identical without welding)
    assert len(out.cells[0].data) == 2


def test_cpp_matches_python():
    core = pytest.importorskip("meshioplusplus._core")
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 0], [5, 5, 5]], float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2], [3, 1, 2]]))])
    res = core.clean(mesh, True, 1e-9, True, True, True)
    ref, _, _, rep = _clean_py(mesh, True, 1e-9, True, True, True)
    assert np.allclose(res["mesh"].points, ref.points)
    assert res["points_welded"] == rep["points_welded"]
    assert res["points_removed_orphan"] == rep["points_removed_orphan"]


def test_sets_remapped():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [9, 9, 9]], float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2]]))])
    mesh.point_sets = {"corner": np.array([0, 3])}  # node 3 is an orphan
    out = clean(mesh, remove_orphans=True)
    # orphan node dropped from the set; node 0 remapped to 0
    assert np.array_equal(out.point_sets["corner"], [0])


def test_roundtrip_write_read(tmp_path):
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 0]], float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2], [3, 1, 2]]))])
    out = clean(mesh, weld=True, atol=1e-9)
    p = tmp_path / "c.vtu"
    meshioplusplus.write(p, out)
    back = meshioplusplus.read(p)
    assert len(back.points) == len(out.points)
