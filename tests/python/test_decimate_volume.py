"""Tests for volume decimation (the ``decimate_volume`` operation)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import decimate, decimate_volume

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(
    _core is None, reason="decimate_volume has no pure-Python fallback"
)


def _cube6():
    # A unit cube split into 6 positively-oriented tets sharing the main
    # diagonal 0-6 -- every vertex is a cube corner (a sharp feature).
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
    conn = np.array(
        [
            [0, 1, 2, 6],
            [0, 2, 3, 6],
            [0, 3, 7, 6],
            [0, 7, 4, 6],
            [0, 4, 5, 6],
            [0, 5, 1, 6],
        ]
    )
    return meshioplusplus.Mesh(pts, [("tetra", conn)])


@needs_core
def test_target_cells_is_honoured():
    mesh = _cube6()
    out, report = decimate_volume(
        mesh, target_cells=1, preserve_features=False, return_report=True
    )
    assert len(out.cells[0].data) <= 1
    assert report["tets_removed"] > 0


@needs_core
def test_exactly_one_stopping_criterion_is_required():
    mesh = _cube6()
    with pytest.raises(ValueError):
        decimate_volume(mesh)
    with pytest.raises(ValueError):
        decimate_volume(mesh, ratio=0.5, target_cells=1)


@needs_core
def test_preserve_boundary_pins_every_vertex_of_an_all_boundary_mesh():
    mesh = _cube6()
    out = decimate_volume(mesh, target_cells=1, preserve_boundary=True)
    assert len(out.cells[0].data) == len(mesh.cells[0].data)


@needs_core
def test_frozen_by_point_set_name_prevents_collapse():
    mesh = _cube6()
    mesh.point_sets = {"pin": np.arange(mesh.points.shape[0])}
    out = decimate_volume(mesh, target_cells=1, preserve_features=False, frozen="pin")
    assert len(out.cells[0].data) == len(mesh.cells[0].data)


@needs_core
def test_no_surviving_tet_inverts():
    mesh = _cube6()
    out = decimate_volume(mesh, ratio=0.4, preserve_features=False)
    assert meshioplusplus.compute_stats(out)["num_inverted"] == 0


@needs_core
def test_rejects_non_tetra_blocks_by_name():
    hexmesh = meshioplusplus.Mesh(
        np.array(
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
        ),
        [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))],
    )
    with pytest.raises(ValueError, match="tet-only"):
        decimate_volume(hexmesh, target_cells=1)


@needs_core
def test_decimate_on_a_volume_mesh_still_raises_pointing_here():
    with pytest.raises(ValueError, match="extract_surface"):
        decimate(_cube6(), ratio=0.5)


@needs_core
def test_data_is_carried_and_regions_survive():
    mesh = _cube6()
    mesh.point_data["scalar"] = np.arange(mesh.points.shape[0], dtype=float)
    mesh.cell_data["tag"] = [np.arange(len(mesh.cells[0].data), dtype=np.int64)]
    mesh.cell_sets = {"all": [np.arange(len(mesh.cells[0].data))]}

    out, report = decimate_volume(
        mesh, target_cells=1, preserve_features=False, return_report=True
    )
    assert "scalar" in out.point_data
    assert "tag" in out.cell_data
    assert "all" in out.cell_sets
    assert len(out.cell_sets["all"][0]) == len(out.cells[0].data)
    assert set(report) == {
        "tets_removed",
        "points_removed",
        "collapses_rejected",
        "max_error_applied",
    }
