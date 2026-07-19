"""Tests for the crop operation (bounding box / half-space)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import crop


def _quad_grid():
    # 2x1 grid of quads (3x2 points), x in [0,2], y in [0,1]
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0]], float
    )
    return meshioplusplus.Mesh(pts, [("quad", np.array([[0, 1, 4, 3], [1, 2, 5, 4]]))])


def test_bbox_keeps_cells_inside():
    mesh = _quad_grid()
    out = crop(mesh, bbox=[-0.1, -0.1, -0.1, 1.1, 1.1, 1.1], mode="all")
    assert len(out.cells[0].data) == 1  # only the first quad is fully inside
    # no unused points
    used = set(out.cells[0].data.reshape(-1).tolist())
    assert used == set(range(len(out.points)))


def test_mode_all_vs_any_on_straddling_cell():
    mesh = _quad_grid()
    box = [-0.1, -0.1, -0.1, 1.1, 1.1, 1.1]
    all_mode = crop(mesh, bbox=box, mode="all")
    any_mode = crop(mesh, bbox=box, mode="any")
    assert len(all_mode.cells[0].data) == 1
    assert len(any_mode.cells[0].data) == 2  # second quad straddles the boundary


def test_halfspace():
    mesh = _quad_grid()
    # keep x >= 1  (normal +x, point at x=1)
    out = crop(mesh, plane=([1, 0, 0], [1, 0, 0]), mode="all")
    assert len(out.cells[0].data) == 1  # only the right quad


def test_record_ids():
    mesh = _quad_grid()
    out = crop(
        mesh, bbox=[-0.1, -0.1, -0.1, 1.1, 1.1, 1.1], mode="all", record_ids=True
    )
    assert "crop:original_point_id" in out.point_data
    assert "crop:original_cell_id" in out.cell_data
    # original point ids are a subset of the input point indices
    assert set(out.point_data["crop:original_point_id"].tolist()).issubset(
        set(range(6))
    )


def test_pruned_mesh_has_no_unused_points():
    mesh = _quad_grid()
    out = crop(mesh, bbox=[-0.1, -0.1, -0.1, 1.1, 1.1, 1.1], mode="all")
    used = np.zeros(len(out.points), dtype=bool)
    for cb in out.cells:
        used[cb.data.reshape(-1)] = True
    assert used.all()


def test_cpp_matches_python():
    core = pytest.importorskip("meshioplusplus._core")
    mesh = _quad_grid()
    got = core.crop_bbox(mesh, [-0.1, -0.1, -0.1], [1.1, 1.1, 1.1], "any", False)
    # against the numpy fallback via the public API (which uses _core) — compare
    # point counts / cell counts as a structural check
    from meshioplusplus._crop import _crop_py, _points_3d

    p3 = _points_3d(mesh)
    mask = np.all((p3 >= [-0.1, -0.1, -0.1]) & (p3 <= [1.1, 1.1, 1.1]), axis=1)
    ref, _, _ = _crop_py(mesh, mask, "any", False)
    assert len(got["mesh"].points) == len(ref.points)
    assert len(got["mesh"].cells[0].data) == len(ref.cells[0].data)


def test_roundtrip_write_read(tmp_path):
    mesh = _quad_grid()
    out = crop(mesh, bbox=[-0.1, -0.1, -0.1, 1.1, 1.1, 1.1], mode="all")
    p = tmp_path / "cr.vtu"
    meshioplusplus.write(p, out)
    back = meshioplusplus.read(p)
    assert len(back.points) == len(out.points)


def test_requires_exactly_one_region():
    mesh = _quad_grid()
    with pytest.raises(ValueError):
        crop(mesh)
    with pytest.raises(ValueError):
        crop(mesh, bbox=[0, 0, 0, 1, 1, 1], plane=([0, 0, 0], [1, 0, 0]))
