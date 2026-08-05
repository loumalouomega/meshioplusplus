"""Tests for the crop operation (bounding box / half-space / data predicate)."""

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


# --------------------------------------------------------------------------- #
# the predicate crop                                                           #
# --------------------------------------------------------------------------- #
def _tagged_quad_grid(with_nan=False):
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
            [3, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [2, 1, 0],
            [3, 1, 0],
        ],
        float,
    )
    m = meshioplusplus.Mesh(
        pts, [("quad", np.array([[0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6]]))]
    )
    m.cell_data["t"] = [np.array([0.0, np.nan if with_nan else 1.0, 2.0])]
    return m


def test_predicate_keeps_the_matching_cells():
    out = crop(_tagged_quad_grid(), where=("t", "<", 1.5))
    assert len(out.cells[0].data) == 2
    assert len(out.points) == 6  # the dropped cell's own points are pruned
    assert np.array_equal(out.cell_data["t"][0], [0.0, 1.0])


@pytest.mark.parametrize(
    "compare, n", [("<", 1), ("<=", 2), (">", 1), (">=", 2), ("==", 1), ("!=", 2)]
)
def test_predicate_honours_every_comparison(compare, n):
    assert len(crop(_tagged_quad_grid(), where=("t", compare, 1.0)).cells[0].data) == n


@pytest.mark.parametrize("compare", ["<", "<=", ">", ">=", "==", "!="])
def test_a_non_finite_cell_value_never_matches(compare):
    """The rule ``refine`` states and this inherits by sharing the evaluator.

    It has to hold for ``!=`` too, which is the case a naive implementation gets
    wrong: ``NaN != 1.0`` is true in IEEE.
    """
    out = crop(_tagged_quad_grid(with_nan=True), where=("t", compare, 1.0))
    assert np.all(np.isfinite(out.cell_data["t"][0]))


def test_predicate_rejects_what_is_not_a_scalar_cell_array():
    m = _tagged_quad_grid()
    with pytest.raises(ValueError, match="no cell_data array"):
        crop(m, where=("nope", "<", 1.0))
    # point_data is refused BY NAME rather than averaged onto the cells.
    m.point_data["p"] = np.zeros(len(m.points))
    with pytest.raises(ValueError, match="point_data"):
        crop(m, where=("p", "<", 1.0))
    m.cell_data["v"] = [np.zeros((3, 3))]
    with pytest.raises(ValueError, match="scalar"):
        crop(m, where=("v", "<", 1.0))


def test_predicate_mode_is_refused_rather_than_ignored():
    with pytest.raises(ValueError, match="mode="):
        crop(_tagged_quad_grid(), where=("t", "<", 1.0), mode="any")
    with pytest.raises(ValueError, match="unknown comparison"):
        crop(_tagged_quad_grid(), where=("t", "~", 1.0))
    with pytest.raises(ValueError, match="exactly one"):
        crop(_tagged_quad_grid(), bbox=[0, 0, 0, 1, 1, 1], where=("t", "<", 1.0))


def test_predicate_cpp_matches_python():
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._crop import _crop_predicate_py

    m = _tagged_quad_grid(with_nan=True)
    for compare in ("<", "<=", ">", ">=", "==", "!="):
        got = core.crop_predicate(m, "t", compare, 1.0, True)
        want, _, _ = _crop_predicate_py(m, "t", compare, 1.0, True)
        a, b = np.asarray(got["mesh"].points), np.asarray(want.points)
        assert a.dtype == b.dtype and a.tobytes() == b.tobytes()
        ca = np.asarray(got["mesh"].cells[0].data)
        cb = np.asarray(want.cells[0].data)
        assert ca.dtype == cb.dtype and ca.tobytes() == cb.tobytes()
        assert np.array_equal(
            got["mesh"].cell_data["crop:original_cell_id"][0],
            want.cell_data["crop:original_cell_id"][0],
        )


def test_inside_outside_composes_from_a_distance_field():
    """The composition the general predicate exists to serve.

    A dedicated crop-by-surface would have served this one case; this serves it
    and every other field a mesh carries.
    """
    box = meshioplusplus.extract_surface(
        meshioplusplus.grid([2, 2, 2], (0.25, 0.25, 0.25), (0.25, 0.25, 0.25))
    )
    domain = meshioplusplus.grid([6, 6, 6], (0, 0, 0), (0.2, 0.2, 0.2))
    field = meshioplusplus.distance_to_surface(domain, box, location="center")
    inside = crop(field, where=("sdf:distance", "<", 0.0))
    outside = crop(field, where=("sdf:distance", ">=", 0.0))
    n_in = len(inside.cells[0].data)
    n_out = len(outside.cells[0].data)
    assert n_in > 0 and n_out > 0
    # The two halves partition the domain: every cell is on exactly one side.
    assert n_in + n_out == len(domain.cells[0].data)
