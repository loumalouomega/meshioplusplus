"""Tests for ODT remeshing (the ``optimize_volume`` operation)."""

import hashlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import optimize_volume

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(
    _core is None, reason="optimize_volume has no pure-Python fallback"
)


def _two_tet_23():
    # (a,b,c,d) and (a,b,c,e) share face abc; the 2-3 flip around edge de raises
    # the worst scaled Jacobian from ~0.066 to ~0.362.
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.3, 0.9, 0.0],
            [0.2252, 0.2808, 0.6977],
            [0.4659, 0.3149, -0.0367],
        ],
        dtype=float,
    )
    cells = [("tetra", np.array([[0, 1, 2, 3], [0, 1, 2, 4]], dtype=np.int64))]
    return meshioplusplus.Mesh(pts, cells)


def _three_tet_32():
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.3, 0.9, 0.0],
            [0.45, 0.35, 0.9],
            [0.45, 0.35, -0.9],
        ],
        dtype=float,
    )
    cells = [
        ("tetra", np.array([[3, 4, 0, 1], [3, 4, 1, 2], [3, 4, 2, 0]], dtype=np.int64))
    ]
    return meshioplusplus.Mesh(pts, cells)


def _jittered_grid(seed=1):
    grid = meshioplusplus.grid([4, 4, 4], origin=[0, 0, 0], spacing=[1, 1, 1])
    tet = meshioplusplus.convert_cells(grid, mode="simplexify")
    pts = tet.points.copy()
    onb = (
        np.isclose(pts[:, 0], 0)
        | np.isclose(pts[:, 0], 4)
        | np.isclose(pts[:, 1], 0)
        | np.isclose(pts[:, 1], 4)
        | np.isclose(pts[:, 2], 0)
        | np.isclose(pts[:, 2], 4)
    )
    rng = np.random.default_rng(seed)
    pts[~onb] += rng.uniform(-0.35, 0.35, size=(np.sum(~onb), 3))
    mesh = meshioplusplus.Mesh(pts, [("tetra", tet.cells[0].data.copy())])
    mesh.point_data["temp"] = pts[:, 0].copy()
    return mesh


def _boundary_facets(mesh):
    surf = meshioplusplus.extract_surface(mesh)
    tris = surf.cells[0].data
    return {tuple(sorted(int(v) for v in t)) for t in tris}


@needs_core
def test_2_3_flip_fires_and_improves():
    _, rep = optimize_volume(_two_tet_23(), relocate=False, return_report=True)
    assert rep["num_23_flips"] == 1
    assert rep["num_32_flips"] == 0
    assert rep["num_tets"] == 3
    assert rep["min_quality_after"] > rep["min_quality_before"] + 0.2


@needs_core
def test_3_2_flip_fires_and_improves():
    _, rep = optimize_volume(_three_tet_32(), relocate=False, return_report=True)
    assert rep["num_32_flips"] == 1
    assert rep["num_23_flips"] == 0
    assert rep["num_tets"] == 2
    assert rep["min_quality_after"] > rep["min_quality_before"] + 0.2


@needs_core
def test_flips_help_where_smoothing_cannot():
    # every vertex of the fixture is on the boundary, so relocation cannot move
    # a single point -- any improvement is proof the connectivity changed.
    _, sm = optimize_volume(_two_tet_23(), flip=False, return_report=True)
    assert sm["num_vertices_moved"] == 0
    assert sm["min_quality_after"] == pytest.approx(sm["min_quality_before"])
    _, full = optimize_volume(_two_tet_23(), return_report=True)
    assert full["num_flips"] > 0
    assert full["min_quality_after"] > sm["min_quality_after"] + 0.2


@needs_core
def test_monotone_quality_no_inversion_and_boundary_invariant():
    mesh = _jittered_grid()
    before = _boundary_facets(mesh)
    out, rep = optimize_volume(mesh, return_report=True)
    assert rep["min_quality_after"] >= rep["min_quality_before"] - 1e-12
    assert _core.compute_stats(out)["num_inverted"] == 0
    assert _boundary_facets(out) == before


@needs_core
def test_point_set_and_point_data_carry():
    mesh = _jittered_grid()
    out = optimize_volume(mesh)
    assert len(out.points) == len(mesh.points)
    assert "temp" in out.point_data
    np.testing.assert_array_equal(out.point_data["temp"], mesh.point_data["temp"])


@needs_core
def test_point_regions_carry_cell_regions_dropped():
    mesh = _jittered_grid()
    mesh.point_sets["corner"] = np.array([0, 1, 2], dtype=np.int64)
    mesh.cell_sets["all"] = [np.arange(len(mesh.cells[0].data), dtype=np.int64)]
    out = optimize_volume(mesh)
    assert "corner" in out.point_sets
    np.testing.assert_array_equal(np.sort(out.point_sets["corner"]), [0, 1, 2])
    assert "all" not in out.cell_sets


@needs_core
def test_deterministic_across_runs():
    def h(m):
        d = hashlib.sha256()
        d.update(m.points.tobytes())
        d.update(m.cells[0].data.tobytes())
        return d.hexdigest()

    a = optimize_volume(_jittered_grid())
    b = optimize_volume(_jittered_grid())
    assert h(a) == h(b)


@needs_core
def test_rejects_non_tet_by_name():
    hexm = meshioplusplus.grid([1, 1, 1], origin=[0, 0, 0], spacing=[1, 1, 1])
    with pytest.raises(ValueError, match="tetra"):
        optimize_volume(hexm)
    tri = meshioplusplus.Mesh(
        np.zeros((3, 3)), [("triangle", np.array([[0, 1, 2]], dtype=np.int64))]
    )
    with pytest.raises(ValueError):
        optimize_volume(tri)


@needs_core
def test_no_op_when_both_halves_disabled():
    _, rep = optimize_volume(
        _two_tet_23(), relocate=False, flip=False, return_report=True
    )
    assert rep["num_flips"] == 0
    assert rep["num_vertices_moved"] == 0
    assert rep["num_tets"] == 2
