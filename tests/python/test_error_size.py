"""Tests for ``meshioplusplus.error_size`` (roadmap 6.10.2): a per-point size
field from a per-cell error indicator by the Zienkiewicz-Zhu adaptive
procedure.

The C++ core (``_core.error_size``) is checked against the numpy twin
``_error_size_py`` (they agree to rounding), the refusals by name, and the
adaptive loop end to end: ``estimate_error`` -> ``error_size`` ->
``mesh_cad_surface(background=...)`` on the CAD JSON sphere refines where the
field is badly resolved.
"""

import math
import pathlib

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core
from meshioplusplus._error_size import _error_size_py

_SPHERE = pathlib.Path(__file__).parent / "meshes" / "cad_json" / "sphere.cad.json"

needs_core = pytest.mark.skipif(
    not hasattr(_core, "error_size"), reason="core predates error_size"
)


def _indicated(n=4):
    """A hexahedron grid carrying the ZZ indicator of a quadratic field."""
    m = mp.box((n, n, n), (n, n, n))
    m.point_data["f"] = (m.points**2).sum(axis=1)
    return mp.estimate_error(m, "f")


def _mixed():
    """Triangles and quads in two blocks, with a hand-made indicator; one cell
    reads NaN and one 0 (an infinite size), both skipped."""
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0], [3, 0, 0]],
        dtype=float,
    )
    tri = np.array([[0, 1, 4], [0, 4, 3], [2, 6, 5]])
    quad = np.array([[1, 2, 5, 4]])
    m = mp.Mesh(pts, [("triangle", tri), ("quad", quad)])
    m.cell_data["error:zz"] = [np.array([0.5, np.nan, 0.0]), np.array([2.0])]
    return m


def _assert_twin(mesh, **kw):
    args = dict(
        target_error=kw.get("target_error", 1e-2),
        order=kw.get("order", 1.0),
        min_size=kw.get("min_size", 0.0),
        max_size=kw.get("max_size", 0.0),
    )
    out, report = mp.error_size(mesh, return_report=True, **args)
    twin, twin_report = _error_size_py(
        mesh,
        "error:zz",
        args["target_error"],
        args["order"],
        args["min_size"],
        args["max_size"],
        "size",
    )
    np.testing.assert_allclose(
        out.point_data["size"], twin.point_data["size"], rtol=1e-12, equal_nan=True
    )
    for key in ("num_cells", "num_skipped", "num_clamped"):
        assert report[key] == twin_report[key], key
    for key in ("min_size", "max_size"):
        assert report[key] == pytest.approx(twin_report[key], rel=1e-12), key
    return out, report


# --------------------------------------------------------------------------- #
# the core and its numpy twin                                                 #
# --------------------------------------------------------------------------- #


@needs_core
def test_core_matches_twin_on_hexahedra():
    out, report = _assert_twin(_indicated())
    # The quadratic field's indicator vanishes on a cell or so by symmetry: an
    # infinite size, skipped.
    assert report["num_cells"] + report["num_skipped"] == 64
    assert report["num_cells"] > 56
    assert np.all(np.isfinite(out.point_data["size"]))


@needs_core
@pytest.mark.parametrize("order", [1.0, 2.0, 0.5])
def test_core_matches_twin_by_order(order):
    _assert_twin(_indicated(3), order=order, target_error=0.1)


@needs_core
def test_core_matches_twin_with_bounds():
    mesh = _indicated(3)
    _, free = _assert_twin(mesh)
    lo = 0.5 * (free["min_size"] + free["max_size"])
    out, report = _assert_twin(mesh, min_size=lo, max_size=2.0 * lo)
    assert report["num_clamped"] > 0
    s = out.point_data["size"]
    assert np.all((s >= lo) & (s <= 2.0 * lo))


@needs_core
def test_core_matches_twin_on_mixed_blocks_with_skipped_cells():
    out, report = _assert_twin(_mixed())
    assert report["num_skipped"] == 2 and report["num_cells"] == 2
    s = out.point_data["size"]
    # Point 3 is only on the NaN triangle and point 6 only on the zero one.
    assert np.isnan(s[3]) and np.isnan(s[6])
    assert np.all(np.isfinite(s[[0, 1, 2, 4, 5]]))
    # A max_size gives the zero-indicator cell a size: no longer skipped.
    _, bounded = _assert_twin(_mixed(), max_size=10.0)
    assert bounded["num_skipped"] == 1 and bounded["num_clamped"] >= 1


def test_formula_on_one_triangle():
    """One cell: h (eta / sqrt(1) / eta_K)^(1/p), h its longest side."""
    m = mp.Mesh(
        np.array([[0, 0, 0], [3, 0, 0], [0, 4, 0]], float), [("triangle", [[0, 1, 2]])]
    )
    m.cell_data["error:zz"] = [np.array([0.2])]
    out, report = mp.error_size(m, target_error=0.05, order=2.0, return_report=True)
    expect = 5.0 * math.sqrt(0.05 / 0.2)
    np.testing.assert_allclose(out.point_data["size"], expect, rtol=1e-14)
    assert report["min_size"] == pytest.approx(expect)
    assert report["num_cells"] == 1


def test_passes_everything_else_through():
    mesh = _indicated(2)
    out = mp.error_size(mesh, target_error=1e-2, output="h")
    assert np.array_equal(out.points, mesh.points)
    assert set(out.point_data) == set(mesh.point_data) | {"h"}
    assert set(out.cell_data) == set(mesh.cell_data)
    assert "h" not in mesh.point_data


# --------------------------------------------------------------------------- #
# refusals                                                                     #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("target", [0.0, -1.0, math.inf, math.nan])
def test_refuses_a_bad_target_error(target):
    with pytest.raises(ValueError, match="target_error"):
        mp.error_size(_mixed(), target_error=target)


@pytest.mark.parametrize("order", [0.0, -1.0, math.nan])
def test_refuses_a_bad_order(order):
    with pytest.raises(ValueError, match="order"):
        mp.error_size(_mixed(), target_error=1.0, order=order)


def test_refuses_bad_bounds():
    with pytest.raises(ValueError, match="not negative"):
        mp.error_size(_mixed(), target_error=1.0, min_size=-1.0)
    with pytest.raises(ValueError, match="exceed"):
        mp.error_size(_mixed(), target_error=1.0, min_size=2.0, max_size=1.0)


def test_refuses_a_missing_or_wrong_indicator():
    m = _mixed()
    with pytest.raises(ValueError, match="nope"):
        mp.error_size(m, "nope", target_error=1.0)
    m.point_data["p"] = np.zeros(len(m.points))
    with pytest.raises(ValueError, match="point_data"):
        mp.error_size(m, "p", target_error=1.0)
    m.cell_data["v"] = [np.zeros((3, 2)), np.zeros((1, 2))]
    with pytest.raises(ValueError, match="one component"):
        mp.error_size(m, "v", target_error=1.0)


def test_refuses_an_output_collision_without_overwrite():
    m = _mixed()
    m.point_data["size"] = np.zeros(len(m.points))
    with pytest.raises(ValueError, match="overwrite"):
        mp.error_size(m, target_error=1.0)
    out = mp.error_size(m, target_error=1.0, overwrite=True)
    assert np.isfinite(out.point_data["size"][0])


# --------------------------------------------------------------------------- #
# the adaptive loop end to end                                                #
# --------------------------------------------------------------------------- #


@needs_core
def test_adaptive_loop_on_the_sphere():
    """A field steep near the north pole: the indicator is large there, so
    error_size asks for small cells there, and mesh_cad_surface with that
    background puts more, smaller triangles in the north than in the south."""
    if not hasattr(_core, "mesh_cad_surface"):
        pytest.skip("core predates mesh_cad_surface")
    model = mp.read(_SPHERE)
    coarse = mp.mesh_cad_surface(model, 0.3, edges=False)
    coarse.point_data["f"] = np.exp(4.0 * coarse.points[:, 2])
    indicated = mp.estimate_error(coarse, "f")
    sized, report = mp.error_size(
        indicated, target_error=0.2, min_size=0.04, max_size=0.4, return_report=True
    )
    assert report["num_cells"] > 0 and report["num_skipped"] == 0
    size = sized.point_data["size"]
    north = coarse.points[:, 2] > 0.5
    south = coarse.points[:, 2] < -0.5
    assert size[north].mean() < size[south].mean()

    adapted, rep = mp.mesh_cad_surface(
        model, 0.4, edges=False, background=sized, return_report=True
    )
    assert rep["num_failed"] == 0
    tri = np.concatenate(
        [np.asarray(b.data) for b in adapted.cells if b.type == "triangle"]
    )
    centroid_z = adapted.points[tri].mean(axis=1)[:, 2]
    assert (centroid_z > 0.5).sum() > 2 * (centroid_z < -0.5).sum()
