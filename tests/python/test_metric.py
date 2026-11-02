"""The metric toolbox (roadmap 7.4.3): the intersection of two fields'
metrics resolves both fronts (the roadmap probe: it dominates each in every
direction), normalisation reaches its target complexity, clamping, gradation,
the refusals, the pipeline step in both engines."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "prepare_metric"), reason="needs the C++ core"
)


def _square(n=30, planar=False):
    sq = mp.convert_cells(mp.rectangle(dims=(n, n)), mode="simplexify")
    if planar:
        sq = mp.Mesh(sq.points[:, :2], sq.cells)
    return sq


def _full(m6):
    """(n, 6) lower-triangle metrics as (n, 3, 3)."""
    a, b, c, d, e, f = m6.T
    return np.stack(
        [np.stack([a, b, d], -1), np.stack([b, c, e], -1), np.stack([d, e, f], -1)], 1
    )


def _front_metrics():
    sq = _square()
    x, y = sq.points[:, 0], sq.points[:, 1]
    sq.point_data["u1"] = np.tanh(40 * (x - 0.5))
    sq.point_data["u2"] = np.tanh(40 * (y - 0.4))
    for name in ("u1", "u2"):
        sq = mp.hessian(sq, name, location="point")
        sq = mp.metric_from_hessian(
            sq, f"{name}:hessian", error=1e-2, hmin=1e-3, hmax=0.5, output=f"m_{name}"
        )
    return sq


def test_intersection_resolves_both_fronts():
    """The probe: the intersection asks, in every direction, for at least what
    each field's metric does, so it is fine across both fronts."""
    sq = _front_metrics()
    out = mp.prepare_metric(sq, ["m_u1", "m_u2"], output="both")
    m = _full(out.point_data["both"])
    a = _full(sq.point_data["m_u1"])
    b = _full(sq.point_data["m_u2"])
    rng = np.random.default_rng(0)
    e = rng.normal(size=(64, 3))
    e /= np.linalg.norm(e, axis=1)[:, None]

    def q(t):
        return np.einsum("ki,nij,kj->nk", e, t, e)

    assert np.all(q(m) >= np.maximum(q(a), q(b)) * (1 - 1e-9) - 1e-9)
    x, y = sq.points[:, 0], sq.points[:, 1]
    # (on a front's line itself the curvature vanishes: the fine sizes sit
    # just beside it, so compare the bands around each front with far away)
    band1, far1 = np.abs(x - 0.5) < 0.05, np.abs(x - 0.5) > 0.2
    band2, far2 = np.abs(y - 0.4) < 0.05, np.abs(y - 0.4) > 0.2
    np.testing.assert_allclose(m[band1 & far2, 0, 0], a[band1 & far2, 0, 0], rtol=1e-6)
    np.testing.assert_allclose(m[band2 & far1, 1, 1], b[band2 & far1, 1, 1], rtol=1e-6)
    assert m[band1, 0, 0].max() > 50 * m[far1 & far2, 0, 0].max()
    assert m[band2, 1, 1].max() > 50 * m[far1 & far2, 1, 1].max()


def test_normalisation_reaches_the_target():
    sq = _square(10, planar=True)
    sq.point_data["m"] = np.tile([400.0, 0.0, 25.0], (len(sq.points), 1))
    for norm in (0.0, 1.0, 2.0):
        out, rep = mp.prepare_metric(
            sq, "m", target_vertices=3000, norm=norm, return_report=True
        )
        assert rep["complexity_before"] == pytest.approx(100.0)
        assert rep["complexity_after"] == pytest.approx(3000.0)
        m = out.point_data["metric"]
        np.testing.assert_allclose(m[:, 0] / m[:, 2], 16.0)


def test_a_3d_metric_on_a_planar_mesh_normalises_in_its_plane():
    sq = _square(10)  # 3-D points with z = 0, as hessian makes 3-D metrics there
    sq.point_data["m"] = np.tile([400.0, 0.0, 25.0, 0.0, 0.0, 4.0], (len(sq.points), 1))
    out, rep = mp.prepare_metric(sq, "m", target_vertices=2000, return_report=True)
    assert rep["complexity_before"] == pytest.approx(100.0)
    assert rep["complexity_after"] == pytest.approx(2000.0)
    bent = mp.Mesh(
        sq.points + np.column_stack([0 * sq.points[:, :2], sq.points[:, 0] ** 2]),
        sq.cells,
    )
    bent.point_data["m"] = sq.point_data["m"]
    with pytest.raises(ValueError, match="cells of the metric's dimension"):
        mp.prepare_metric(bent, "m", target_vertices=2000)


def test_clamp_and_gradation():
    sq = _square(20, planar=True)
    m = np.tile([1.0, 0.0, 1.0], (len(sq.points), 1))
    centre = int(np.argmin(np.linalg.norm(sq.points - 0.5, axis=1)))
    m[centre] = [1e6, 0.0, 1e6]
    sq.point_data["m"] = m
    out = mp.prepare_metric(sq, "m", hmin=0.01, hmax=0.5)
    size = 1 / np.sqrt(out.point_data["metric"][:, 0])
    assert size.min() == pytest.approx(0.01) and size.max() == pytest.approx(0.5)
    graded, rep = mp.prepare_metric(
        sq, "m", hmin=0.01, hmax=0.5, gradation=1.3, return_report=True
    )
    assert rep["num_graded"] > 0
    gsize = 1 / np.sqrt(graded.point_data["metric"][:, 0])
    dist = np.linalg.norm(sq.points - sq.points[centre], axis=1)
    # sizes grow away from the fine point, and gradation only ever refines
    assert np.all(gsize <= size + 1e-12)
    near = (dist > 0) & (dist < 0.1)
    assert gsize[near].max() < 0.5
    order = np.argsort(dist)
    assert np.corrcoef(dist[order][:200], gsize[order][:200])[0, 1] > 0.9


def test_refusals():
    sq = _square(4, planar=True)
    sq.point_data["m"] = np.tile([1.0, 0.0, 1.0], (len(sq.points), 1))
    sq.point_data["w"] = np.zeros((len(sq.points), 6))
    with pytest.raises(ValueError, match="at least one"):
        mp.prepare_metric(sq, [])
    with pytest.raises(ValueError, match="no point data"):
        mp.prepare_metric(sq, "missing")
    with pytest.raises(ValueError, match="like the others"):
        mp.prepare_metric(sq, ["m", "w"])
    with pytest.raises(ValueError, match="norm"):
        mp.prepare_metric(sq, "m", norm=0.5)
    with pytest.raises(ValueError, match="gradation"):
        mp.prepare_metric(sq, "m", gradation=0.9)


def test_pipeline_step_in_both_engines(tmp_path):
    sq = _square(12, planar=True)
    x = sq.points[:, 0]
    sq.point_data["a"] = np.column_stack(
        [1 + 400 * np.exp(-((x - 0.5) ** 2) / 0.01), np.zeros(len(x)), np.ones(len(x))]
    )
    sq.point_data["b"] = sq.point_data["a"][:, [2, 1, 0]]
    src = tmp_path / "s.vtu"
    mp.write(str(src), sq, compression=None)
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {
                "Op": "PrepareMetric",
                "Arrays": ["a", "b"],
                "TargetVertices": 800,
                "Hmin": 0.005,
                "Hmax": 0.3,
                "Gradation": 1.4,
            }
        ],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep_py = mp.run_pipeline(doc)
    py = mp.read(str(tmp_path / "py.vtu"))
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    rep_cpp = _core.run_pipeline_json(json.dumps(doc))
    cpp = mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(cpp.point_data["metric"], py.point_data["metric"])
    assert rep_py["steps"][0] == rep_cpp["steps"][0]
