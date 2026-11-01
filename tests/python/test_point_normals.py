"""Point-cloud normals by local plane fits (roadmap 7.5.2): radial and
outward on spheres, two pieces oriented apart, a plane facing +z or the
viewpoint, degenerate neighbourhoods, the PCD round trip, the pipeline step in
both engines."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "estimate_normals"), reason="needs the C++ core"
)


def _sphere_cloud(n, centre=(0.0, 0.0, 0.0), radius=1.0, seed=0):
    d = np.random.default_rng(seed).normal(size=(n, 3))
    return np.asarray(centre) + radius * d / np.linalg.norm(d, axis=1)[:, None]


def test_sphere_normals_are_radial_and_outward():
    pts = _sphere_cloud(4000, radius=3.0)
    out, report = mp.estimate_normals(
        pts, k=12, record_variation=True, return_report=True
    )
    assert report == {"num_degenerate": 0, "num_components": 1}
    n = out.point_data["normals"]
    np.testing.assert_allclose(np.linalg.norm(n, axis=1), 1.0, atol=1e-12)
    cos = np.einsum("ij,ij->i", n, pts / 3.0)
    assert cos.min() > 0.99
    var = out.point_data["normals:variation"]
    assert np.all((var >= 0) & (var < 0.02))


def test_torus_is_consistent():
    rng = np.random.default_rng(1)
    u, v = rng.uniform(0, 2 * np.pi, (2, 6000))
    R, r = 2.0, 0.6
    pts = np.column_stack(
        [
            (R + r * np.cos(v)) * np.cos(u),
            (R + r * np.cos(v)) * np.sin(u),
            r * np.sin(v),
        ]
    )
    exact = np.column_stack([np.cos(v) * np.cos(u), np.cos(v) * np.sin(u), np.sin(v)])
    n = mp.estimate_normals(pts, k=10).point_data["normals"]
    assert np.einsum("ij,ij->i", n, exact).min() > 0.95


def test_two_spheres_are_two_pieces_each_outward():
    a = _sphere_cloud(1000, seed=2)
    b = _sphere_cloud(1000, centre=(10, 0, 0), seed=3)
    out, report = mp.estimate_normals(np.vstack([a, b]), return_report=True)
    assert report["num_components"] == 2
    n = out.point_data["normals"]
    assert np.einsum("ij,ij->i", n[:1000], a).min() > 0.9
    assert np.einsum("ij,ij->i", n[1000:], b - [10, 0, 0]).min() > 0.9


def test_plane_faces_up_or_the_viewpoint():
    g = np.stack(np.meshgrid(np.arange(15.0), np.arange(15.0)), -1).reshape(-1, 2)
    pts = np.column_stack([g, np.zeros(len(g))])
    up = mp.estimate_normals(pts).point_data["normals"]
    np.testing.assert_array_equal(up, [[0.0, 0.0, 1.0]] * len(pts))
    down = mp.estimate_normals(pts, orientation="viewpoint", viewpoint=(7, 7, -3))
    np.testing.assert_array_equal(
        down.point_data["normals"], [[0.0, 0.0, -1.0]] * len(pts)
    )
    flat_2d = mp.estimate_normals(g).point_data["normals"]
    np.testing.assert_array_equal(flat_2d, up)


def test_mesh_passes_through_and_pcd_round_trip(tmp_path):
    sphere = mp.icosa_sphere(subdivisions=2)
    sphere.point_data["t"] = np.arange(len(sphere.points), dtype=float)
    out = mp.estimate_normals(sphere, k=8)
    np.testing.assert_array_equal(out.points, sphere.points)
    np.testing.assert_array_equal(out.cells[0].data, sphere.cells[0].data)
    np.testing.assert_array_equal(out.point_data["t"], sphere.point_data["t"])
    assert "normals" not in sphere.point_data
    path = str(tmp_path / "cloud.pcd")
    mp.write(
        path, mp.Mesh(out.points, [], point_data={"normals": out.point_data["normals"]})
    )
    back = mp.read(path)
    np.testing.assert_allclose(
        back.point_data["normals"], out.point_data["normals"], atol=1e-6
    )


def test_degenerate_and_refusals():
    line = np.column_stack([np.arange(8.0), 2 * np.arange(8.0), np.zeros(8)])
    out, report = mp.estimate_normals(line, return_report=True)
    assert report["num_degenerate"] == 8 and report["num_components"] == 0
    assert np.isnan(out.point_data["normals"]).all()
    with pytest.raises(ValueError, match="at least three points"):
        mp.estimate_normals(np.zeros((2, 3)))
    with pytest.raises(ValueError, match="k must be at least 2"):
        mp.estimate_normals(line, k=1)
    with pytest.raises(ValueError, match="unknown orientation"):
        mp.estimate_normals(line, orientation="inward")


def test_pipeline_step_in_both_engines(tmp_path):
    src = tmp_path / "cloud.vtu"
    pts = _sphere_cloud(500, seed=4)
    mp.write(str(src), mp.Mesh(pts, [("vertex", np.arange(500)[:, None])]))
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "EstimateNormals", "K": 10, "RecordVariation": True}],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep_py = mp.run_pipeline(doc)
    py = mp.read(str(tmp_path / "py.vtu"))
    assert "normals:variation" in py.point_data
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    rep_cpp = _core.run_pipeline_json(json.dumps(doc))
    cpp = mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(cpp.point_data["normals"], py.point_data["normals"])
    assert rep_py["steps"][0] == rep_cpp["steps"][0]
