"""Convex hull and planar Delaunay triangulation of scattered points (roadmap
7.5.1): extreme points only, closed outward surfaces, empty circumcircles, the
terrain keeping its z, the refusals, the pipeline step in both engines."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "convex_hull"), reason="needs the C++ core"
)


def _closed_and_outward(hull):
    tris = hull.cells_dict["triangle"]
    directed = {
        (a, b) for t in tris for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))
    }
    assert len(directed) == 3 * len(tris)
    assert all((b, a) in directed for a, b in directed)
    p = hull.points
    centre = p.mean(axis=0)
    n = np.cross(p[tris[:, 1]] - p[tris[:, 0]], p[tris[:, 2]] - p[tris[:, 0]])
    assert np.all(np.einsum("ij,ij->i", n, p[tris[:, 0]] - centre) > 0)


def test_hull_of_a_cloud_contains_it():
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(3000, 3))
    hull, report = mp.convex_hull(pts, return_report=True)
    _closed_and_outward(hull)
    assert report["dimension"] == 3
    assert report["num_facets"] == 2 * report["num_vertices"] - 4
    ids = hull.point_data["hull:point_id"]
    np.testing.assert_array_equal(hull.points, pts[ids])
    assert np.all(np.diff(ids) > 0)
    tris = hull.cells_dict["triangle"]
    p = hull.points
    n = np.cross(p[tris[:, 1]] - p[tris[:, 0]], p[tris[:, 2]] - p[tris[:, 0]])
    height = np.einsum("tk,tnk->tn", n, pts[None, :, :] - p[tris[:, 0]][:, None, :])
    assert height.max() <= 1e-12
    vol = np.einsum("ij,ij->i", np.cross(p[tris[:, 1]], p[tris[:, 2]]), p[tris[:, 0]])
    assert report["volume"] == pytest.approx(vol.sum() / 6, rel=1e-12)


def test_lattice_keeps_only_the_corners():
    g = np.stack(np.meshgrid(*[np.arange(5.0)] * 3, indexing="ij"), -1).reshape(-1, 3)
    mesh = mp.Mesh(g, [])
    mesh.point_data["id"] = np.arange(len(g))
    hull, report = mp.convex_hull(mesh, return_report=True)
    assert report["num_vertices"] == 8 and report["num_facets"] == 12
    assert report["volume"] == 64.0 and report["area"] == 96.0
    np.testing.assert_array_equal(
        hull.point_data["id"], hull.point_data["hull:point_id"]
    )
    _closed_and_outward(hull)


def test_planar_hull_is_a_counter_clockwise_loop():
    rng = np.random.default_rng(1)
    pts = np.column_stack([rng.random((500, 2)), np.full(500, 2.5)])
    hull, report = mp.convex_hull(pts, return_report=True)
    assert report["dimension"] == 2
    lines = hull.cells_dict["line"]
    assert len(lines) == report["num_vertices"] == report["num_facets"]
    np.testing.assert_array_equal(lines[:, 1], np.roll(lines[:, 0], -1))
    loop = hull.points[lines[:, 0], :2]
    x, y = loop.T
    area = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
    assert area == pytest.approx(report["area"]) and area > 0
    assert np.all(hull.points[:, 2] == 2.5)
    # dimension=2 on a 3-D cloud takes the x, y footprint
    cloud = rng.normal(size=(400, 3))
    flat = mp.convex_hull(cloud, dimension=2, return_report=True)[1]
    assert flat["dimension"] == 2 and flat["volume"] == 0.0


def test_delaunay_circumcircles_are_empty_and_z_rides_along():
    rng = np.random.default_rng(2)
    xy = rng.random((400, 2))
    z = np.sin(4 * xy[:, 0]) * xy[:, 1]
    terrain = mp.Mesh(np.column_stack([xy, z]), [("vertex", np.arange(400)[:, None])])
    terrain.point_data["z"] = z
    out, report = mp.delaunay_2d(terrain, return_report=True)
    assert [b.type for b in out.cells] == ["triangle"]
    np.testing.assert_array_equal(out.points, terrain.points)
    np.testing.assert_array_equal(out.point_data["z"], z)
    tris = out.cells_dict["triangle"]
    assert report["num_triangles"] == len(tris) == 2 * 400 - report["num_boundary"] - 2
    a, b, c = (xy[tris[:, k]] for k in range(3))

    def cross(u, v):
        return u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]

    assert np.all(cross(b - a, c - a) > 0)  # counter-clockwise
    # circumcentres; no point strictly inside any circumcircle
    d = 2 * cross(b - a, c - a)
    ux = (
        (a**2).sum(1) * (b[:, 1] - c[:, 1])
        + (b**2).sum(1) * (c[:, 1] - a[:, 1])
        + (c**2).sum(1) * (a[:, 1] - b[:, 1])
    ) / d
    uy = (
        (a**2).sum(1) * (c[:, 0] - b[:, 0])
        + (b**2).sum(1) * (a[:, 0] - c[:, 0])
        + (c**2).sum(1) * (b[:, 0] - a[:, 0])
    ) / d
    centre = np.column_stack([ux, uy])
    r2 = ((a - centre) ** 2).sum(1)
    dist2 = ((xy[None, :, :] - centre[:, None, :]) ** 2).sum(-1)
    assert np.all(dist2 >= r2[:, None] * (1 - 1e-9))


def test_delaunay_duplicates_and_grid():
    g = np.stack(np.meshgrid(np.arange(4.0), np.arange(4.0)), -1).reshape(-1, 2)
    pts = np.vstack([g, g[5:6]])
    out, report = mp.delaunay_2d(pts, return_report=True)
    assert report == {"num_triangles": 18, "num_duplicates": 1, "num_boundary": 12}
    assert 16 not in out.cells_dict["triangle"]
    again = mp.delaunay_2d(pts)
    np.testing.assert_array_equal(again.cells[0].data, out.cells[0].data)


def test_refusals():
    with pytest.raises(ValueError, match="coplanar"):
        mp.convex_hull(
            np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0.0]]), dimension=3
        )
    with pytest.raises(ValueError, match="collinear"):
        mp.convex_hull(np.array([[0, 0], [1, 1], [2, 2.0]]))
    with pytest.raises(ValueError, match="dimension"):
        mp.convex_hull(np.eye(3), dimension=4)
    with pytest.raises(ValueError, match="non-finite"):
        mp.convex_hull(np.array([[0, 0, 0], [1, 0, 0], [0, np.nan, 1.0]]))
    with pytest.raises(ValueError, match="collinear"):
        mp.delaunay_2d(np.array([[0, 0], [1, 1], [2, 2.0]]))
    with pytest.raises(ValueError, match="exact range"):
        mp.delaunay_2d(np.array([[0, 0], [1e300, 0], [0, 1.0]]))
    with pytest.raises(ValueError, match="array of points"):
        mp.delaunay_2d(np.zeros(5))


@pytest.mark.parametrize(
    "step",
    [
        {"Op": "ConvexHull"},
        {"Op": "ConvexHull", "Dimension": 2},
        {"Op": "Delaunay2D"},
    ],
)
def test_pipeline_step_in_both_engines(tmp_path, step):
    rng = np.random.default_rng(3)
    src = tmp_path / "cloud.vtu"
    pts = rng.random((300, 3))
    mp.write(str(src), mp.Mesh(pts, [("vertex", np.arange(300)[:, None])]))
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [step],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep_py = mp.run_pipeline(doc)
    py = mp.read(str(tmp_path / "py.vtu"))
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    rep_cpp = _core.run_pipeline_json(json.dumps(doc))
    cpp = mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(cpp.points, py.points)
    np.testing.assert_array_equal(cpp.cells[0].data, py.cells[0].data)
    assert rep_py["steps"][0] == rep_cpp["steps"][0]
