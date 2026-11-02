"""Surface reconstruction from an oriented point cloud (roadmap 7.5.2): a
sphere within one sample spacing, closed and outward; Hoppe's field with
``k=1``; normals estimated or given as an array; an open scan left open; the
PLY round trip of a scan; the length unit; refusals; the pipeline step in both
engines."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core, _units

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "reconstruct_surface"), reason="needs the C++ core"
)


def _cloud(n, radius=1.0):
    return mp.sphere_points(n, radius).points


def _open_edges(mesh):
    tri = mesh.cells[0].data
    edges = np.sort(np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), 1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    return int((counts == 1).sum())


def _volume(mesh):
    p = mesh.points[mesh.cells[0].data]
    return float(np.einsum("ij,ij->i", p[:, 0], np.cross(p[:, 1], p[:, 2])).sum() / 6)


def test_sphere_within_one_spacing():
    pts = _cloud(4000)
    out, report = mp.reconstruct_surface(pts, normals=pts, return_report=True)
    assert [c.type for c in out.cells] == ["triangle"]
    assert report["num_samples"] == 4000 and not report["estimated_normals"]
    assert report["num_components"] == 1
    assert report["num_triangles"] == len(out.cells[0].data)
    assert _open_edges(out) == 0
    assert _volume(out) == pytest.approx(4 / 3 * np.pi, rel=0.02)  # outward
    h = mp.hausdorff_distance(out, mp.icosa_sphere(1.0, 5))
    assert h["distance"] < report["spacing"]
    assert not out.point_data and not out.cell_data


def test_one_neighbour_is_hoppe_and_mesh_input():
    pts = _cloud(3000, 2.0)
    cloud = mp.Mesh(
        pts, [("vertex", np.arange(len(pts))[:, None])], point_data={"n": pts / 2}
    )
    a = mp.reconstruct_surface(cloud, normals="n", k=1)
    assert _open_edges(a) == 0
    r = np.linalg.norm(a.points, axis=1)
    assert np.abs(r - 2.0).max() < 0.1
    b = mp.reconstruct_surface(cloud, normals="n")
    assert not np.array_equal(a.points, b.points)


def test_estimated_normals_and_refusal():
    pts = _cloud(3000)
    out, report = mp.reconstruct_surface(pts, return_report=True)
    assert report["estimated_normals"] and report["num_components"] == 1
    assert _volume(out) > 0
    with pytest.raises(ValueError, match="no 'normals' point data"):
        mp.reconstruct_surface(pts, estimate_normals=False)


def test_open_scan_stays_open():
    pts = _cloud(6000)
    half = pts[pts[:, 2] > 0]
    out, report = mp.reconstruct_surface(half, normals=half, return_report=True)
    assert report["num_undefined"] > 0 and report["num_components"] == 1
    assert _open_edges(out) > 0
    assert out.points[:, 2].min() > -(report["max_distance"] + report["cell_size"])


def test_ply_scan_round_trip_and_unit(tmp_path):
    pts = _cloud(2000)
    scan = mp.Mesh(pts, [], point_data={"normals": pts.copy()})
    _units.set_length_unit(scan, 1e-3)
    out = mp.reconstruct_surface(scan, resolution=24)
    assert _units.length_unit(out) == 1e-3
    path = str(tmp_path / "scan.ply")
    columns = {name: pts[:, i].copy() for i, name in enumerate(("nx", "ny", "nz"))}
    mp.write(path, mp.Mesh(pts, [], point_data=columns))
    back = mp.read(path)
    assert {"nx", "ny", "nz"} <= set(back.point_data)  # PLY's vertex properties
    again, report = mp.reconstruct_surface(back, resolution=24, return_report=True)
    assert not report["estimated_normals"]
    np.testing.assert_allclose(
        np.sort(np.linalg.norm(again.points, axis=1)),
        np.sort(np.linalg.norm(out.points, axis=1)),
        atol=1e-5,
    )


def test_refusals():
    pts = _cloud(200)
    with pytest.raises(ValueError, match="k must be at least 1"):
        mp.reconstruct_surface(pts, normals=pts, k=0)
    with pytest.raises(ValueError, match="not both"):
        mp.reconstruct_surface(pts, normals=pts, cell_size=0.1, resolution=10)
    with pytest.raises(ValueError, match="max_cells"):
        mp.reconstruct_surface(pts, normals=pts, cell_size=1e-3, max_cells=1000)
    with pytest.raises(ValueError, match="three coordinates"):
        mp.reconstruct_surface(pts[:, :2], normals=pts)
    with pytest.raises(ValueError, match=r"\(n, 3\) array"):
        mp.reconstruct_surface(pts, normals=pts[:10])
    with pytest.raises(ValueError, match="at least three samples"):
        mp.reconstruct_surface(pts, normals=np.zeros_like(pts))


def test_pipeline_step_in_both_engines(tmp_path):
    src = tmp_path / "cloud.vtu"
    pts = _cloud(800)
    mp.write(
        str(src),
        mp.Mesh(
            pts, [("vertex", np.arange(800)[:, None])], point_data={"normals": pts}
        ),
        compression=None,  # the build may have no zlib
    )
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "ReconstructSurface", "K": 4, "Resolution": 20}],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep_py = mp.run_pipeline(doc)
    py = mp.read(str(tmp_path / "py.vtu"))
    assert py.cells[0].type == "triangle"
    assert rep_py["steps"][0]["NumComponents"] == 1.0
    with pytest.raises(ValueError, match="K must be an integer"):
        mp.apply_operations(py, [{"Op": "ReconstructSurface", "K": 1.5}])
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    rep_cpp = _core.run_pipeline_json(json.dumps(doc))
    cpp = mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(cpp.points, py.points)
    np.testing.assert_array_equal(cpp.cells[0].data, py.cells[0].data)
    assert rep_py["steps"][0] == rep_cpp["steps"][0]
