"""Point-cloud conditioning (roadmap 9.9): the Python surface of
``filter_points``, its CLI verb and its pipeline step."""

import json

import numpy as np
import pytest

import meshioplusplus as mio

_core = pytest.importorskip("meshioplusplus._core")
pytestmark = pytest.mark.skipif(
    not hasattr(_core, "filter_points"), reason="needs the compiled core"
)


def _scan(n_in=1500, n_out=12, seed=5):
    rng = np.random.default_rng(seed)
    inl = np.column_stack(
        [rng.uniform(0, 1, n_in), rng.uniform(0, 1, n_in), rng.normal(0, 0.002, n_in)]
    )
    out = rng.uniform(0, 1, (n_out, 3)) + [3.0, 3.0, 2.0]
    return np.vstack([inl, out])


def test_statistical_removes_the_far_points():
    pts = _scan()
    mesh, rep = mio.filter_points(
        pts, statistical_k=12, std_ratio=2.0, return_report=True
    )
    assert rep["num_input"] == len(pts)
    assert rep["num_removed_statistical"] >= 12
    assert (rep["point_map"] < 1500).all()
    assert len(mesh.points) == len(rep["point_map"])
    assert np.array_equal(np.asarray(mesh.points), pts[rep["point_map"]])


def test_radius_voxel_and_poisson():
    pts = _scan(3000, 8)
    r = mio.filter_points(pts, radius=0.08, min_neighbors=5, return_report=True)[1]
    assert r["num_removed_radius"] >= 8
    v = mio.filter_points(pts[:3000], voxel_size=0.1, return_report=True)[1]
    assert 0 < v["num_removed_voxel"] < 3000
    mesh, p = mio.filter_points(pts[:3000], min_distance=0.05, return_report=True)
    xyz = np.asarray(mesh.points)
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    np.fill_diagonal(d, np.inf)
    assert d.min() >= 0.05
    assert p["point_map"][0] == 0  # greedy in index order


def test_data_and_point_regions_follow():
    pts = np.array([[i, 0.0, 0.0] for i in range(6)] + [[100.0, 100.0, 100.0]])
    mesh = mio.Mesh(pts, [("vertex", np.arange(7).reshape(-1, 1))])
    mesh.point_data["id"] = np.arange(7, dtype=float)
    mesh.regions.append(mio.Region("keep", "point", [1, 6, 3]))
    out, rep = mio.filter_points(mesh, radius=1.5, min_neighbors=1, return_report=True)
    assert len(out.points) == 6
    assert out.point_data["id"].tolist() == [0, 1, 2, 3, 4, 5]
    assert [r.name for r in out.regions] == ["keep"]
    assert out.regions[0].entries.tolist() == [1, 3]


def test_refusals():
    pts = _scan(100, 0)
    with pytest.raises(ValueError):
        mio.filter_points(pts)  # no step
    with pytest.raises(ValueError):
        mio.filter_points(pts, radius=-1.0)
    with pytest.raises(ValueError):
        mio.filter_points(pts, voxel_size=1e-30)


def test_pipeline_step_matches_the_function():
    pts = _scan(500, 5)
    mesh = mio.Mesh(pts, [("vertex", np.arange(len(pts)).reshape(-1, 1))])
    direct = mio.filter_points(mesh, statistical_k=10, std_ratio=1.5, voxel_size=0.05)
    piped = mio.apply_operations(
        mesh,
        [
            {
                "Op": "FilterPoints",
                "StatisticalK": 10,
                "StdRatio": 1.5,
                "VoxelSize": 0.05,
            }
        ],
    )
    out = piped[0] if isinstance(piped, tuple) else piped
    assert np.array_equal(np.asarray(out.points), np.asarray(direct.points))


def test_cli_json(tmp_path, capsys):
    from meshioplusplus._cli._main import main

    src = str(tmp_path / "scan.vtu")
    pts = _scan(800, 6)
    mio.write(src, mio.Mesh(pts, [("vertex", np.arange(len(pts)).reshape(-1, 1))]))
    out = str(tmp_path / "clean.vtu")
    assert main(
        ["filter-points", src, out, "--outliers", "statistical", "--k", "10", "--std-ratio", "2",
         "--voxel", "0.05", "--json"]
    ) in (0, None)  # fmt: skip
    rep = json.loads(capsys.readouterr().out)
    assert rep["num_input"] == len(pts)
    assert rep["num_kept"] == len(mio.read(out).points)
    assert rep["num_removed_statistical"] >= 6
