"""Geodesic distance by the heat method (roadmap 7.3.1): accuracy on spheres
and planes, sources by id, region and point set, the scope refusals, the
pipeline step in both engines."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core
from meshioplusplus._regions import Region

from . import helpers

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "geodesic_distance"), reason="needs the C++ core"
)


def _sphere(subdivisions):
    m = mp.icosa_sphere(subdivisions=subdivisions)
    p = m.points / np.linalg.norm(m.points, axis=1)[:, None]
    return mp.Mesh(p, m.cells)


def _distance(mesh):
    return np.asarray(mesh.point_data["geodesic:distance"])


def test_sphere_pole_to_equator():
    """The probe: pole to equator is pi/2 within about 1e-2, and the error
    shrinks as the sphere is refined."""
    errors = []
    for sub in (3, 4, 5):
        mesh = _sphere(sub)
        pole = int(np.argmax(mesh.points[:, 2]))
        out, report = mp.geodesic_distance(mesh, pole, return_report=True)
        d = _distance(out)
        equator = np.abs(mesh.points[:, 2]) < 1e-9
        errors.append(np.abs(d[equator] - np.pi / 2).max())
        exact = np.arccos(np.clip(mesh.points @ mesh.points[pole], -1, 1))
        assert np.abs(d - exact).mean() < 0.025
        assert d[pole] == 0.0
        assert report["converged"] and report["num_sources"] == 1
        assert report["num_unreached"] == 0
    assert errors[-1] < 1.3e-2
    assert errors[0] > errors[1] > errors[2]


def test_larger_time_factor_is_closer_on_a_sphere():
    mesh = _sphere(4)
    pole = int(np.argmax(mesh.points[:, 2]))
    exact = np.arccos(np.clip(mesh.points @ mesh.points[pole], -1, 1))
    e1 = np.abs(_distance(mp.geodesic_distance(mesh, pole)) - exact).mean()
    e4 = np.abs(
        _distance(mp.geodesic_distance(mesh, pole, time_factor=4.0)) - exact
    ).mean()
    assert e4 < e1


def test_plane_distance_is_euclidean():
    m, _ = mp.apply_operations(
        None,
        [
            {
                "Op": "Generate",
                "Shape": "rectangle",
                "Dims": [40, 40],
                "CellType": "triangle",
            }
        ],
    )
    centre = int(np.argmin(np.linalg.norm(m.points[:, :2] - 0.5, axis=1)))
    d = _distance(mp.geodesic_distance(m, centre))
    exact = np.linalg.norm(m.points[:, :2] - m.points[centre, :2], axis=1)
    assert np.abs(d - exact).mean() < 0.02


def test_two_sources_take_the_nearer():
    mesh = _sphere(4)
    north = int(np.argmax(mesh.points[:, 2]))
    south = int(np.argmin(mesh.points[:, 2]))
    d = _distance(mp.geodesic_distance(mesh, [north, south]))
    assert d[north] == pytest.approx(0.0, abs=2e-2)
    assert d[south] == pytest.approx(0.0, abs=2e-2)
    assert d.max() == pytest.approx(np.pi / 2, abs=3e-2)


def test_sources_by_region_and_point_set():
    mesh = _sphere(3)
    pole = int(np.argmax(mesh.points[:, 2]))
    by_id = _distance(mp.geodesic_distance(mesh, [pole]))
    tagged = mp.Mesh(mesh.points, mesh.cells)
    tagged.regions = [Region("top", "point", [pole])]
    np.testing.assert_allclose(_distance(mp.geodesic_distance(tagged, "top")), by_id)
    np.testing.assert_allclose(
        _distance(mp.geodesic_distance(tagged, source_region="top")), by_id
    )


def test_unreached_pieces_and_isolated_points_are_nan():
    a = _sphere(2)
    b = mp.Mesh(a.points + [5.0, 0.0, 0.0], a.cells)
    n = len(a.points)
    pts = np.vstack([a.points, b.points, [[9.0, 9.0, 9.0]]])
    cells = [("triangle", np.vstack([a.cells[0].data, b.cells[0].data + n]))]
    out, report = mp.geodesic_distance(mp.Mesh(pts, cells), 0, return_report=True)
    d = _distance(out)
    assert np.all(np.isfinite(d[:n]))
    assert np.all(np.isnan(d[n:]))
    assert report["num_unreached"] == n + 1


def test_input_is_not_modified_and_data_rides_along():
    mesh = _sphere(2)
    mesh.point_data["t"] = np.arange(len(mesh.points), dtype=float)
    out = mp.geodesic_distance(mesh, 0)
    assert "geodesic:distance" not in mesh.point_data
    np.testing.assert_array_equal(out.point_data["t"], mesh.point_data["t"])
    np.testing.assert_array_equal(out.points, mesh.points)


def test_refusals():
    mesh = _sphere(2)
    with pytest.raises(ValueError, match="at least one source"):
        mp.geodesic_distance(mesh, [])
    with pytest.raises(ValueError, match="outside"):
        mp.geodesic_distance(mesh, [10**6])
    with pytest.raises(ValueError, match="time factor"):
        mp.geodesic_distance(mesh, 0, time_factor=0.0)
    with pytest.raises(ValueError, match="extract_surface"):
        mp.geodesic_distance(helpers.tet_mesh, 0)
    with pytest.raises(ValueError, match="simplexify"):
        mp.geodesic_distance(helpers.quad_mesh, 0)
    lonely = mp.Mesh(np.vstack([mesh.points, [[3, 3, 3]]]), mesh.cells)
    with pytest.raises(ValueError, match="in no triangle"):
        mp.geodesic_distance(lonely, len(mesh.points))


def test_pipeline_step_in_both_engines(tmp_path):
    mesh = _sphere(3)
    src = tmp_path / "s.vtu"
    mp.write(str(src), mesh, compression=None)
    pole = int(np.argmax(mesh.points[:, 2]))
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "GeodesicDistance", "Sources": [pole]}],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep_py = mp.run_pipeline(doc)
    d_py = _distance(mp.read(str(tmp_path / "py.vtu")))
    expected = _distance(mp.geodesic_distance(mesh, pole))
    np.testing.assert_allclose(d_py, expected, rtol=0, atol=1e-12)
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    rep_cpp = _core.run_pipeline_json(json.dumps(doc))
    np.testing.assert_array_equal(_distance(mp.read(str(tmp_path / "cpp.vtu"))), d_py)
    assert rep_py["steps"][0] == rep_cpp["steps"][0]
