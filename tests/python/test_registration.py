"""Rigid registration by iterative closest points (roadmap 7.5.3): a rotated
and shifted copy registers back within 1e-9 (the roadmap probe), both
methods, onto a surface and onto a cloud, outlier rejection, the refusals, the
pipeline step in both engines."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "register_rigid"), reason="needs the C++ core"
)


def _rotation(axis, angle):
    a = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K


def _motion(angle=0.12, shift=(0.15, -0.08, 0.05)):
    m = np.eye(4)
    m[:3, :3] = _rotation((1, 2, 3), angle)
    m[:3, 3] = shift
    return m


def _apply(m, pts):
    return pts @ m[:3, :3].T + m[:3, 3]


def _ellipsoid(subdivisions=3):
    s = mp.icosa_sphere(subdivisions=subdivisions)
    p = s.points / np.linalg.norm(s.points, axis=1)[:, None]
    p = p * [1.6, 1.0, 0.7] + 0.1 * np.sin(3 * p[:, [1, 2, 0]])
    return mp.Mesh(p, s.cells)


@pytest.mark.parametrize("method", ["point_to_point", "point_to_plane"])
def test_rotated_shifted_copy_registers_back(method):
    target = _ellipsoid()
    source = mp.Mesh(_apply(_motion(), target.points), [])
    cloud = mp.estimate_normals(target, k=10)  # a cloud target, with normals
    cloud = mp.Mesh(cloud.points, [], point_data=dict(cloud.point_data))
    # Point-to-point onto a surface converges only linearly (the closest
    # points slide along the faces); onto the cloud of the same points, and
    # point-to-plane onto either, it is exact within a few iterations.
    targets = (cloud,) if method == "point_to_point" else (target, cloud)
    for tgt in targets:
        m, report = mp.register_rigid(source, tgt, method=method)
        assert report["converged"], report
        np.testing.assert_allclose(m @ _motion(), np.eye(4), atol=1e-9)
        back = mp.transform(source, matrix=m)
        np.testing.assert_allclose(back.points, target.points, atol=1e-9)
        assert report["rms"] < 1e-9 < report["initial_rms"]
        assert report["num_pairs"] == len(target.points)


def test_point_to_plane_needs_fewer_iterations_on_a_surface():
    target = _ellipsoid(4)
    rng = np.random.default_rng(0)
    # samples off the vertices: face points of the target, moved
    tris = target.cells_dict["triangle"][rng.integers(0, 2000, 800)]
    w = rng.dirichlet([1, 1, 1], size=800)
    samples = np.einsum("ij,ijk->ik", w, target.points[tris])
    source = _apply(_motion(0.05, (0.03, 0.02, -0.02)), samples)
    _, p2p = mp.register_rigid(source, target, max_iterations=500, tolerance=1e-8)
    _, p2l = mp.register_rigid(source, target, method="point_to_plane", tolerance=1e-8)
    assert p2l["converged"] and p2l["rms"] < 1e-9
    assert p2l["iterations"] < p2p["iterations"]


def test_max_distance_drops_outliers():
    target = _ellipsoid()
    pts = _apply(_motion(0.02, (0.01, 0, 0)), target.points)
    outliers = np.array([[10.0, 10.0, 10.0], [-9.0, 4.0, 2.0]])
    source = np.vstack([pts, outliers])
    m, report = mp.register_rigid(
        source, target, method="point_to_plane", max_distance=0.5
    )
    assert report["converged"] and report["num_pairs"] == len(pts)
    np.testing.assert_allclose(m @ _motion(0.02, (0.01, 0, 0)), np.eye(4), atol=1e-8)
    m_all, _ = mp.register_rigid(source, target, method="point_to_plane")
    assert np.abs(m_all @ _motion(0.02, (0.01, 0, 0)) - np.eye(4)).max() > 1e-4


def test_initial_guess_is_honoured():
    target = _ellipsoid()
    motion = _motion(0.9, (1.0, 0.5, -0.3))  # too far for ICP alone
    source = _apply(motion, target.points)
    m, report = mp.register_rigid(source, target, initial=np.linalg.inv(motion))
    assert report["initial_rms"] < 1e-12
    np.testing.assert_allclose(m @ motion, np.eye(4), atol=1e-9)


def test_refusals(capfd):
    target = _ellipsoid(1)
    with pytest.raises(ValueError, match="unknown method"):
        mp.register_rigid(target, target, method="affine")
    with pytest.raises(ValueError, match="estimate_normals"):
        mp.register_rigid(target, target.points, method="point_to_plane")
    with pytest.raises(ValueError, match="extract_surface"):
        from . import helpers

        mp.register_rigid(target, helpers.tet_mesh)
    with pytest.raises(ValueError, match="max_distance"):
        mp.register_rigid(target.points + 100.0, target, max_distance=1.0)
    _, report = mp.register_rigid(
        _apply(_motion(), target.points), target, max_iterations=1
    )
    assert not report["converged"]
    assert "no convergence" in capfd.readouterr().err


def test_pipeline_step_in_both_engines(tmp_path):
    target = _ellipsoid()
    tgt = tmp_path / "target.vtu"
    src = tmp_path / "source.vtu"
    # Uncompressed: the build may have no zlib.
    mp.write(str(tgt), target, compression=None)
    mp.write(
        str(src),
        mp.Mesh(_apply(_motion(), target.points), target.cells),
        compression=None,
    )
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {"Op": "RegisterRigid", "Target": str(tgt), "Method": "point_to_plane"}
        ],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep_py = mp.run_pipeline(doc)
    py = mp.read(str(tmp_path / "py.vtu"))
    np.testing.assert_allclose(py.points, target.points, atol=1e-9)
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    rep_cpp = _core.run_pipeline_json(json.dumps(doc))
    cpp = mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(cpp.points, py.points)
    assert rep_py["steps"][0] == rep_cpp["steps"][0]
