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


# --- alignment beyond local ICP (roadmap 7.3.1 - 7.3.3, 7.3.5, 7.3.6) ---------


def _box_surface(lo, hi, k):
    """An axis-aligned box's surface: k x k quads per face, two outward
    triangles each."""
    pts, tris = [], []
    for a in range(3):
        b, c = (a + 1) % 3, (a + 2) % 3
        for side in (0, 1):
            base = len(pts)
            for i in range(k + 1):
                for j in range(k + 1):
                    p = [0.0, 0.0, 0.0]
                    p[a] = hi[a] if side else lo[a]
                    p[b] = lo[b] + (hi[b] - lo[b]) * i / k
                    p[c] = lo[c] + (hi[c] - lo[c]) * j / k
                    pts.append(p)
            for i in range(k):
                for j in range(k):
                    v = [base + i * (k + 1) + j, base + (i + 1) * (k + 1) + j]
                    v += [base + (i + 1) * (k + 1) + j + 1, base + i * (k + 1) + j + 1]
                    t1, t2 = [v[0], v[1], v[2]], [v[0], v[2], v[3]]
                    if not side:
                        t1, t2 = t1[::-1], t2[::-1]
                    tris += [t1, t2]
    return np.array(pts, float), np.array(tris, np.int64)


def _surface(boxes):
    pts, tris, base = [], [], 0
    for lo, hi, k in boxes:
        p, t = _box_surface(lo, hi, k)
        pts.append(p)
        tris.append(t + base)
        base += len(p)
    return mp.Mesh(np.vstack(pts), [("triangle", np.vstack(tris))])


def _bracket():
    """An L bracket with a rib in one corner: no proper symmetry."""
    return _surface(
        [
            ((0, 0, 0), (4, 1, 0.5), 8),
            ((0, 0, 0.5), (0.5, 1, 3), 8),
            ((3, 0.6, 0.5), (4, 1, 1.2), 4),
        ]
    )


def test_initial_auto_registers_a_large_rotation_back():
    target = _bracket()
    motion = _motion(np.radians(150.0), (2.0, -1.0, 0.5))
    source = mp.Mesh(_apply(motion, target.points), target.cells)
    m, report = mp.register_rigid(
        source,
        target,
        method="point_to_plane",
        initial="auto",
        tolerance=1e-14,
        max_iterations=200,
    )
    np.testing.assert_allclose(m @ motion, np.eye(4), atol=1e-9)
    assert len(report["candidate_rms"]) == 24
    assert report["chosen_candidate"] >= 0
    assert report["warnings"] == []
    with pytest.raises(ValueError, match="initial"):
        mp.register_rigid(source, target, initial="nope")


def test_max_distance_auto_and_regions():
    target = _bracket()
    m, report = mp.register_rigid(
        mp.Mesh(_apply(_motion(0.05), target.points), target.cells),
        target,
        method="point_to_plane",
        max_distance="auto",
    )
    assert report["num_rejected"] >= 0
    np.testing.assert_allclose(m @ _motion(0.05), np.eye(4), atol=1e-6)
    with pytest.raises(ValueError, match="nope"):
        mp.register_rigid(target, target, source_region="nope")


def test_register_landmarks():
    motion = _motion(2.0, (1.0, -2.0, 3.0))
    src = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0.0]])
    m, report = mp.register_landmarks(src, _apply(motion, src))
    np.testing.assert_allclose(m, motion, atol=1e-12)
    assert report["rms"] < 1e-12 and len(report["residuals"]) == 3
    with pytest.raises(ValueError, match="three"):
        mp.register_landmarks(src[:2], src[:2])
    with pytest.raises(ValueError, match="collinear"):
        mp.register_landmarks([[0, 0, 0], [1, 0, 0], [2, 0, 0]], [[0, 0, 0]] * 3)


def test_find_symmetry_plane_and_mirror_auto():
    surf = _surface(
        [
            ((-2, 0, 0), (2, 1, 0.5), 8),
            ((-0.25, 0, 0.5), (0.25, 1, 3), 8),
            ((1.5, 0.6, 0.5), (2, 1, 1.2), 4),
            ((-2, 0.6, 0.5), (-1.5, 1, 1.2), 4),
        ]
    )
    motion = _motion(0.6, (0.3, -0.2, 0.1))
    mesh = mp.Mesh(_apply(motion, surf.points), surf.cells)
    plane = mp.find_symmetry_plane(mesh)
    assert plane["share"] > 0.9
    normals = motion[:3, :3].T  # rows: the images of the axes
    assert max(abs(np.dot(plane["normal"], n)) for n in normals) > 1 - 1e-8
    out, report = mp.mirror(mesh, plane="auto", return_report=True)
    assert report["plane"]["normal"] == plane["normal"]
    assert len(out.points) == len(mesh.points)
    with pytest.raises(ValueError, match="plane"):
        mp.mirror(mesh, plane="xy")


def test_auto_align_lays_a_block_on_its_largest_face():
    block = _surface([((0, 0, 0), (3, 2, 1), 6)])
    moved = mp.Mesh(_apply(_motion(0.9, (4, 5, 6)), block.points), block.cells)
    out, report = mp.auto_align(moved, return_report=True)
    assert report["mode"] == "face" and report["share"] > 0.9
    ext = out.points.max(axis=0) - out.points.min(axis=0)
    np.testing.assert_allclose(ext, [2, 3, 1], atol=1e-9)
    assert abs(out.points[:, 2].min()) < 1e-12
    np.testing.assert_allclose(report["matrix"].shape, (4, 4))


def test_deviation_register():
    reference = _bracket()
    scan = mp.Mesh(_apply(_motion(2.5, (1, 2, 3)), reference.points), reference.cells)
    _, before = mp.deviation(scan, reference, tolerance=1e-6, return_report=True)
    assert before["share_within"] < 0.5
    out, report = mp.deviation(
        scan, reference, tolerance=1e-6, register=True, return_report=True
    )
    assert report["registered"] and report["share_within"] > 0.999
    np.testing.assert_allclose(out.points, reference.points, atol=1e-8)
    with pytest.raises(ValueError, match="primitive"):
        mp.deviation(scan, ("plane", [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]), register=True)


def test_alignment_steps_in_both_engines(tmp_path):
    reference = _bracket()
    ref = tmp_path / "ref.vtu"
    scan = tmp_path / "scan.vtu"
    mp.write(str(ref), reference, compression=None)
    motion = _motion(np.radians(150.0), (2.0, -1.0, 0.5))
    mp.write(
        str(scan),
        mp.Mesh(_apply(motion, reference.points), reference.cells),
        compression=None,
    )
    src = [0, 0, 0, 1, 0, 0, 0, 2, 0]
    dst = list(_apply(motion, np.reshape(src, (3, 3))).ravel())
    doc = {
        "Version": 1,
        "Input": {"Path": str(scan)},
        "Operations": [
            {
                "Op": "RegisterRigid",
                "Target": str(ref),
                "Method": "point_to_plane",
                "Initial": "auto",
                "AdaptiveTrim": 3.0,
                "MaxIterations": 200,
            },
            {
                "Op": "Deviation",
                "Reference": str(ref),
                "Register": True,
                "Tolerance": 1e-6,
            },
            {"Op": "RegisterLandmarks", "SourcePoints": dst, "TargetPoints": src},
            {"Op": "RegisterLandmarks", "SourcePoints": src, "TargetPoints": dst},
            {"Op": "AutoAlign"},
            {"Op": "Mirror", "Plane": "auto"},
        ],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    with pytest.raises(ValueError, match="no symmetry plane"):
        mp.run_pipeline(doc)
    doc["Operations"].pop()
    rep_py = mp.run_pipeline(doc)
    py = mp.read(str(tmp_path / "py.vtu"))
    assert abs(py.points[:, 2].min()) < 1e-9
    assert rep_py["steps"][1]["ShareWithin"] > 0.999
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    rep_cpp = _core.run_pipeline_json(json.dumps(doc))
    cpp = mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(cpp.points, py.points)
    assert rep_py["steps"] == rep_cpp["steps"]
