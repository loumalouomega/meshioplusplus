"""Boundary-driven morphing by harmonic extension (roadmap 7.4.1): lifting a
hexahedral box's top face inverts no cell (the roadmap probe), a rigid motion
of the boundary moves everything rigidly, stiffening and steps keep a large
rotation valid where plain Laplace inverts cells, the constrained-point
choices, the refusals, the pipeline step in both engines."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "morph"), reason="needs the C++ core"
)


def _box(n=6):
    return mp.box(dims=(n, n, n))


def _annulus():
    return mp.convert_cells(
        mp.annulus(inner_radius=0.3, outer_radius=1.5, ntheta=48, nr=12),
        mode="simplexify",
    )


def _rotate_inner(mesh, degrees):
    r = np.linalg.norm(mesh.points[:, :2], axis=1)
    inner = np.isclose(r, 0.3)
    th = np.radians(degrees)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    d = np.zeros((len(mesh.points), 3))
    d[inner, :2] = mesh.points[inner, :2] @ R.T - mesh.points[inner, :2]
    return d


def test_lifting_a_hex_box_top_inverts_nothing():
    """The probe: lift the top face by 0.1 (a bump), the interior follows."""
    box = _box()
    x, y, z = box.points.T
    d = np.zeros((len(box.points), 3))
    top = z == 1.0
    d[top, 2] = 0.1 * 16 * x[top] * (1 - x[top]) * y[top] * (1 - y[top])
    out, report = mp.morph(box, d, return_report=True)
    assert report["num_inverted"] == 0 and report["converged"]
    assert report["num_free"] == 5**3
    assert report["max_displacement"] == pytest.approx(0.1)
    np.testing.assert_array_equal(out.cells[0].data, box.cells[0].data)
    lift = out.points[:, 2] - z
    assert lift.min() >= -1e-12  # nothing moves down
    np.testing.assert_allclose(out.points[top, 2], 1.0 + d[top, 2], atol=1e-15)


def test_rigid_boundary_motion_moves_everything_rigidly():
    box = _box(4)
    shift = np.tile([0.3, -0.2, 0.1], (len(box.points), 1))
    out = mp.morph(box, shift)
    np.testing.assert_allclose(out.points - box.points, shift, atol=1e-9)


def test_stiffening_and_steps_keep_a_large_rotation_valid():
    mesh = _annulus()
    d = _rotate_inner(mesh, 90)
    _, plain = mp.morph(
        mesh, d, stiffening=0.0, steps=4, allow_inverted=True, return_report=True
    )
    assert plain["num_inverted"] > 0
    with pytest.raises(ValueError, match="inside out"):
        mp.morph(mesh, d, stiffening=0.0, steps=4)
    out, stiff = mp.morph(mesh, d, stiffening=1.0, steps=4, return_report=True)
    assert stiff["num_inverted"] == 0


def test_constrained_points_and_array_names():
    box = _box(4)
    z = box.points[:, 2]
    d = np.zeros((len(box.points), 3))
    d[:, 0] = 0.05 * z  # a shear, prescribed on the bottom and top only
    pinned = np.flatnonzero((z == 0.0) | (z == 1.0))
    box.point_sets["ends"] = pinned
    box.point_data["d"] = d
    out, rep = mp.morph(
        box, "d", fix_boundary=False, constrained="ends", return_report=True
    )
    assert rep["num_constrained"] == len(pinned)
    # the harmonic extension of a linear field is that field
    np.testing.assert_allclose(out.points[:, 0] - box.points[:, 0], d[:, 0], atol=1e-9)
    mask = np.zeros(len(box.points), dtype=bool)
    mask[pinned] = True
    again = mp.morph(box, "d", fix_boundary=False, constrained=mask)
    np.testing.assert_array_equal(again.points, out.points)
    assert "morph:displacement" not in mp.morph(box, d).point_data


def test_refusals():
    box = _box(2)
    with pytest.raises(ValueError, match="no point data"):
        mp.morph(box, "missing")
    with pytest.raises(ValueError, match="one row per point"):
        mp.morph(box, np.zeros((3, 3)))
    with pytest.raises(ValueError, match="steps"):
        mp.morph(box, np.zeros((len(box.points), 3)), steps=0)
    with pytest.raises(ValueError, match="stiffening"):
        mp.morph(box, np.zeros((len(box.points), 3)), stiffening=-1.0)
    line = mp.Mesh(np.eye(3), [("line", np.array([[0, 1], [1, 2]]))])
    with pytest.raises(ValueError, match="no surface or volume"):
        mp.morph(line, np.zeros((3, 3)))


def test_pipeline_step_in_both_engines(tmp_path):
    mesh = _annulus()
    mesh.point_data["d"] = _rotate_inner(mesh, 45)
    src = tmp_path / "a.vtu"
    mp.write(str(src), mesh, compression=None)
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "Morph", "Array": "d", "Steps": 2}],
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
    assert rep_py["steps"][0] == rep_cpp["steps"][0]
