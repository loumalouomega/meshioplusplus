"""Untangling inverted tetrahedra and hexahedra (roadmap 7.4.2): a
tetrahedral cube with a node pushed through a face ends with every signed
volume positive (the roadmap probe), a heavily jittered mesh, hexahedra, the
boundary kept, the connectivity kept, the pipeline step in both engines."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "untangle"), reason="needs the C++ core"
)


def _inverted(mesh):
    return mp.check_quality(mesh)["num_inverted"]


def _tet_cube(n=4):
    return mp.convert_cells(mp.box(dims=(n, n, n)), mode="simplexify")


def _jittered(mesh, fraction, scale, seed=0):
    rng = np.random.default_rng(seed)
    p = mesh.points.copy()
    interior = np.flatnonzero(np.all((p > 0) & (p < 1), axis=1))
    pick = rng.choice(interior, max(1, int(len(interior) * fraction)), replace=False)
    p[pick] += rng.normal(scale=scale, size=(len(pick), 3))
    return mp.Mesh(p, mesh.cells)


def test_node_pushed_through_a_face():
    """The probe."""
    tets = _tet_cube()
    p = tets.points.copy()
    p[(1 * 5 + 1) * 5 + 1] = [0.62, 0.58, 0.6]
    tangled = mp.Mesh(p, tets.cells)
    assert _inverted(tangled) > 0
    out, report = mp.untangle(tangled, return_report=True)
    assert _inverted(out) == 0
    assert report["num_inverted_before"] > 0 and report["num_inverted_after"] == 0
    assert report["min_quality_after"] > 0
    np.testing.assert_array_equal(out.cells[0].data, tangled.cells[0].data)
    boundary = np.any((p == 0) | (p == 1), axis=1)
    np.testing.assert_array_equal(out.points[boundary], p[boundary])


def test_heavy_jitter_and_only_near_tangles():
    mesh = _jittered(_tet_cube(10), 0.2, 0.08)
    assert _inverted(mesh) > 100
    out, report = mp.untangle(mesh, return_report=True)
    assert report["num_inverted_after"] == 0 and _inverted(out) == 0
    # a single tangle moves only its neighbourhood
    one = _tet_cube(10)
    p = one.points.copy()
    centre = int(np.argmin(np.linalg.norm(p - 0.5, axis=1)))
    p[centre] += [0.0, 0.0, 0.13]
    fixed, rep = mp.untangle(mp.Mesh(p, one.cells), return_report=True)
    moved = np.flatnonzero(np.abs(fixed.points - p).max(axis=1) > 0)
    assert rep["num_inverted_after"] == 0
    assert 0 < len(moved) < 0.1 * len(p)


def test_hexahedra():
    box = mp.box(dims=(6, 6, 6))
    p = box.points.copy()
    centre = int(np.argmin(np.linalg.norm(p - 0.5, axis=1)))
    p[centre] += [0.0, 0.0, 0.3]
    out, report = mp.untangle(mp.Mesh(p, box.cells), return_report=True)
    assert report["num_inverted_before"] > 0
    assert report["num_inverted_after"] == 0 and report["min_quality_after"] > 0


def test_smooth_improves_quality_of_a_valid_mesh():
    mesh = _jittered(_tet_cube(6), 0.5, 0.02, seed=3)
    assert _inverted(mesh) == 0
    same, rep0 = mp.untangle(mesh, return_report=True)
    assert rep0["sweeps"] == 0
    np.testing.assert_array_equal(same.points, mesh.points)
    better, rep = mp.untangle(mesh, smooth=True, return_report=True)
    assert rep["min_quality_after"] > rep["min_quality_before"]


def test_frozen_and_refusals():
    tets = _tet_cube()
    p = tets.points.copy()
    node = (1 * 5 + 1) * 5 + 1
    p[node] = [0.62, 0.58, 0.6]
    tangled = mp.Mesh(p, tets.cells)
    tangled.point_sets["stuck"] = np.arange(len(p))
    out, rep = mp.untangle(tangled, frozen="stuck", return_report=True)
    assert rep["num_inverted_after"] == rep["num_inverted_before"] > 0
    with pytest.raises(ValueError, match="no tetrahedron or hexahedron"):
        mp.untangle(mp.icosa_sphere(subdivisions=1))
    with pytest.raises(ValueError, match="no point_set"):
        mp.untangle(tangled, frozen="nothing")


def test_pipeline_step_in_both_engines(tmp_path):
    src = tmp_path / "t.vtu"
    mp.write(str(src), _jittered(_tet_cube(6), 0.3, 0.07, seed=1), compression=None)
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "Untangle"}],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep_py = mp.run_pipeline(doc)
    py = mp.read(str(tmp_path / "py.vtu"))
    assert _inverted(py) == 0
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    rep_cpp = _core.run_pipeline_json(json.dumps(doc))
    cpp = mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(cpp.points, py.points)
    assert rep_py["steps"][0] == rep_cpp["steps"][0]
