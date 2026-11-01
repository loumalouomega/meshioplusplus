"""Cylindrical and spherical components about an axis (roadmap 7.6.1): the
recorded coordinates, vectors and tensors in local components, the inverse,
singular points, the core against the numpy twin, and the pipeline step in
both engines."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _coordinates, _core

from . import helpers


def _ring(n=12, radius=2.0):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pts = np.column_stack([radius * np.cos(t), radius * np.sin(t), 0.3 * np.arange(n)])
    lines = np.array([[i, (i + 1) % n] for i in range(n)])
    mesh = mp.Mesh(pts, [("line", lines)])
    # A rigid rotation about z: purely tangential, |v| = r.
    mesh.point_data["v"] = np.column_stack([-pts[:, 1], pts[:, 0], np.zeros(n)])
    return mesh


def test_rigid_rotation_is_tangential():
    out = mp.convert_coordinates(_ring(), "cylindrical", arrays=["v"])
    np.testing.assert_allclose(out.point_data["v"][:, 0], 0.0, atol=1e-12)
    np.testing.assert_allclose(out.point_data["v"][:, 1], 2.0, atol=1e-12)
    c = out.point_data["coords:cylindrical"]
    np.testing.assert_allclose(c[:, 0], 2.0, atol=1e-12)
    np.testing.assert_allclose(c[1, 1], 2 * np.pi / 12, atol=1e-12)
    np.testing.assert_allclose(out.points, _ring().points)  # points stay put


def test_sphere_normals_are_radial_in_spherical_components():
    s = mp.icosa_sphere(subdivisions=2)
    s.point_data["n"] = s.points / np.linalg.norm(s.points, axis=1)[:, None]
    out, report = mp.convert_coordinates(
        s, "spherical", arrays=["n"], return_report=True
    )
    np.testing.assert_allclose(
        out.point_data["n"], [[1.0, 0.0, 0.0]] * len(s.points), atol=1e-12
    )
    coords = out.point_data["coords:spherical"]
    np.testing.assert_allclose(coords[:, 0], 1.0, atol=1e-12)
    assert np.all((coords[:, 1] >= 0) & (coords[:, 1] <= np.pi))
    # The two poles lie on the axis: their azimuthal direction is undefined
    # (counted singular), yet e_r there is still the outward normal.
    on_axis = np.hypot(s.points[:, 0], s.points[:, 1]) < 1e-12
    assert on_axis.sum() == 2
    assert report["num_singular"] == 2


def test_tensors_and_cell_data_round_trip():
    rng = np.random.default_rng(1)
    mesh = mp.Mesh(helpers.tet_mesh.points + [0.3, -0.2, 0.1], helpers.tet_mesh.cells)
    mesh.point_data["S"] = rng.normal(size=(len(mesh.points), 3, 3))
    mesh.point_data["F"] = rng.normal(size=(len(mesh.points), 9))
    mesh.cell_data["w"] = [rng.normal(size=(len(b.data), 3)) for b in mesh.cells]
    for system in ("cylindrical", "spherical"):
        kw = dict(axis=(1.0, 2.0, 0.5), origin=(0.1, 0.2, 0.0), arrays=["S", "F", "w"])
        fwd = mp.convert_coordinates(mesh, system, **kw)
        back = mp.convert_coordinates(
            fwd, system, inverse=True, record_coordinates=False, **kw
        )
        for name in ("S", "F"):
            np.testing.assert_allclose(
                back.point_data[name], mesh.point_data[name], atol=1e-12
            )
        np.testing.assert_allclose(
            back.cell_data["w"][0], mesh.cell_data["w"][0], atol=1e-12
        )
        # The trace of a tensor is invariant under the rotation.
        np.testing.assert_allclose(
            np.trace(fwd.point_data["S"], axis1=1, axis2=2),
            np.trace(mesh.point_data["S"], axis1=1, axis2=2),
            atol=1e-12,
        )


@pytest.mark.skipif(
    not hasattr(_core, "convert_coordinates"), reason="needs the C++ core"
)
@pytest.mark.parametrize("system", ["cylindrical", "spherical"])
def test_core_matches_numpy_twin(system):
    rng = np.random.default_rng(2)
    pts = rng.normal(size=(40, 3))
    pts[0] = [0.0, 0.0, 3.0]  # on the axis
    pts[1] = [0.0, 0.0, 0.0]  # at the centre
    tris = np.array([[i, i + 1, i + 2] for i in range(0, 39, 3)])
    mesh = mp.Mesh(pts, [("triangle", tris)])
    mesh.point_data["u"] = rng.normal(size=(40, 3))
    mesh.cell_data["T"] = [rng.normal(size=(len(tris), 9))]
    kw = dict(
        axis=(0, 0, 2.0), origin=(0, 0, 0), reference=(1, 1, 0), arrays=["u", "T"]
    )
    core, rc = mp.convert_coordinates(mesh, system, return_report=True, **kw)
    twin, rt = _coordinates._convert_py(
        mesh,
        system,
        kw["axis"],
        kw["origin"],
        kw["reference"],
        kw["arrays"],
        False,
        True,
    )
    assert rc == rt
    assert rc["num_singular"] == 2  # the centre lies on the axis too
    for name in ("u", f"coords:{system}"):
        np.testing.assert_allclose(
            core.point_data[name], twin.point_data[name], atol=1e-13
        )
    np.testing.assert_allclose(
        core.cell_data["T"][0], twin.cell_data["T"][0], atol=1e-13
    )


def test_refusals():
    ring = _ring()
    with pytest.raises(ValueError, match="unknown system"):
        mp.convert_coordinates(ring, "polar")
    with pytest.raises(ValueError, match="nonzero"):
        mp.convert_coordinates(ring, axis=(0, 0, 0))
    with pytest.raises(ValueError, match="no point or cell data"):
        mp.convert_coordinates(ring, arrays=["missing"])
    ring.point_data["s"] = np.zeros(len(ring.points))
    with pytest.raises(ValueError, match="neither a vector"):
        mp.convert_coordinates(ring, arrays=["s"])


def test_pipeline_step_in_both_engines(tmp_path):
    src = tmp_path / "ring.vtu"
    mp.write(str(src), _ring(), compression=None)
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {"Op": "ConvertCoordinates", "System": "cylindrical", "Arrays": ["v"]}
        ],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep_py = mp.run_pipeline(doc)
    py = mp.read(str(tmp_path / "py.vtu"))
    np.testing.assert_allclose(py.point_data["v"][:, 1], 2.0, atol=1e-12)
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    rep_cpp = _core.run_pipeline_json(json.dumps(doc))
    cpp = mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_allclose(cpp.point_data["v"], py.point_data["v"], atol=1e-13)
    assert rep_py["steps"][0] == rep_cpp["steps"][0]
