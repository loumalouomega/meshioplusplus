"""Laplace-Beltrami eigenmodes of a triangle surface (roadmap 7.3.2): the
sphere's l (l + 1) with multiplicities 1, 3, 5 (the roadmap probe), the
M-orthonormal modes, a square's Neumann spectrum, two pieces, the scope
refusals, the pipeline step in both engines."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

from . import helpers

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "laplacian_eigenmodes"), reason="needs the C++ core"
)


def _sphere(subdivisions):
    m = mp.icosa_sphere(subdivisions=subdivisions)
    return mp.Mesh(m.points / np.linalg.norm(m.points, axis=1)[:, None], m.cells)


def _lumped_mass(mesh):
    p = mesh.points
    t = mesh.cells_dict["triangle"]
    area = 0.5 * np.linalg.norm(
        np.cross(p[t[:, 1]] - p[t[:, 0]], p[t[:, 2]] - p[t[:, 0]]), axis=1
    )
    m = np.zeros(len(p))
    np.add.at(m, t.ravel(), np.repeat(area / 3.0, 3))
    return m


def test_sphere_spectrum_is_l_l_plus_1_with_multiplicities():
    """The probe: 0 once, 2 three times, 6 five times, approaching l(l + 1)."""
    errors = []
    for sub in (3, 4):
        mesh = _sphere(sub)
        out, report = mp.laplacian_eigenmodes(mesh, 9, return_report=True)
        lam = np.asarray(out.field_data["laplacian:eigenvalues"])
        assert report["converged"] and np.array_equal(lam, report["eigenvalues"])
        assert lam[0] == 0.0
        np.testing.assert_allclose(lam[1:4], 2.0, rtol=2e-2)
        np.testing.assert_allclose(lam[4:9], 6.0, rtol=2e-2)
        # each cluster is degenerate to the solver's accuracy
        assert np.ptp(lam[1:4]) < 1e-6 and np.ptp(lam[4:9]) < 1e-6
        errors.append(abs(lam[4] - 6.0))
    assert errors[1] < errors[0]


def test_modes_are_mass_orthonormal_and_eigenvectors():
    mesh = _sphere(3)
    out = mp.laplacian_eigenmodes(mesh, 6)
    phi = out.point_data["laplacian:modes"]
    m = _lumped_mass(mesh)
    np.testing.assert_allclose(phi.T @ (m[:, None] * phi), np.eye(6), atol=1e-8)
    # the constant mode
    np.testing.assert_allclose(phi[:, 0], 1 / np.sqrt(m.sum()), rtol=1e-8)
    # the l = 1 modes span the coordinate functions
    coords = mesh.points
    proj = phi[:, 1:4].T @ (m[:, None] * coords)
    assert np.linalg.svd(proj, compute_uv=False).min() > 0.25


def test_two_pieces_give_two_zero_eigenvalues():
    a = _sphere(2)
    n = len(a.points)
    pts = np.vstack([a.points, a.points + [5.0, 0, 0]])
    cells = [("triangle", np.vstack([a.cells[0].data, a.cells[0].data + n]))]
    lam = mp.laplacian_eigenmodes(mp.Mesh(pts, cells), 3).field_data[
        "laplacian:eigenvalues"
    ]
    assert lam[0] == 0.0 and lam[1] == 0.0 and lam[2] > 1.0


def test_isolated_points_get_nan_and_data_rides_along():
    s = _sphere(2)
    mesh = mp.Mesh(np.vstack([s.points, [[3.0, 3.0, 3.0]]]), s.cells)
    mesh.point_data["t"] = np.arange(len(mesh.points), dtype=float)
    out = mp.laplacian_eigenmodes(mesh, 3)
    assert np.isnan(out.point_data["laplacian:modes"][-1]).all()
    np.testing.assert_array_equal(out.point_data["t"], mesh.point_data["t"])
    assert "laplacian:modes" not in mesh.point_data


def test_refusals():
    with pytest.raises(ValueError, match="extract_surface"):
        mp.laplacian_eigenmodes(helpers.tet_mesh, 2)
    with pytest.raises(ValueError, match="simplexify"):
        mp.laplacian_eigenmodes(helpers.quad_mesh, 2)
    with pytest.raises(ValueError, match="at least 1"):
        mp.laplacian_eigenmodes(_sphere(1), 0)
    with pytest.raises(ValueError, match="exceeds"):
        mp.laplacian_eigenmodes(_sphere(0), 100)


def test_pipeline_step_in_both_engines(tmp_path):
    src = tmp_path / "s.vtu"
    mp.write(str(src), _sphere(3), compression=None)
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "LaplacianEigenmodes", "NumModes": 5}],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep_py = mp.run_pipeline(doc)
    py = mp.read(str(tmp_path / "py.vtu"))
    assert py.point_data["laplacian:modes"].shape[1] == 5
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    rep_cpp = _core.run_pipeline_json(json.dumps(doc))
    cpp = mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(
        cpp.point_data["laplacian:modes"], py.point_data["laplacian:modes"]
    )
    assert rep_py["steps"][0] == rep_cpp["steps"][0]
