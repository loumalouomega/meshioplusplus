"""The mesh's P1 operators as sparse matrices (roadmap 15.2.1)."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import mesh_operator

scipy_sparse = pytest.importorskip("scipy.sparse")


def _square():
    pts = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    return mp.Mesh(pts, [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))])


def test_square_reference_values():
    mesh = _square()
    lumped = mesh_operator(mesh, "mass_lumped")
    # corner 0 and 2 sit in both triangles (2 * 0.5/3), 1 and 3 in one
    np.testing.assert_allclose(lumped.diagonal(), [1 / 3, 1 / 6, 1 / 3, 1 / 6])
    consistent = mesh_operator(mesh, "mass_consistent")
    assert consistent.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(consistent.toarray(), consistent.toarray().T)
    k = mesh_operator(mesh, "laplacian")
    # the right-angle cotangents: K_01 = -1/2, K_02 = 0
    assert k[0, 1] == pytest.approx(-0.5) and k[0, 2] == pytest.approx(0.0)
    np.testing.assert_allclose(k @ np.ones(4), 0.0, atol=1e-14)


def test_gradient_is_exact_on_linear_fields_and_factors_the_stiffness():
    mesh = mp.icosa_sphere(1.0, 2)
    pts = mesh.points
    g = mesh_operator(mesh, "gradient")
    ntri = len(mesh.cells[0].data)
    assert g.shape == (3 * ntri, len(pts))
    f = 2 * pts[:, 0] - 3 * pts[:, 1] + 0.5 * pts[:, 2]
    grad = (g @ f).reshape(ntri, 3)
    # the gradient of a linear field restricted to a plane is the tangential part
    tris = np.asarray(mesh.cells[0].data)
    a, b, c = (pts[tris[:, i]] for i in range(3))
    normal = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(normal, axis=1)
    unit = normal / (2 * area)[:, None]
    full = np.array([2.0, -3.0, 0.5])
    expect = full - (unit @ full)[:, None] * unit
    np.testing.assert_allclose(grad, expect, atol=1e-12)
    # K = G^T diag(area) G
    weights = scipy_sparse.diags(np.repeat(area, 3))
    k = mesh_operator(mesh, "laplacian")
    diff = (g.T @ weights @ g - k).toarray()
    assert np.abs(diff).max() < 1e-12


def test_masses_sum_to_the_area_and_the_laplacian_is_semidefinite():
    mesh = mp.icosa_sphere(1.0, 2)
    lumped = mesh_operator(mesh, "mass_lumped")
    consistent = mesh_operator(mesh, "mass_consistent")
    assert lumped.sum() == pytest.approx(consistent.sum())
    k = mesh_operator(mesh, "laplacian").toarray()
    np.testing.assert_allclose(k, k.T, atol=1e-13)
    assert np.linalg.eigvalsh(k).min() > -1e-10


def test_the_tuple_form_without_scipy_and_the_report():
    (indptr, indices, values, shape), report = mesh_operator(
        _square(), "mass_lumped", sparse=False, return_report=True
    )
    assert shape == (4, 4) and report == {"num_degenerate": 0}
    assert len(indptr) == 5 and len(indices) == len(values) == 4


def test_a_degenerate_triangle_has_empty_gradient_rows():
    mesh = _square()
    mesh.points = np.vstack([mesh.points, [[2.0, 0, 0]]])
    mesh.cells[0] = mp.CellBlock(
        "triangle", np.array([[0, 1, 2], [0, 2, 3], [0, 1, 4]])  # collinear
    )
    g, report = mesh_operator(mesh, "gradient", return_report=True)
    assert report["num_degenerate"] == 1 and g.shape == (9, 5)
    assert g[6:].nnz == 0


@pytest.mark.parametrize(
    "build,needle",
    [
        (lambda: (mp.box(), "laplacian"), "extract_surface"),
        (lambda: (_square(), "bogus"), "unknown kind"),
    ],
)
def test_refusals(build, needle):
    mesh, kind = build()
    with pytest.raises(ValueError, match=needle):
        mesh_operator(mesh, kind)


def test_cli_verb(tmp_path, capsys):
    from meshioplusplus._cli._main import main

    src = str(tmp_path / "m.vtu")
    mp.write(src, _square())
    out = str(tmp_path / "k.npz")
    rc = main(["mesh-operator", src, out, "--kind", "laplacian", "--json"])
    assert rc in (0, None)
    info = json.loads(capsys.readouterr().out)
    assert info["shape"] == [4, 4]
    z = np.load(out)
    k = scipy_sparse.csr_matrix(
        (z["values"], z["indices"], z["indptr"]), shape=tuple(z["shape"])
    )
    np.testing.assert_allclose(k @ np.ones(4), 0.0, atol=1e-14)


def test_smooth_field_matches_the_screened_poisson_solve():
    from scipy.sparse.linalg import spsolve

    mesh = mp.icosa_sphere(1.0, 2)
    rng = np.random.default_rng(3)
    f = rng.normal(size=len(mesh.points))
    mesh.point_data["f"] = f
    mesh.point_data["v"] = np.column_stack([f, -2 * f])
    out, report = mp.smooth_field(mesh, "f", 0.4, return_report=True)
    k = mesh_operator(mesh, "laplacian")
    m = mesh_operator(mesh, "mass_lumped")
    expect = spsolve((m + 0.16 * k).tocsc(), m @ f)
    np.testing.assert_allclose(out.point_data["f:smooth"], expect, atol=1e-9)
    assert report["converged"] and "f" in out.point_data
    # the mass-weighted mean is kept, and the variance drops
    w = m.diagonal()
    assert w @ out.point_data["f:smooth"] == pytest.approx(w @ f, abs=1e-9)
    assert out.point_data["f:smooth"].std() < f.std()
    # vector arrays are smoothed per component; length 0 is the identity
    v = mp.smooth_field(mesh, "v", 0.4, output="vs").point_data["vs"]
    np.testing.assert_allclose(v[:, 0], expect, atol=1e-9)
    np.testing.assert_allclose(v[:, 1], -2 * expect, atol=1e-9)
    same = mp.smooth_field(mesh, "f", 0.0).point_data["f:smooth"]
    np.testing.assert_allclose(same, f, atol=1e-9)


def test_smooth_field_refusals_and_the_cli(tmp_path, capsys):
    from meshioplusplus._cli._main import main

    mesh = mp.icosa_sphere(1.0, 1)
    mesh.point_data["f"] = np.arange(len(mesh.points), dtype=float)
    with pytest.raises(ValueError, match="nothing"):
        mp.smooth_field(mesh, "nothing", 0.1)
    with pytest.raises(ValueError, match="length"):
        mp.smooth_field(mesh, "f", -1.0)
    src, dst = str(tmp_path / "s.vtu"), str(tmp_path / "o.vtu")
    mp.write(src, mesh)
    rc = main(["smooth-field", src, dst, "--array", "f", "--length", "0.3", "--json"])
    assert rc in (0, None)
    assert json.loads(capsys.readouterr().out)["converged"]
    assert "f:smooth" in mp.read(dst).point_data
