"""Tests for roadmap items 9.11–9.15 and 15.1–15.6."""

from __future__ import annotations

import numpy as np

import meshioplusplus as mio
from meshioplusplus import physicsnemo
from meshioplusplus.mcp import _tools as mcp_tools


def _unit_tet_mesh() -> mio.Mesh:
    return mio.box(dims=(2, 2, 2), cell_type="tetra")


def _unit_sphere_mesh(subdivisions: int = 2) -> mio.Mesh:
    return mio.icosa_sphere(radius=1.0, subdivisions=subdivisions)


# ---------------------------------------------------------------------------
# 9.11: FVM face-pair quality + quality_gate fvm:* rules
# ---------------------------------------------------------------------------


def test_fvm_quality_and_gate():
    m = mio.box(dims=(2, 2, 2))
    res = mio.compute_fvm_quality(m)
    assert res["num_interior_faces"] == 12
    assert "fvm:non_orthogonality" in res["cell_arrays"]
    assert "fvm:face_skewness" in res["cell_arrays"]
    assert "fvm:volume_ratio" in res["cell_arrays"]

    # Structured Cartesian box has ~0 non-orthogonality and ~0 skewness.
    non_ortho = res["cell_arrays"]["fvm:non_orthogonality"][0]
    skew = res["cell_arrays"]["fvm:face_skewness"][0]
    assert np.max(non_ortho) < 1e-6
    assert np.max(skew) < 1e-6

    attached = mio.attach_fvm_quality(m)
    assert "fvm:non_orthogonality" in attached.cell_data

    # quality_gate lazy fvm:* rule
    gate = mio.check_quality(m, require="fvm:non_orthogonality <= 70")
    assert gate["passed"] is True

    gate_fail = mio.check_quality(m, require="fvm:non_orthogonality <= -1")
    assert gate_fail["passed"] is False


# ---------------------------------------------------------------------------
# 9.12: Feature-preserving surface denoising
# ---------------------------------------------------------------------------


def test_denoise_surface_reduces_noise():
    sph = _unit_sphere_mesh(3)
    rng = np.random.default_rng(42)
    pts = sph.points.copy()
    radii_clean = np.linalg.norm(pts, axis=1)
    radial_noise = rng.normal(0.0, 0.02, size=(len(pts), 1))
    noisy_pts = pts + radial_noise * (pts / radii_clean[:, None])
    noisy = mio.Mesh(noisy_pts, sph.cells)

    for method in ("bilateral", "guided"):
        denoised, rep = mio.denoise_surface(
            noisy,
            method=method,
            normal_iterations=8,
            vertex_iterations=10,
            return_report=True,
        )
        assert rep["num_nodes_moved"] > 0
        err_before = np.std(np.linalg.norm(noisy.points, axis=1) - radii_clean)
        err_after = np.std(np.linalg.norm(denoised.points, axis=1) - radii_clean)
        assert err_after < err_before


# ---------------------------------------------------------------------------
# 9.13: Curve-skeleton / centerline extraction
# ---------------------------------------------------------------------------


def test_extract_skeleton():
    tor = mio.torus(
        major_radius=1.0, minor_radius=0.25, nu=20, nv=12, cell_type="triangle"
    )
    skel, rep = mio.extract_skeleton(tor, max_iterations=8, return_report=True)
    assert rep["num_nodes"] > 4
    assert any(cb.type == "line" for cb in skel.cells)
    assert "skeleton:radius" in skel.point_data


# ---------------------------------------------------------------------------
# 9.14: Shortest geodesic path & geodesic Voronoi
# ---------------------------------------------------------------------------


def test_shortest_path_and_geodesic_voronoi():
    sph = _unit_sphere_mesh(2)
    z = sph.points[:, 2]
    src = int(np.argmax(z))
    dst = int(np.argmin(z))

    for method in ("dijkstra", "heat_trace"):
        path_mesh, rep = mio.shortest_path(
            sph, source=src, target=dst, method=method, return_report=True
        )
        assert len(path_mesh.points) >= 2
        assert rep["total_length"] > 2.5  # great circle ~ pi

    vor = mio.geodesic_voronoi(sph, sources=[src, dst])
    assert "geodesic:nearest_source" in vor.point_data
    assert "geodesic:distance" in vor.point_data
    labels = set(np.asarray(vor.point_data["geodesic:nearest_source"]).tolist())
    assert labels == {0, 1}


# ---------------------------------------------------------------------------
# 9.15: FFD & RBF deformation
# ---------------------------------------------------------------------------


def test_deform_ffd_and_rbf():
    m = mio.box(dims=(2, 2, 2))

    # FFD rigid translation
    disp_ffd = np.zeros((2 * 2 * 2, 3), dtype=np.float64)
    disp_ffd[:, 0] = 0.25
    ffd_mesh, ffd_rep = mio.deform_ffd(
        m,
        displacements=disp_ffd,
        grid_dims=(2, 2, 2),
        degrees=(1, 1, 1),
        return_report=True,
    )
    assert np.allclose(ffd_mesh.points[:, 0], m.points[:, 0] + 0.25)
    assert ffd_rep["num_inverted"] == 0

    # RBF TPS + Wendland C2
    ctrl = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]],
        dtype=np.float64,
    )
    disp = np.zeros_like(ctrl)
    disp[:, 1] = 0.1
    for kernel, radius in (("tps", 0.0), ("wendland_c2", 3.0)):
        rbf_mesh, _ = mio.deform_rbf(
            m,
            control_points=ctrl,
            displacements=disp,
            kernel=kernel,
            radius=radius,
            return_report=True,
        )
        assert np.allclose(rbf_mesh.points[:, 1], m.points[:, 1] + 0.1, atol=1e-5)


# ---------------------------------------------------------------------------
# 15.1: Quadrature mode in data_integrate
# ---------------------------------------------------------------------------


def test_data_integrate_quadrature():
    m = mio.box(dims=(2, 2, 2), order=2)
    # Integrate x^2 over [0, 1]^3 -> exact value 1/3 on quadratic hexes
    m.point_data["x2"] = m.points[:, 0] ** 2
    res = mio.data_integrate(m, arrays=["quad:x2"])
    assert abs(float(res[0]["domain"]["total_per_component"][0]) - 1.0 / 3.0) < 1e-10


# ---------------------------------------------------------------------------
# 15.2: Tet mesh_operator & smooth_field
# ---------------------------------------------------------------------------


def test_tet_mesh_operator_and_smooth_field():
    tet = _unit_tet_mesh()
    L = mio.mesh_operator(tet, kind="laplacian", sparse=True)
    assert L.shape == (len(tet.points), len(tet.points))
    assert L.nnz > 0

    M = mio.mesh_operator(tet, kind="mass_lumped", sparse=True)
    assert abs(float(M.sum()) - 1.0) < 1e-9

    tet.point_data["u"] = tet.points[:, 0] ** 2 + np.sin(5.0 * tet.points[:, 1])
    out = mio.smooth_field(tet, array="u", length=0.2)
    assert "u" in out.point_data


# ---------------------------------------------------------------------------
# 15.3: High-order evaluate_cells, change_basis, critical_points, gradient, isosurface
# ---------------------------------------------------------------------------


def test_high_order_operations():
    base = mio.rectangle(dims=(2, 2), cell_type="triangle")
    quad = mio.convert_cells(base, mode="elevate")  # triangle6
    x, y = quad.points[:, 0], quad.points[:, 1]
    quad.point_data["f"] = x**2 + y**2 - x - y

    # evaluate_cells
    ev = mio.evaluate_cells(
        quad,
        cells=[0],
        reference_points=[[1.0 / 3.0, 1.0 / 3.0, 0.0]],
        fields=["f"],
        derivatives=True,
    )
    assert ev["num_evaluated"] == 1
    assert "f" in ev["values"]
    assert "f" in ev["gradients"]

    # gradient on quadratic mesh should be finite (not NaN)
    grad_mesh = mio.gradient(quad, array="f", location="cell")
    grad_key = next(iter(grad_mesh.cell_data.keys()))
    grad_vals = grad_mesh.cell_data[grad_key][0]
    assert np.all(np.isfinite(grad_vals))

    # critical_points finds interior minimum near (0.5, 0.5)
    cp_mesh, cp_rep = mio.critical_points(quad, array="f", return_report=True)
    assert cp_rep["num_minima"] >= 1
    dists = np.linalg.norm(cp_mesh.points[:, :2] - np.array([0.5, 0.5]), axis=1)
    assert np.min(dists) < 1e-4

    # isosurface on quadratic mesh
    iso = mio.isosurface(quad, array="f", isovalues=[-0.2])
    assert len(iso.points) > 0

    # change_basis Lagrange -> Bezier -> Lagrange
    pts_curve = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.5, 0.0]])
    lag_mesh = mio.Mesh(pts_curve, [("VTK_LAGRANGE_CURVE", np.array([[0, 1, 2]]))])
    lag_mesh.point_data["s"] = np.array([0.0, 1.0, 0.5])
    bez_mesh, bez_rep = mio.change_basis(lag_mesh, target="bezier", return_report=True)
    assert bez_rep["num_converted"] == 1
    back_mesh, back_rep = mio.change_basis(
        bez_mesh, target="lagrange", return_report=True
    )
    assert back_rep["num_converted"] == 1
    assert np.allclose(back_mesh.points, pts_curve, atol=1e-10)


# ---------------------------------------------------------------------------
# 15.4: Streamlines pipeline step
# ---------------------------------------------------------------------------


def test_streamlines_pipeline_step(tmp_path):
    m = mio.box(dims=(2, 2, 2))
    vel = np.zeros_like(m.points)
    vel[:, 0] = 1.0
    m.point_data["v"] = vel
    in_path = str(tmp_path / "flow.vtu")
    out_path = str(tmp_path / "lines.vtu")
    mio.write(in_path, m)

    mio.run_pipeline(
        {
            "Input": {"Path": in_path},
            "Operations": [
                {
                    "Op": "Streamlines",
                    "Array": "v",
                    "Seeds": 4,
                    "MaxLength": 0.5,
                }
            ],
            "Output": {"Path": out_path},
        }
    )
    out = mio.read(out_path)
    assert len(out.points) > 0


# ---------------------------------------------------------------------------
# 15.5: Inter-mesh transfer_operator + physicsnemo.transfer_matrix
# ---------------------------------------------------------------------------


def test_transfer_operator_and_physicsnemo():
    coarse = mio.rectangle(dims=(2, 2), cell_type="triangle")
    fine = mio.refine(coarse, levels=1, record_parent_ids=True)

    P = mio.transfer_operator(coarse, fine, mode="auto", sparse=True)
    assert P.shape == (len(fine.points), len(coarse.points))
    # Linear field f(x, y) = 2x + 3y transfers exactly
    fc = 2.0 * coarse.points[:, 0] + 3.0 * coarse.points[:, 1]
    ff_exact = 2.0 * fine.points[:, 0] + 3.0 * fine.points[:, 1]

    sp = physicsnemo.transfer_matrix(coarse, fine, format="scipy")
    ff_interp = sp @ fc
    assert np.allclose(ff_interp, ff_exact, atol=1e-10)

    # Transpose restriction matrix
    R = mio.transfer_operator(coarse, fine, transpose=True, sparse=True)
    assert R.shape == (len(coarse.points), len(fine.points))


# ---------------------------------------------------------------------------
# 15.6: Spectral shape descriptors (HKS / WKS)
# ---------------------------------------------------------------------------


def test_spectral_descriptors():
    sph = _unit_sphere_mesh(2)
    out, rep = mio.spectral_descriptors(
        sph, kind="hks", num_modes=12, num_scales=6, return_report=True
    )
    assert "spectral:hks" in out.point_data
    assert out.point_data["spectral:hks"].shape == (len(sph.points), 6)
    assert rep["reused_eigenmodes"] is False

    # Reusing attached eigenmodes for WKS
    out2, rep2 = mio.spectral_descriptors(
        out, kind="wks", num_modes=12, num_scales=6, return_report=True
    )
    assert "spectral:wks" in out2.point_data
    assert out2.point_data["spectral:wks"].shape == (len(sph.points), 6)
    assert rep2["reused_eigenmodes"] is True


# ---------------------------------------------------------------------------
# MCP tools smoke check
# ---------------------------------------------------------------------------


def test_mcp_new_tools_smoke(tmp_path):
    sph_in = str(tmp_path / "sph.vtu")
    mio.write(sph_in, _unit_sphere_mesh(2))

    # fvm_quality on box
    box_in = str(tmp_path / "box.vtu")
    mio.write(box_in, mio.box(dims=(2, 2, 2)))
    out_fvm = mcp_tools.tool_fvm_quality(box_in, output_path=str(tmp_path / "fvm.vtu"))
    assert out_fvm["num_interior_faces"] == 12

    # shortest_path
    out_sp = mcp_tools.tool_shortest_path(
        sph_in, str(tmp_path / "sp.vtu"), source=0, target=10
    )
    assert out_sp["total_length"] > 0.0

    # spectral_descriptors
    out_sd = mcp_tools.tool_spectral_descriptors(
        sph_in, str(tmp_path / "sd.vtu"), num_modes=8, num_scales=4
    )
    assert out_sd["num_modes_used"] >= 1
