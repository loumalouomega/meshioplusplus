"""Tests for surface remeshing (the ``remesh`` operation)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import remesh

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(
    _core is None, reason="remesh has no pure-Python fallback"
)


def _icosahedron():
    t = (1.0 + np.sqrt(5.0)) / 2.0
    pts = np.array(
        [
            [-1, t, 0],
            [1, t, 0],
            [-1, -t, 0],
            [1, -t, 0],
            [0, -1, t],
            [0, 1, t],
            [0, -1, -t],
            [0, 1, -t],
            [t, 0, -1],
            [t, 0, 1],
            [-t, 0, -1],
            [-t, 0, 1],
        ],
        dtype=float,
    )
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    faces = np.array(
        [
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ],
        dtype=np.int64,
    )
    return meshioplusplus.Mesh(pts, [("triangle", faces)])


@needs_core
def test_produces_the_requested_number_of_clusters():
    out, report = remesh(_icosahedron(), 150, return_report=True)
    assert out.points.shape[0] == 150
    assert report["num_clusters"] == 150
    assert out.cells[0].type == "triangle"
    assert out.cells[0].data.shape[0] > 0


@needs_core
def test_output_is_watertight():
    out = remesh(_icosahedron(), 150)
    conn = out.cells[0].data
    edges = {}
    for tri in conn:
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            key = (a, b) if a < b else (b, a)
            edges[key] = edges.get(key, 0) + 1
    assert set(edges.values()) == {2}


@needs_core
def test_quadric_metric_is_accepted_and_deterministic():
    a, ra = remesh(_icosahedron(), 80, metric="quadric", return_report=True)
    b, rb = remesh(_icosahedron(), 80, metric="quadric", return_report=True)
    assert a.points.shape == b.points.shape
    assert ra["num_iterations"] == rb["num_iterations"]
    np.testing.assert_array_equal(a.points, b.points)


@needs_core
def test_unknown_metric_raises():
    with pytest.raises(Exception):
        remesh(_icosahedron(), 50, metric="bogus")


@needs_core
def test_too_few_clusters_raises():
    with pytest.raises(Exception):
        remesh(_icosahedron(), 3)


@needs_core
def test_rejects_volume_cells_by_name():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    conn = np.array([[0, 1, 2, 3]], dtype=np.int64)
    mesh = meshioplusplus.Mesh(pts, [("tetra", conn)])
    with pytest.raises(Exception, match="extract_surface"):
        remesh(mesh, 10)


@needs_core
def test_point_data_is_dropped_field_data_is_carried():
    mesh = _icosahedron()
    mesh.point_data["temperature"] = np.full(mesh.points.shape[0], 3.0)
    mesh.field_data["solver"] = np.array([1.5])
    out = remesh(mesh, 40)
    assert "temperature" not in out.point_data
    assert "solver" in out.field_data


@needs_core
def test_subdivide_zero_disables_refinement():
    out, report = remesh(_icosahedron(), 10, subdivide=0, return_report=True)
    assert report["subdivide_applied"] == 0


def _gaussian_bump(n=40, half=1.0, amplitude=0.6, sigma=0.25):
    i, j = np.meshgrid(np.arange(n + 1), np.arange(n + 1))
    x = -half + 2.0 * half * i / n
    y = -half + 2.0 * half * j / n
    z = amplitude * np.exp(-(x * x + y * y) / (sigma * sigma))
    pts = np.column_stack([x.ravel(), y.ravel(), z.ravel()])

    def vid(ii, jj):
        return jj * (n + 1) + ii

    faces = []
    for jj in range(n):
        for ii in range(n):
            b = vid(ii, jj)
            faces.append([b, b + 1, b + n + 2])
            faces.append([b, b + n + 2, b + n + 1])
    return meshioplusplus.Mesh(pts, [("triangle", np.array(faces, dtype=np.int64))])


@needs_core
def test_curvature_gradation_concentrates_clusters_near_high_curvature():
    bump = _gaussian_bump()

    def count_near_center(mesh, radius):
        x, y = mesh.points[:, 0], mesh.points[:, 1]
        return int(np.count_nonzero(x * x + y * y < radius * radius))

    uniform = remesh(bump, 150, subdivide=0, gradation=0.0)
    curved = remesh(bump, 150, subdivide=0, gradation=2.0)

    radius = 0.4
    n_uniform = count_near_center(uniform, radius)
    n_curved = count_near_center(curved, radius)
    assert n_curved > n_uniform, f"{n_curved} vs {n_uniform} (uniform)"


def _open_square_patch(n=24, half=1.0):
    i, j = np.meshgrid(np.arange(n + 1), np.arange(n + 1))
    x = -half + 2.0 * half * i / n
    y = -half + 2.0 * half * j / n
    pts = np.column_stack([x.ravel(), y.ravel(), np.zeros(x.size)])

    def vid(ii, jj):
        return jj * (n + 1) + ii

    faces = []
    for jj in range(n):
        for ii in range(n):
            b = vid(ii, jj)
            faces.append([b, b + 1, b + n + 2])
            faces.append([b, b + n + 2, b + n + 1])
    return meshioplusplus.Mesh(pts, [("triangle", np.array(faces, dtype=np.int64))])


@needs_core
def test_preserves_boundary_of_an_open_patch():
    patch = _open_square_patch()
    out, report = remesh(patch, 80, subdivide=0, return_report=True)

    assert (
        len(out.cells) == 2
    ), "an open input should leave a second, boundary line block"
    assert out.cells[0].type == "triangle"
    assert out.cells[1].type == "line"
    assert out.cells[1].data.shape[0] > 0
    assert report["num_non_manifold_vertices"] >= 0

    total_length = 0.0
    for a, b in out.cells[1].data:
        total_length += float(np.linalg.norm(out.points[a] - out.points[b]))
    true_perimeter = 4.0 * 2.0 * 1.0
    assert true_perimeter * 0.5 < total_length < true_perimeter * 1.5


@needs_core
def test_no_preserve_boundary_still_runs():
    patch = _open_square_patch()
    out = remesh(patch, 80, subdivide=0, preserve_boundary=False)
    assert out.cells[0].type == "triangle"


def _cylinder_mesh(n_circ=24, n_axial=48, radius=1.0, height=8.0):
    theta = 2.0 * np.pi * np.arange(n_circ) / n_circ
    ring_x, ring_y = radius * np.cos(theta), radius * np.sin(theta)
    j = np.arange(n_axial + 1)
    x = np.tile(ring_x, n_axial + 1)
    y = np.tile(ring_y, n_axial + 1)
    z = np.repeat(height * j / n_axial, n_circ)
    pts = np.column_stack([x, y, z])

    def vid(i, jj):
        return jj * n_circ + (i % n_circ)

    faces = []
    for jj in range(n_axial):
        for i in range(n_circ):
            v0, v1, v2, v3 = (
                vid(i, jj),
                vid(i + 1, jj),
                vid(i + 1, jj + 1),
                vid(i, jj + 1),
            )
            faces.append([v0, v1, v2])
            faces.append([v0, v2, v3])
    return meshioplusplus.Mesh(pts, [("triangle", np.array(faces, dtype=np.int64))])


def _axial_over_circumferential_ratio(mesh, radius):
    pts = mesh.points
    conn = mesh.cells[0].data
    sum_axial = sum_arc = 0.0
    for tri in conn:
        for k in range(3):
            a, b = tri[k], tri[(k + 1) % 3]
            xa, ya, za = pts[a]
            xb, yb, zb = pts[b]
            dtheta = np.arctan2(yb, xb) - np.arctan2(ya, xa)
            dtheta = (dtheta + np.pi) % (2 * np.pi) - np.pi
            sum_arc += abs(dtheta) * radius
            sum_axial += abs(za - zb)
    return sum_axial / sum_arc


@needs_core
def test_anisotropic_metric_elongates_along_the_low_curvature_axis():
    radius = 1.0
    cyl = _cylinder_mesh(radius=radius)

    iso = remesh(cyl, 200, subdivide=0, preserve_boundary=False)
    iso_ratio = _axial_over_circumferential_ratio(iso, radius)

    aniso = remesh(
        cyl,
        200,
        subdivide=0,
        preserve_boundary=False,
        metric="anisotropic",
        max_anisotropy=8.0,
    )
    aniso_ratio = _axial_over_circumferential_ratio(aniso, radius)

    assert aniso_ratio > iso_ratio * 1.3, f"{aniso_ratio} vs isotropic {iso_ratio}"


@needs_core
def test_anisotropic_metric_keeps_clusters_on_a_curved_surface():
    amplitude, sigma = 0.5, 0.3
    bump = _gaussian_bump(30, 1.0, amplitude, sigma)

    out = remesh(
        bump,
        150,
        subdivide=0,
        metric="anisotropic",
        max_anisotropy=6.0,
        preserve_boundary=True,
    )
    x, y, z = out.points[:, 0], out.points[:, 1], out.points[:, 2]
    target = amplitude * np.exp(-(x * x + y * y) / (sigma * sigma))
    dev = float(np.max(np.abs(z - target)))
    assert (
        dev < amplitude * 0.5
    ), f"max vertical deviation from the true bump surface: {dev}"


@needs_core
def test_max_anisotropy_validation():
    ico = _icosahedron()
    with pytest.raises(Exception):
        remesh(ico, 60, metric="anisotropic", max_anisotropy=0.5)
    with pytest.raises(Exception):
        remesh(ico, 60, metric="isotropic", max_anisotropy=5.0)
    # the default is always accepted, regardless of metric
    remesh(ico, 60, metric="isotropic", max_anisotropy=4.0)
