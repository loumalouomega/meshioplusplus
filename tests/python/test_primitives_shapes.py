"""The second tier of parametric generators -- the curves and point clouds (line,
circle, helix, sphere_points, box_points), the surfaces (tetra_sphere,
octa_sphere, cube_sphere, geo_sphere, torus, moebius, tube, cone, capsule,
box_surface, platonic) and the solids (ball, spherical_shell, pipe, fichera,
solid_torus) -- and the full-Lagrange (``complete``) and VTK Lagrange
(``order`` 3 and up) cells every generator builds.
"""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _primitives as prim

NEW = (
    "line",
    "circle",
    "helix",
    "sphere_points",
    "box_points",
    "tetra_sphere",
    "octa_sphere",
    "cube_sphere",
    "geo_sphere",
    "torus",
    "moebius",
    "tube",
    "cone",
    "capsule",
    "box_surface",
    "platonic",
    "ball",
    "spherical_shell",
    "pipe",
    "fichera",
    "solid_torus",
)

# Index arithmetic and sqrt only: the twin is byte-identical to the core.
POLYNOMIAL = {
    "line",
    "box_points",
    "tetra_sphere",
    "octa_sphere",
    "cube_sphere",
    "geo_sphere",
    "box_surface",
    "platonic",
    "ball",
    "fichera",
}


def _twin(monkeypatch, shape, **kwargs):
    """The numpy twin, through the public function with the core taken away."""

    def no_core(op, *args):
        raise ImportError("no core")

    with monkeypatch.context() as m:
        m.setattr(prim, "_core_call", no_core)
        return getattr(mp, shape)(**kwargs)


def _same(a, b, exact):
    assert [(c.type, np.asarray(c.data).shape) for c in a.cells] == [
        (c.type, np.asarray(c.data).shape) for c in b.cells
    ]
    for ca, cb in zip(a.cells, b.cells):
        assert np.array_equal(np.asarray(ca.data), np.asarray(cb.data))
    if exact:
        assert np.asarray(a.points).tobytes() == np.asarray(b.points).tobytes()
    else:
        np.testing.assert_allclose(a.points, b.points, rtol=1e-13, atol=1e-14)
    ra = {(r.name, r.kind): np.asarray(r.entries) for r in a.regions}
    rb = {(r.name, r.kind): np.asarray(r.entries) for r in b.regions}
    assert ra.keys() == rb.keys()
    for key in ra:
        assert np.array_equal(ra[key], rb[key]), key


def _region(mesh, name):
    return next(r for r in mesh.regions if r.name == name)


def _triangles(mesh):
    out = []
    for block in mesh.cells:
        d = np.asarray(block.data)
        if block.type.startswith("triangle"):
            out.append(d[:, :3])
        elif block.type.startswith("quad"):
            out += [d[:, [0, 1, 2]], d[:, [0, 2, 3]]]
        elif block.type == "polygon":
            out += [d[:, [0, k, k + 1]] for k in range(1, d.shape[1] - 1)]
    return np.concatenate(out)


def _enclosed_volume(mesh):
    """Signed volume a closed surface encloses (positive when outward)."""
    t = _triangles(mesh)
    p = mesh.points
    return np.einsum("ij,ij->i", p[t[:, 0]], np.cross(p[t[:, 1]], p[t[:, 2]])).sum() / 6


def _euler(mesh):
    """V - E + F over the corner nodes of the surface cells."""
    corners = {"triangle": 3, "quad": 4}
    verts, edges, faces = set(), set(), 0
    for block in mesh.cells:
        d = np.asarray(block.data)
        k = corners.get(block.type.rstrip("0123456789"), d.shape[1])
        faces += d.shape[0]
        for row in d[:, :k]:
            verts.update(row.tolist())
            for i in range(k):
                a, b = int(row[i]), int(row[(i + 1) % k])
                edges.add((min(a, b), max(a, b)))
    return len(verts) - len(edges) + faces


_CORNER_TETS = [(0, 1, 3, 4), (1, 2, 0, 5), (2, 3, 1, 6), (3, 0, 2, 7)]
_CORNER_TETS += [(4, 7, 5, 0), (5, 4, 6, 1), (6, 5, 7, 2), (7, 6, 4, 3)]


def _min_jacobian(mesh):
    """The smallest corner Jacobian of the tetrahedra and hexahedra."""
    p = mesh.points
    worst = np.inf
    for block in mesh.cells:
        d = np.asarray(block.data)
        if block.type.startswith("tetra"):
            corners = [(0, 1, 2, 3)]
        elif block.type.startswith("hexahedron"):
            corners = _CORNER_TETS
        else:
            continue
        for a, b, c, e in corners:
            u, v, w = (
                p[d[:, b]] - p[d[:, a]],
                p[d[:, c]] - p[d[:, a]],
                p[d[:, e]] - p[d[:, a]],
            )
            worst = min(worst, float(np.einsum("ij,ij->i", u, np.cross(v, w)).min()))
    return worst


def _welded(mesh):
    """No two points coincide (seams and poles are shared, not duplicated)."""
    p = np.round(np.asarray(mesh.points), 9)
    return np.unique(p, axis=0).shape[0] == p.shape[0]


# --------------------------------------------------------------------------- #
# contract                                                                    #
# --------------------------------------------------------------------------- #


def test_public_names():
    for name in NEW:
        assert name in mp.__all__
        assert callable(getattr(mp, name))


@pytest.mark.parametrize(
    "call",
    [
        lambda: mp.line(n=0),
        lambda: mp.circle(n=0),
        lambda: mp.helix(n=0),
        lambda: mp.sphere_points(n=0),
        lambda: mp.box_points(n=0),
        lambda: mp.cube_sphere(n=0),
        lambda: mp.geo_sphere(n=0),
        lambda: mp.torus(nu=0),
        lambda: mp.moebius(nv=0),
        lambda: mp.tube(nz=0),
        lambda: mp.cone(ntheta=0),
        lambda: mp.capsule(ncap=0),
        lambda: mp.box_surface(dims=(2, 0, 2)),
        lambda: mp.ball(n=0),
        lambda: mp.spherical_shell(nr=0),
        lambda: mp.pipe(nz=0),
        lambda: mp.fichera(n=0),
        lambda: mp.solid_torus(nu=0),
    ],
)
def test_zero_count_is_an_empty_mesh(call):
    out = call()
    assert out.points.shape == (0, 3)
    assert len(out.cells) == 0


@pytest.mark.parametrize(
    "call, text",
    [
        (lambda: mp.line(start=(1, 1, 1), end=(1, 1, 1)), "start and end must differ"),
        (lambda: mp.circle(n=2), "n must be at least 3"),
        (lambda: mp.helix(pitch=0.0), "pitch must be finite and nonzero"),
        (lambda: mp.helix(n=4, turns=2.0), "n must exceed 2 \\* turns"),
        (lambda: mp.sphere_points(radius=0.0), "radius must be positive"),
        (lambda: mp.box_points(size=(1, 0, 1)), "size must be positive"),
        (lambda: mp.tetra_sphere(subdivisions=-1), "must not be negative"),
        (lambda: mp.cube_sphere(cell_type="hexahedron"), "unknown cell_type"),
        (lambda: mp.torus(minor_radius=2.0), "minor_radius must be less than"),
        (lambda: mp.torus(nv=2), "nu and nv must be at least 3"),
        (lambda: mp.moebius(width=3.0), "width must be less than twice"),
        (lambda: mp.moebius(num_twists=-1), "must not be negative"),
        (lambda: mp.tube(ntheta=2), "ntheta must be at least 3"),
        (lambda: mp.cone(height=0.0), "height must be positive"),
        (lambda: mp.capsule(radius=-1.0), "radius must be positive"),
        (lambda: mp.platonic("sphere"), "unknown name 'sphere'"),
        (lambda: mp.ball(cell_type="wedge"), "unknown cell_type 'wedge'"),
        (lambda: mp.spherical_shell(outer_radius=0.25), "must exceed inner_radius"),
        (lambda: mp.spherical_shell(ratio=0.0), "ratio must be positive"),
        (lambda: mp.pipe(ntheta=2), "ntheta must be at least 3"),
        (lambda: mp.fichera(size=0.0), "size must be positive"),
        (lambda: mp.solid_torus(minor_radius=1.5), "minor_radius must be less than"),
        (lambda: mp.torus(order=11), "order must be between 1 and 10"),
    ],
)
def test_bad_arguments_raise_by_name(call, text):
    with pytest.raises(ValueError, match=r"meshio\+\+: .*" + text):
        call()


@pytest.mark.parametrize(
    "call",
    [
        lambda: mp.line(n=100, max_cells=50),
        lambda: mp.sphere_points(n=100, max_cells=50),
        lambda: mp.box_points(n=100, max_cells=50),
        lambda: mp.octa_sphere(subdivisions=3, max_cells=50),
        lambda: mp.geo_sphere(n=4, max_cells=50),
        lambda: mp.torus(max_cells=50),
        lambda: mp.cone(max_cells=50),
        lambda: mp.box_surface(dims=(4, 4, 4), max_cells=50),
        lambda: mp.ball(n=3, max_cells=50),
        lambda: mp.pipe(max_cells=50),
        lambda: mp.fichera(n=3, max_cells=50),
        lambda: mp.solid_torus(max_cells=50),
    ],
)
def test_max_cells_budget(call):
    with pytest.raises(ValueError, match="above the limit of 50"):
        call()


# --------------------------------------------------------------------------- #
# twin vs core                                                                #
# --------------------------------------------------------------------------- #

_CASES = [
    ("line", dict(n=5, start=(0.1, 0.2, 0.3), end=(1.0, -2.0, 3.0))),
    ("circle", dict(radius=1.5, n=7, center=(1.0, 2.0, 3.0), order=2)),
    ("helix", dict(radius=0.5, pitch=0.3, turns=1.5, n=13, order=3)),
    ("sphere_points", dict(n=37, radius=2.0, center=(1.0, 0.0, 0.0))),
    ("box_points", dict(n=40, size=(1.0, 2.0, 3.0), origin=(-1, 0, 1), seed=12345)),
    ("tetra_sphere", dict(radius=2.0, subdivisions=2, order=2)),
    ("octa_sphere", dict(subdivisions=2, center=(1.0, 2.0, 3.0))),
    ("cube_sphere", dict(n=3, order=2, complete=True)),
    ("cube_sphere", dict(n=2, cell_type="triangle", order=3)),
    ("geo_sphere", dict(n=3, order=2)),
    ("geo_sphere", dict(n=5)),
    ("torus", dict(nu=8, nv=5, order=2, complete=True)),
    ("torus", dict(nu=6, nv=4, cell_type="triangle", order=3)),
    ("moebius", dict(nu=9, nv=2, order=2)),
    ("moebius", dict(nu=9, nv=3, num_twists=2, cell_type="triangle", order=3)),
    ("tube", dict(ntheta=7, nz=3, order=2, complete=True)),
    ("cone", dict(ntheta=7, nz=3, order=2)),
    ("cone", dict(ntheta=5, nz=1, cell_type="triangle", order=3)),
    ("capsule", dict(ntheta=7, nz=2, ncap=3, order=2, complete=True)),
    ("box_surface", dict(dims=(2, 3, 1), size=(1.0, 2.0, 0.5), order=2)),
    ("box_surface", dict(dims=(1, 2, 2), cell_type="triangle", order=3)),
    ("platonic", dict(name="dodecahedron", radius=2.0)),
    ("platonic", dict(name="tetrahedron", center=(1.0, 1.0, 1.0))),
    ("ball", dict(n=2, order=2, complete=True)),
    ("ball", dict(n=1, cell_type="tetra", order=3)),
    ("spherical_shell", dict(n=2, nr=2, ratio=1.5, order=2, complete=True)),
    ("spherical_shell", dict(n=2, nr=1)),
    ("pipe", dict(ntheta=6, nr=2, nz=2, ratio=1.3, order=2, complete=True)),
    ("fichera", dict(n=1, order=3)),
    ("fichera", dict(n=2, cell_type="tetra", order=2)),
    ("solid_torus", dict(n=2, nu=6, order=2, complete=True)),
]


@pytest.mark.parametrize("name, kwargs", _CASES)
def test_twin_matches_core(monkeypatch, name, kwargs):
    exact = name in POLYNOMIAL or (name == "spherical_shell" and "ratio" not in kwargs)
    _same(getattr(mp, name)(**kwargs), _twin(monkeypatch, name, **kwargs), exact)


# --------------------------------------------------------------------------- #
# curves and point clouds                                                     #
# --------------------------------------------------------------------------- #


def test_line_runs_from_start_to_end():
    m = mp.line(n=4, start=(1.0, 0.0, 0.0), end=(1.0, 2.0, 4.0))
    np.testing.assert_array_equal(m.points[0], [1.0, 0.0, 0.0])
    np.testing.assert_array_equal(m.points[-1], [1.0, 2.0, 4.0])
    assert m.cells[0].type == "line"
    assert _region(m, "start").kind == "point"
    assert np.asarray(_region(m, "end").entries).tolist() == [4]


def test_circle_is_closed_and_on_the_circle():
    m = mp.circle(radius=2.0, n=9, order=2)
    assert m.cells[0].type == "line3"
    conn = np.asarray(m.cells[0].data)
    assert conn[-1, 1] == conn[0, 0]
    np.testing.assert_allclose(np.linalg.norm(m.points[:, :2], axis=1), 2.0)


def test_helix_rises_by_the_pitch():
    m = mp.helix(radius=1.0, pitch=0.5, turns=3.0, n=30, order=2)
    np.testing.assert_allclose(m.points[30], [1.0, 0.0, 1.5], atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(m.points[:, :2], axis=1), 1.0)
    # The mid-edge nodes are on the helix: z follows the angle.
    ang = np.unwrap(np.arctan2(m.points[31:, 1], m.points[31:, 0]))
    np.testing.assert_allclose(
        m.points[31:, 2], 0.5 * (ang - ang[0]) / (2 * np.pi) + 0.025
    )


def test_sphere_points_are_on_the_sphere():
    m = mp.sphere_points(n=200, radius=3.0)
    np.testing.assert_allclose(np.linalg.norm(m.points, axis=1), 3.0)
    assert m.cells[0].type == "vertex"
    assert abs(m.points[:, 2].mean()) < 1e-12


def test_box_points_follow_splitmix64():
    m = mp.box_points(n=2, size=(1.0, 1.0, 1.0), seed=0)
    # splitmix64 seeded with 0 first returns 0xE220A8397B1DCDAF.
    assert m.points[0, 0] == (0xE220A8397B1DCDAF >> 11) * 2.0**-53
    other = mp.box_points(n=2, seed=1)
    assert not np.array_equal(m.points, other.points)
    assert mp.box_points(n=50, seed=7).points.tobytes() == (
        mp.box_points(n=50, seed=7).points.tobytes()
    )


# --------------------------------------------------------------------------- #
# surfaces                                                                    #
# --------------------------------------------------------------------------- #

_CLOSED = [
    (lambda: mp.tetra_sphere(), 4 / 3 * np.pi, 0.07),
    (lambda: mp.octa_sphere(), 4 / 3 * np.pi, 0.1),
    (lambda: mp.cube_sphere(n=6), 4 / 3 * np.pi, 0.04),
    (lambda: mp.cube_sphere(n=4, cell_type="triangle"), 4 / 3 * np.pi, 0.08),
    (lambda: mp.geo_sphere(n=5), 4 / 3 * np.pi, 0.03),
    (lambda: mp.torus(), 2 * np.pi**2 * 0.25**2, 0.04),
    (lambda: mp.torus(cell_type="triangle"), 2 * np.pi**2 * 0.25**2, 0.04),
    (lambda: mp.cone(ntheta=64), np.pi / 3, 0.01),
    (lambda: mp.cone(ntheta=64, cell_type="triangle"), np.pi / 3, 0.01),
    (lambda: mp.capsule(ntheta=64, ncap=16), np.pi * 0.25 + np.pi / 6, 0.01),
    (lambda: mp.box_surface(dims=(2, 3, 4), size=(1, 2, 3)), 6.0, 1e-12),
    (lambda: mp.box_surface(dims=(1, 1, 2), cell_type="triangle"), 1.0, 1e-12),
]


@pytest.mark.parametrize("call, volume, rtol", _CLOSED)
def test_closed_surfaces_are_outward_welded_spheres(call, volume, rtol):
    m = call()
    assert _welded(m)
    np.testing.assert_allclose(_enclosed_volume(m), volume, rtol=rtol)


@pytest.mark.parametrize(
    "call, chi",
    [
        (lambda: mp.geo_sphere(n=3), 2),
        (lambda: mp.cube_sphere(n=3), 2),
        (lambda: mp.capsule(), 2),
        (lambda: mp.cone(), 2),
        (lambda: mp.box_surface(dims=(2, 2, 2)), 2),
        (lambda: mp.torus(), 0),
        (lambda: mp.tube(), 0),
        (lambda: mp.moebius(), 0),
        (lambda: mp.moebius(num_twists=2), 0),
    ],
)
def test_euler_characteristic(call, chi):
    assert _euler(call()) == chi


@pytest.mark.parametrize("n", [1, 2, 4])
def test_geo_sphere_counts(n):
    m = mp.geo_sphere(n=n)
    assert len(m.cells[0].data) == 20 * n * n
    assert m.points.shape[0] == 10 * n * n + 2
    np.testing.assert_allclose(np.linalg.norm(m.points, axis=1), 1.0)


def test_moebius_has_no_consistent_orientation():
    _, report = mp.compute_normals(mp.moebius(nu=12, nv=1), return_report=True)
    assert report["quality"]["inconsistent_pairs"] > 0
    _, report = mp.compute_normals(
        mp.moebius(nu=12, nv=1, num_twists=2), return_report=True
    )
    assert report["quality"]["inconsistent_pairs"] == 0


def test_moebius_boundary_is_one_loop_for_an_odd_twist():
    m = mp.moebius(nu=10, nv=2)
    assert len(_region(m, "boundary").entries) == 2 * 10
    m = mp.moebius(nu=10, nv=2, num_twists=2)
    assert len(_region(m, "boundary").entries) == 2 * 10


def test_tube_rims_and_cone_regions():
    m = mp.tube(ntheta=8, nz=3)
    assert len(_region(m, "bottom").entries) == 8
    assert len(_region(m, "top").entries) == 8
    m = mp.cone(ntheta=8, nz=3)
    side, base = _region(m, "side"), _region(m, "base")
    assert side.kind == "cell" and base.kind == "cell"
    assert len(side.entries) + len(base.entries) == sum(len(c.data) for c in m.cells)
    assert len(side.entries) == len(base.entries) == 8 * 3


def test_box_surface_regions_are_its_faces():
    m = mp.box_surface(dims=(2, 3, 4), size=(1.0, 2.0, 3.0))
    sizes = {r.name: len(r.entries) for r in m.regions}
    assert sizes == {"xmin": 12, "xmax": 12, "ymin": 8, "ymax": 8, "zmin": 6, "zmax": 6}
    p = m.points
    for name, axis, value in [("xmin", 0, 0.0), ("ymax", 1, 2.0), ("zmax", 2, 3.0)]:
        cells = np.asarray(m.cells[0].data)[np.asarray(_region(m, name).entries)]
        np.testing.assert_array_equal(p[cells.reshape(-1), axis], value)


@pytest.mark.parametrize(
    "name, faces, corners",
    [
        ("tetrahedron", 4, 4),
        ("cube", 6, 8),
        ("octahedron", 8, 6),
        ("dodecahedron", 12, 20),
        ("icosahedron", 20, 12),
    ],
)
def test_platonic_solids(name, faces, corners):
    m = mp.platonic(name, radius=2.0)
    assert len(m.cells[0].data) == faces
    assert m.points.shape[0] == corners
    np.testing.assert_allclose(np.linalg.norm(m.points, axis=1), 2.0)
    assert _euler(m) == 2
    assert _enclosed_volume(m) > 0
    edge = np.linalg.norm(
        m.points[np.asarray(m.cells[0].data)[:, 0]]
        - m.points[np.asarray(m.cells[0].data)[:, 1]],
        axis=1,
    )
    np.testing.assert_allclose(edge, edge[0])


# --------------------------------------------------------------------------- #
# solids                                                                      #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "call",
    [
        lambda: mp.ball(n=3),
        lambda: mp.ball(n=2, cell_type="tetra"),
        lambda: mp.spherical_shell(ratio=2.0),
        lambda: mp.spherical_shell(ratio=0.5, nr=3),
        lambda: mp.pipe(ratio=0.5),
        lambda: mp.fichera(),
        lambda: mp.fichera(cell_type="tetra"),
        lambda: mp.solid_torus(),
    ],
)
def test_solids_are_positive_and_welded(call):
    m = call()
    assert _min_jacobian(m) > 0
    assert _welded(m)


def test_ball_surface_is_the_sphere():
    m = mp.ball(radius=2.0, n=2)
    hexes = np.asarray(m.cells[0].data)
    from meshioplusplus._skin import _CELL_FACES

    faces = _CELL_FACES["hexahedron"]
    surface = np.asarray(_region(m, "surface").entries)
    assert len(surface) == 6 * 4 * 4
    for cell, facet in surface:
        nodes = hexes[cell][list(faces[facet][2])]
        np.testing.assert_allclose(np.linalg.norm(m.points[nodes], axis=1), 2.0)


def test_spherical_shell_grading():
    m = mp.spherical_shell(inner_radius=1.0, outer_radius=2.0, n=1, nr=3, ratio=2.0)
    # Layers thicken by the ratio: 1/7, 2/7, 4/7 of the wall.
    radii = np.unique(np.round(np.linalg.norm(m.points, axis=1), 12))
    np.testing.assert_allclose(radii, [1.0, 1 + 1 / 7, 1 + 3 / 7, 2.0])


def test_pipe_regions():
    m = mp.pipe(ntheta=8, nr=2, nz=3)
    sizes = {r.name: len(r.entries) for r in m.regions}
    assert sizes == {"inner": 24, "outer": 24, "bottom": 16, "top": 16}


def test_fichera_loses_an_octant():
    m = mp.fichera(n=2, size=2.0)
    assert len(m.cells[0].data) == 7 * 8
    centres = m.points[np.asarray(m.cells[0].data)].mean(axis=1)
    assert not ((centres > 0).all(axis=1)).any()
    sizes = {r.name: len(r.entries) for r in m.regions}
    assert sizes["notch_x"] == sizes["notch_y"] == sizes["notch_z"] == 4
    assert sizes["xmin"] == 16 and sizes["xmax"] == 12


# --------------------------------------------------------------------------- #
# complete and higher orders                                                  #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "call, types",
    [
        (lambda **k: mp.rectangle(dims=(2, 2), **k), ["quad9"]),
        (lambda **k: mp.box(dims=(2, 1, 1), **k), ["hexahedron27"]),
        (lambda **k: mp.box(dims=(1, 1, 1), cell_type="wedge", **k), ["wedge18"]),
        (lambda **k: mp.cylinder(nxy=2, nz=1, cell_type="wedge", **k), ["wedge18"]),
        (
            lambda **k: mp.uv_sphere(nu=6, nv=3, cell_type="quad", **k),
            ["quad9", "triangle6"],
        ),
        (lambda **k: mp.disk(n=2, cell_type="triangle", **k), ["triangle6"]),
    ],
)
def test_complete_builds_full_lagrange_types(call, types):
    assert [c.type for c in call(order=2, complete=True).cells] == types


def test_complete_centres_are_shared_and_on_the_geometry():
    m = mp.cube_sphere(n=2, order=2, complete=True)
    np.testing.assert_allclose(np.linalg.norm(m.points, axis=1), 1.0)
    m = mp.box(dims=(2, 2, 2), order=2, complete=True)
    # 27 lattice points per unit, shared: a 5 x 5 x 5 lattice.
    assert m.points.shape[0] == 125
    np.testing.assert_allclose(np.unique(m.points), [0.0, 0.25, 0.5, 0.75, 1.0])


@pytest.mark.parametrize("order", [3, 4, 6])
def test_vtk_lagrange_nodes_follow_the_vtk_lattice(order):
    from meshioplusplus._lagrange import lattice_weights, vtk_lattice

    for call, shape in [
        (lambda: mp.box(dims=(2, 1, 1), order=order), "hexahedron"),
        (lambda: mp.box(cell_type="tetra", order=order), "tetra"),
        (lambda: mp.box(cell_type="wedge", order=order), "wedge"),
        (lambda: mp.rectangle(dims=(2, 2), order=order), "quad"),
        (lambda: mp.reference_triangle(n=2, order=order), "triangle"),
        (lambda: mp.line(n=3, order=order), "line"),
    ]:
        m = call()
        conn = np.asarray(m.cells[0].data)
        weights = [
            lattice_weights(shape, order, ijk) for ijk in vtk_lattice(shape, order)
        ]
        corners = len([w for w in weights if len(w) == 1])
        for row in conn:
            expect = np.array(
                [sum(v * m.points[row[c]] for c, v in w) / order**3 for w in weights]
            )
            np.testing.assert_allclose(m.points[row], expect, atol=1e-14)
            assert len(set(row[:corners].tolist())) == corners


@pytest.mark.parametrize("order", [3, 5])
def test_vtk_lagrange_shares_nodes_and_curves(order):
    m = mp.geo_sphere(n=2, order=order)
    assert m.cells[0].type == "VTK_LAGRANGE_TRIANGLE"
    np.testing.assert_allclose(np.linalg.norm(m.points, axis=1), 1.0)
    # Shared: an order-p triangle mesh of a closed surface has F p^2 / 2 + 2 nodes.
    faces = len(m.cells[0].data)
    assert m.points.shape[0] == faces * order * order // 2 + 2
    m = mp.ball(n=1, order=order)
    assert m.cells[0].type == "VTK_LAGRANGE_HEXAHEDRON"
    assert m.points.shape[0] == (2 * order + 1) ** 3
    assert np.linalg.norm(m.points, axis=1).max() == pytest.approx(1.0)


def test_vtk_lagrange_keeps_side_regions():
    linear = mp.box(dims=(2, 2, 2))
    cubic = mp.box(dims=(2, 2, 2), order=3)
    for name in ("xmin", "zmax"):
        assert np.array_equal(
            np.asarray(_region(linear, name).entries),
            np.asarray(_region(cubic, name).entries),
        )


@pytest.mark.parametrize(
    "call",
    [
        lambda: mp.ball(n=1, order=3),
        lambda: mp.disk(n=2, order=4),
        lambda: mp.box(cell_type="tetra", order=3),
        lambda: mp.helix(order=3),
        lambda: mp.box(order=2, complete=True),
        lambda: mp.cylinder(nxy=2, nz=1, cell_type="wedge", order=2, complete=True),
    ],
)
def test_high_order_round_trips_through_vtu(tmp_path, call):
    m = call()
    path = tmp_path / "m.vtu"
    mp.write(path, m)
    r = mp.read(path)
    np.testing.assert_array_equal(r.points, m.points)
    for a, b in zip(m.cells, r.cells):
        assert a.type == b.type
        assert np.array_equal(np.asarray(a.data), np.asarray(b.data))
