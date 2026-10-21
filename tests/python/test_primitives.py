"""Parametric generators: rectangle, box, disk, uv_sphere, icosa_sphere, cylinder,
the planar domains (reference_triangle, ngon, annulus, lshape, plate_with_hole,
honeycomb) and reference_cells (``meshioplusplus._primitives`` and
``operations/primitives.cpp``).
"""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core
from meshioplusplus import _primitives as prim

_KUHN = [
    [0, 1, 2, 6],
    [0, 2, 3, 6],
    [0, 3, 7, 6],
    [0, 7, 4, 6],
    [0, 4, 5, 6],
    [0, 5, 1, 6],
]


def _tet_volume(points, conn):
    a, b, c, d = (points[conn[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def _cell_volumes(mesh):
    """Signed volumes of the volume blocks (corner nodes only)."""
    out = {}
    for block in mesh.cells:
        conn = np.asarray(block.data)
        if block.type in ("tetra", "tetra10"):
            out[block.type] = _tet_volume(mesh.points, conn)
        elif block.type in ("hexahedron", "hexahedron20"):
            out[block.type] = sum(_tet_volume(mesh.points, conn[:, t]) for t in _KUHN)
        elif block.type in ("wedge", "wedge15"):
            tets = [[0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]]
            out[block.type] = sum(_tet_volume(mesh.points, conn[:, t]) for t in tets)
    return out


def _planar_area(mesh):
    out = {}
    for block in mesh.cells:
        conn = np.asarray(block.data)
        p = mesh.points
        if block.type in ("triangle", "triangle6"):
            out[block.type] = (
                np.cross(p[conn[:, 1]] - p[conn[:, 0]], p[conn[:, 2]] - p[conn[:, 0]])[
                    :, 2
                ]
                / 2
            )
        elif block.type in ("quad", "quad8"):
            diag = np.cross(
                p[conn[:, 2]] - p[conn[:, 0]], p[conn[:, 3]] - p[conn[:, 1]]
            )
            out[block.type] = diag[:, 2] / 2
    return out


def _outward(mesh):
    """Number of surface cells whose normal points toward the centre."""
    bad = 0
    p = mesh.points
    for block in mesh.cells:
        conn = np.asarray(block.data)
        if block.type.startswith("triangle"):
            a, b, c = p[conn[:, 0]], p[conn[:, 1]], p[conn[:, 2]]
            n, cen = np.cross(b - a, c - a), (a + b + c) / 3
        else:
            a, b, c, d = (p[conn[:, i]] for i in range(4))
            n, cen = np.cross(c - a, d - b), (a + b + c + d) / 4
        bad += int((np.einsum("ij,ij->i", n, cen) <= 0).sum())
    return bad


def _region(mesh, name):
    return next(r for r in mesh.regions if r.name == name)


def _same_mesh(a, b, exact):
    assert [(c.type, np.asarray(c.data).shape) for c in a.cells] == [
        (c.type, np.asarray(c.data).shape) for c in b.cells
    ]
    for ca, cb in zip(a.cells, b.cells):
        assert np.array_equal(np.asarray(ca.data), np.asarray(cb.data))
    if exact:
        assert np.asarray(a.points).tobytes() == np.asarray(b.points).tobytes()
    else:
        np.testing.assert_allclose(a.points, b.points, rtol=1e-13, atol=1e-14)
    ra = {r.name: np.asarray(r.entries) for r in a.regions}
    rb = {r.name: np.asarray(r.entries) for r in b.regions}
    assert ra.keys() == rb.keys()
    for name in ra:
        assert np.array_equal(ra[name], rb[name]), name


# --------------------------------------------------------------------------- #
# contract                                                                    #
# --------------------------------------------------------------------------- #


def test_public_names():
    for name in (
        "rectangle",
        "box",
        "disk",
        "uv_sphere",
        "icosa_sphere",
        "cylinder",
        "reference_triangle",
        "ngon",
        "annulus",
        "lshape",
        "plate_with_hole",
        "honeycomb",
        "reference_cells",
    ):
        assert name in mp.__all__
        assert callable(getattr(mp, name))


@pytest.mark.parametrize(
    "call",
    [
        lambda: mp.rectangle(dims=(0, 3)),
        lambda: mp.box(dims=(2, 0, 2)),
        lambda: mp.disk(n=0),
        lambda: mp.uv_sphere(nu=0),
        lambda: mp.cylinder(nz=0),
        lambda: mp.reference_triangle(n=0),
        lambda: mp.ngon(sides=0),
        lambda: mp.ngon(rings=0),
        lambda: mp.annulus(ntheta=0),
        lambda: mp.annulus(nr=0),
        lambda: mp.lshape(n=0),
        lambda: mp.plate_with_hole(nt=0),
        lambda: mp.plate_with_hole(nr=0),
        lambda: mp.honeycomb(dims=(0, 3)),
        lambda: mp.honeycomb(dims=(3, 0)),
    ],
)
def test_zero_count_is_an_empty_mesh(call):
    out = call()
    assert out.points.shape == (0, 3)
    assert len(out.cells) == 0


@pytest.mark.parametrize(
    "call, text",
    [
        (lambda: mp.rectangle(dims=(-1, 2)), "must not be negative"),
        (lambda: mp.rectangle(size=(0.0, 1.0)), "size must be positive"),
        (lambda: mp.rectangle(cell_type="hex"), "unknown cell_type 'hex'"),
        (lambda: mp.rectangle(split="diag"), "unknown split 'diag'"),
        (lambda: mp.rectangle(order=3), "order must be 1 or 2"),
        (lambda: mp.box(cell_type="quad"), "unknown cell_type 'quad'"),
        (lambda: mp.disk(radius=-1.0), "radius must be positive"),
        (lambda: mp.uv_sphere(nu=2), "nu must be at least 3"),
        (lambda: mp.uv_sphere(nv=1), "nv must be at least 2"),
        (lambda: mp.icosa_sphere(subdivisions=-1), "must not be negative"),
        (lambda: mp.cylinder(height=0.0), "height must be positive"),
        (lambda: mp.reference_triangle(n=-1), "must not be negative"),
        (lambda: mp.reference_triangle(size=(1.0, 0.0)), "size must be positive"),
        (lambda: mp.ngon(sides=2), "sides must be at least 3"),
        (lambda: mp.ngon(radius=-1.0), "radius must be positive"),
        (lambda: mp.annulus(inner_radius=1.0, outer_radius=0.5), "must exceed"),
        (lambda: mp.annulus(inner_radius=0.0), "inner_radius must be positive"),
        (lambda: mp.annulus(ntheta=2), "at least 3 for a full ring"),
        (lambda: mp.annulus(angle=7.0), r"angle must be in \(0, 2\*pi\]"),
        (lambda: mp.annulus(ntheta=1, angle=4.0), "must span less than pi"),
        (lambda: mp.annulus(cell_type="hexahedron"), "unknown cell_type"),
        (lambda: mp.lshape(size=0.0), "size must be positive"),
        (lambda: mp.lshape(cell_type="wedge"), "unknown cell_type 'wedge'"),
        (lambda: mp.plate_with_hole(radius=0.6, size=1.0), "less than half the size"),
        (lambda: mp.plate_with_hole(order=3), "order must be 1 or 2"),
        (lambda: mp.honeycomb(dims=(2, -1)), "must not be negative"),
        (lambda: mp.honeycomb(edge=0.0), "edge must be positive"),
        (lambda: mp.reference_cells(["blob"]), "unknown cell type 'blob'"),
        (lambda: mp.reference_cells(["quad", "quad"]), "listed twice"),
        (lambda: mp.reference_cells(["tetra10"], pairs=True), "pairs=True supports"),
    ],
)
def test_bad_arguments_raise_by_name(call, text):
    with pytest.raises(ValueError, match=r"meshio\+\+: .*" + text):
        call()


@pytest.mark.parametrize(
    "call",
    [
        lambda: mp.rectangle(dims=(100, 100), max_cells=50),
        lambda: mp.box(dims=(10, 10, 10), max_cells=50),
        lambda: mp.disk(n=20, max_cells=50),
        lambda: mp.uv_sphere(nu=64, nv=64, max_cells=50),
        lambda: mp.icosa_sphere(subdivisions=6, max_cells=50),
        lambda: mp.cylinder(nxy=10, nz=10, max_cells=50),
        lambda: mp.reference_triangle(n=20, max_cells=50),
        lambda: mp.ngon(rings=10, max_cells=50),
        lambda: mp.annulus(ntheta=20, nr=20, max_cells=50),
        lambda: mp.lshape(n=10, max_cells=50),
        lambda: mp.plate_with_hole(nt=8, nr=8, max_cells=50),
        lambda: mp.honeycomb(dims=(10, 10), max_cells=50),
    ],
)
def test_max_cells_budget(call):
    with pytest.raises(ValueError, match="above the limit of 50"):
        call()


def test_max_cells_zero_disables_budget():
    assert len(mp.rectangle(dims=(8, 8), max_cells=0).cells) == 1


# --------------------------------------------------------------------------- #
# rectangle and box                                                           #
# --------------------------------------------------------------------------- #


def test_rectangle_quads():
    m = mp.rectangle(dims=(3, 2), size=(3.0, 1.0), origin=(1.0, 2.0))
    assert m.points.shape == (12, 3)
    np.testing.assert_array_equal(m.points[:4, 0], [1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(m.points[::4, 1], [2.0, 2.5, 3.0])
    assert (m.points[:, 2] == 0.0).all()
    assert m.cells[0].type == "quad"
    np.testing.assert_array_equal(m.cells[0].data[0], [0, 1, 5, 4])
    assert (_planar_area(m)["quad"] > 0).all()
    assert abs(_planar_area(m)["quad"].sum() - 3.0) < 1e-14


@pytest.mark.parametrize(
    "split, n", [("zigzag", 24), ("up", 24), ("down", 24), ("center", 48)]
)
def test_rectangle_triangles(split, n):
    m = mp.rectangle(dims=(4, 3), cell_type="triangle", split=split)
    assert len(m.cells[0].data) == n
    area = _planar_area(m)["triangle"]
    assert (area > 0).all()
    assert abs(area.sum() - 1.0) < 1e-13


def test_rectangle_zigzag_alternates_diagonals():
    m = mp.rectangle(dims=(2, 1), cell_type="triangle", split="zigzag")
    conn = np.asarray(m.cells[0].data)
    assert [0, 2, 3] in conn.tolist() or [0, 1, 4] in conn.tolist()
    assert conn[:2].tolist() == [[0, 1, 4], [0, 4, 3]]
    assert conn[2:].tolist() == [[1, 2, 4], [2, 5, 4]]


def test_rectangle_regions():
    m = mp.rectangle(dims=(3, 2))
    assert {r.name for r in m.regions} == {"xmin", "xmax", "ymin", "ymax"}
    assert all(r.kind == "side" and r.dim == 1 for r in m.regions)
    assert [len(_region(m, n).entries) for n in ("xmin", "xmax", "ymin", "ymax")] == [
        2,
        2,
        3,
        3,
    ]
    # quad edges: 0 = ymin, 1 = xmax, 2 = ymax, 3 = xmin
    assert sorted(map(tuple, _region(m, "xmin").entries.tolist())) == [(0, 3), (3, 3)]
    assert sorted(map(tuple, _region(m, "ymax").entries.tolist())) == [
        (3, 2),
        (4, 2),
        (5, 2),
    ]


def test_box_hexahedra():
    m = mp.box(dims=(2, 3, 4), size=(2.0, 3.0, 4.0), origin=(-1.0, 0.0, 5.0))
    assert m.points.shape == (60, 3)
    assert m.cells[0].type == "hexahedron" and len(m.cells[0].data) == 24
    np.testing.assert_array_equal(m.cells[0].data[0], [0, 1, 4, 3, 12, 13, 16, 15])
    vol = _cell_volumes(m)["hexahedron"]
    assert (vol > 0).all() and abs(vol.sum() - 24.0) < 1e-12
    assert m.points.min(axis=0).tolist() == [-1.0, 0.0, 5.0]
    assert m.points.max(axis=0).tolist() == [1.0, 3.0, 9.0]


def test_box_hexahedron_regions_use_cell_faces():
    m = mp.box(dims=(2, 2, 2))
    assert [r.name for r in sorted(m.regions, key=lambda r: r.name)] == [
        "xmax",
        "xmin",
        "ymax",
        "ymin",
        "zmax",
        "zmin",
    ]
    assert all(r.dim == 2 for r in m.regions)
    xmin = _region(m, "xmin").entries
    assert len(xmin) == 4 and (xmin[:, 1] == 0).all()
    assert (_region(m, "zmax").entries[:, 1] == 5).all()
    # every boundary facet is in exactly one region: 6 faces x 4
    assert sum(len(r.entries) for r in m.regions) == 24


@pytest.mark.parametrize("cell_type, per", [("tetra", 6), ("wedge", 2)])
def test_box_split_cells(cell_type, per):
    m = mp.box(dims=(2, 2, 3), cell_type=cell_type)
    assert len(m.cells[0].data) == 12 * per
    vol = _cell_volumes(m)[cell_type]
    assert (vol > 0).all() and abs(vol.sum() - 1.0) < 1e-12


def test_box_tetra_split_conforms():
    """Kuhn tetrahedra share whole facets: every interior facet has two owners."""
    m = mp.box(dims=(2, 2, 2), cell_type="tetra")
    conn = np.asarray(m.cells[0].data)
    faces = np.sort(
        np.concatenate(
            [conn[:, f] for f in ([0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3])]
        ),
        axis=1,
    )
    _, counts = np.unique(faces, axis=0, return_counts=True)
    assert set(counts.tolist()) <= {1, 2}
    assert (
        counts == 1
    ).sum() == 6 * 4 * 2  # the boundary: 6 faces x 4 quads x 2 triangles


# --------------------------------------------------------------------------- #
# disk, cylinder                                                              #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("cell_type, n, cells", [("quad", 4, 16), ("triangle", 3, 54)])
def test_disk_counts_and_rim(cell_type, n, cells):
    m = mp.disk(radius=2.0, n=n, center=(1.0, -1.0, 0.5), cell_type=cell_type)
    assert len(m.cells[0].data) == cells
    assert (_planar_area(m)[cell_type] > 0).all()
    rim = _region(m, "rim")
    assert rim.dim == 1
    assert len(rim.entries) == (4 * n if cell_type == "quad" else 6 * n)
    conn = np.asarray(m.cells[0].data)
    edges = {
        "quad": [(0, 1), (1, 2), (2, 3), (3, 0)],
        "triangle": [(0, 1), (1, 2), (2, 0)],
    }[cell_type]
    for c, f in rim.entries.tolist():
        for node in (conn[c, edges[f][0]], conn[c, edges[f][1]]):
            r = np.linalg.norm(m.points[node, :2] - [1.0, -1.0])
            assert np.isclose(r, 2.0)
    assert (m.points[:, 2] == 0.5).all()


def test_disk_quad_rim_nodes_lie_on_the_circle():
    m = mp.disk(radius=1.5, n=6)
    r = np.linalg.norm(m.points[:, :2], axis=1)
    on_rim = np.isclose(r, 1.5)
    assert on_rim.sum() == 4 * 6  # the boundary of the 7 x 7 lattice
    assert r.max() <= 1.5 + 1e-12


def test_disk_triangle_rings():
    m = mp.disk(n=3, cell_type="triangle")
    assert m.points.shape == (1 + 3 * 3 * 4, 3)
    r = np.round(np.linalg.norm(m.points[:, :2], axis=1), 12)
    assert sorted(set(r.tolist())) == [0.0, round(1 / 3, 12), round(2 / 3, 12), 1.0]
    area = _planar_area(m)["triangle"]
    assert (area > 0).all()
    # a hexagonal-ring triangulation of the unit circle: between the inscribed
    # hexagon and the disk
    assert 3 * np.sqrt(3) / 2 < area.sum() < np.pi


def test_cylinder_hexahedra():
    m = mp.cylinder(radius=1.0, height=4.0, nxy=4, nz=3, center=(0.0, 0.0, 10.0))
    assert len(m.cells[0].data) == 4 * 4 * 3
    vol = _cell_volumes(m)["hexahedron"]
    assert (vol > 0).all()
    assert m.points[:, 2].min() == 8.0 and m.points[:, 2].max() == 12.0
    names = {r.name: len(r.entries) for r in m.regions}
    assert names == {"side": 4 * 4 * 3, "top": 16, "bottom": 16}


def test_cylinder_wedges():
    m = mp.cylinder(nxy=2, nz=2, cell_type="wedge")
    assert m.cells[0].type == "wedge" and len(m.cells[0].data) == 6 * 4 * 2
    assert (_cell_volumes(m)["wedge"] > 0).all()
    names = {r.name: len(r.entries) for r in m.regions}
    assert names == {"side": 12 * 2, "top": 24, "bottom": 24}


def test_cylinder_single_cell_keeps_caps_apart_from_the_side():
    m = mp.cylinder(nxy=1, nz=1)
    assert {r.name: len(r.entries) for r in m.regions} == {
        "side": 4,
        "top": 1,
        "bottom": 1,
    }


# --------------------------------------------------------------------------- #
# spheres                                                                     #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "kwargs, cells",
    [
        (dict(cell_type="triangle"), {"triangle": 2 * 12 * 5}),
        (dict(cell_type="quad"), {"quad": 12 * 4, "triangle": 24}),
    ],
)
def test_uv_sphere(kwargs, cells):
    m = mp.uv_sphere(radius=2.0, nu=12, nv=6, center=(1.0, 2.0, 3.0), **kwargs)
    assert {c.type: len(c.data) for c in m.cells} == cells
    assert m.points.shape == (2 + 12 * 5, 3)
    np.testing.assert_allclose(m.points[0], [1.0, 2.0, 5.0])
    np.testing.assert_allclose(m.points[-1], [1.0, 2.0, 1.0])
    assert np.allclose(np.linalg.norm(m.points - [1.0, 2.0, 3.0], axis=1), 2.0)
    centred = mp.Mesh(m.points - [1.0, 2.0, 3.0], m.cells)
    assert _outward(centred) == 0


def test_uv_sphere_poles_are_exact_and_welded():
    m = mp.uv_sphere(nu=5, nv=4)
    assert m.points[0].tolist() == [0.0, 0.0, 1.0]
    assert m.points[-1].tolist() == [0.0, 0.0, -1.0]
    conn = np.asarray(m.cells[0].data)
    assert np.unique(m.points, axis=0).shape[0] == m.points.shape[0]
    assert conn.max() == m.points.shape[0] - 1


def test_icosa_sphere_counts_and_outward():
    for s in range(4):
        m = mp.icosa_sphere(radius=3.0, subdivisions=s)
        assert len(m.cells[0].data) == 20 * 4**s
        assert m.points.shape[0] == 10 * 4**s + 2
        assert np.allclose(np.linalg.norm(m.points, axis=1), 3.0)
        assert _outward(m) == 0


def test_icosa_sphere_is_closed():
    m = mp.icosa_sphere(subdivisions=2)
    conn = np.asarray(m.cells[0].data)
    edges = np.sort(
        np.concatenate([conn[:, [0, 1]], conn[:, [1, 2]], conn[:, [2, 0]]]), axis=1
    )
    _, counts = np.unique(edges, axis=0, return_counts=True)
    assert (counts == 2).all()
    # outward orientation: the enclosed volume is positive and close to 4/3 pi
    p = m.points
    vol = (
        np.einsum(
            "ij,ij->i", p[conn[:, 0]], np.cross(p[conn[:, 1]], p[conn[:, 2]])
        ).sum()
        / 6
    )
    assert 0.9 * 4 / 3 * np.pi < vol < 4 / 3 * np.pi


# --------------------------------------------------------------------------- #
# order=2                                                                     #
# --------------------------------------------------------------------------- #

_QUADRATIC = {
    "triangle": "triangle6",
    "quad": "quad8",
    "tetra": "tetra10",
    "hexahedron": "hexahedron20",
    "wedge": "wedge15",
}


@pytest.mark.parametrize(
    "call, base",
    [
        (lambda: mp.rectangle(dims=(3, 2), order=2), "quad"),
        (lambda: mp.rectangle(dims=(3, 2), cell_type="triangle", order=2), "triangle"),
        (lambda: mp.box(dims=(2, 2, 2), order=2), "hexahedron"),
        (lambda: mp.box(dims=(2, 2, 2), cell_type="tetra", order=2), "tetra"),
        (lambda: mp.box(dims=(2, 2, 2), cell_type="wedge", order=2), "wedge"),
        (lambda: mp.reference_triangle(n=3, order=2), "triangle"),
        (lambda: mp.ngon(sides=5, rings=2, order=2), "triangle"),
        (lambda: mp.lshape(n=2, order=2), "quad"),
        (lambda: mp.lshape(n=2, cell_type="triangle", order=2), "triangle"),
    ],
)
def test_straight_shapes_get_exact_midpoints(call, base):
    m = call()
    block = m.cells[0]
    assert block.type == _QUADRATIC[base]
    from meshioplusplus._convert_cells import _ELEVATE

    _, corners, edges = _ELEVATE[base]
    conn = np.asarray(block.data)
    assert conn.shape[1] == corners + len(edges)
    for k, (a, b) in enumerate(edges):
        mid = (m.points[conn[:, a]] + m.points[conn[:, b]]) * 0.5
        np.testing.assert_allclose(m.points[conn[:, corners + k]], mid, atol=1e-15)
    # mid nodes are shared between neighbours: no duplicate positions
    assert np.unique(m.points.round(12), axis=0).shape[0] == m.points.shape[0]
    # first-appearance order: the new ids rise monotonically through the cells
    new = conn[:, corners:].reshape(-1)
    firsts = new[np.sort(np.unique(new, return_index=True)[1])]
    assert (np.diff(firsts) > 0).all()


def test_disk_and_cylinder_rims_are_curved():
    for m in (
        mp.disk(radius=2.0, n=3, order=2),
        mp.disk(radius=2.0, n=2, cell_type="triangle", order=2),
    ):
        conn = np.asarray(m.cells[0].data)
        rim = _region(m, "rim").entries
        for c, f in rim.tolist():
            row = conn[c]
            if m.cells[0].type == "quad8":
                nodes = row[[f, (f + 1) % 4, 4 + f]]
            else:
                nodes = row[[f, (f + 1) % 3, 3 + f]]
            assert np.allclose(np.linalg.norm(m.points[nodes][:, :2], axis=1), 2.0)
    m = mp.cylinder(radius=2.0, nxy=2, nz=2, cell_type="wedge", order=2)
    on_side = np.isclose(np.linalg.norm(m.points[:, :2], axis=1), 2.0)
    assert on_side.sum() > 0
    assert (np.linalg.norm(m.points[:, :2], axis=1) <= 2.0 + 1e-12).all()


@pytest.mark.parametrize(
    "call",
    [
        lambda: mp.uv_sphere(radius=3.0, nu=7, nv=5, order=2),
        lambda: mp.uv_sphere(radius=3.0, nu=7, nv=5, cell_type="quad", order=2),
        lambda: mp.icosa_sphere(radius=3.0, subdivisions=1, order=2),
    ],
)
def test_sphere_midpoints_lie_on_the_sphere(call):
    """A naive average puts a pole mid-node off its meridian, and a seam mid-node
    on the far side of the sphere; both would leave the radius here."""
    m = call()
    assert {c.type for c in m.cells} <= {"triangle6", "quad8"}
    assert np.allclose(np.linalg.norm(m.points, axis=1), 3.0)
    assert np.unique(m.points.round(10), axis=0).shape[0] == m.points.shape[0]


def test_uv_sphere_pole_midnode_stays_on_its_meridian():
    m = mp.uv_sphere(nu=6, nv=4, order=2)
    conn = np.asarray(m.cells[0].data)
    cap = conn[:6]  # north cap: pole, ring i, ring i + 1, then three mid-nodes
    for row in cap:
        # edge 2 of triangle6 is (2, 0): ring i + 1 -> pole; node 5 is its midpoint
        mid, ring = m.points[row[5]], m.points[row[2]]
        assert np.isclose(np.arctan2(mid[1], mid[0]), np.arctan2(ring[1], ring[0]))


# --------------------------------------------------------------------------- #
# planar domains                                                              #
# --------------------------------------------------------------------------- #


def _euler(mesh):
    """V - E + F of the corner graph (corner nodes only)."""
    corners = {"triangle": 3, "triangle6": 3, "quad": 4, "quad8": 4}
    edges = set()
    faces = 0
    for block in mesh.cells:
        conn = np.asarray(block.data)
        k = corners.get(block.type, conn.shape[1])
        faces += len(conn)
        for row in conn:
            for i in range(k):
                a, b = int(row[i]), int(row[(i + 1) % k])
                edges.add((min(a, b), max(a, b)))
    return len(mesh.points) - len(edges) + faces


def _total_area(mesh):
    return sum(float(a.sum()) for a in _planar_area(mesh).values())


@pytest.mark.parametrize("n", [1, 2, 5])
def test_reference_triangle_counts_area_and_regions(n):
    m = mp.reference_triangle(n=n, size=(2.0, 3.0), origin=(1.0, -1.0, 4.0))
    assert m.points.shape == ((n + 1) * (n + 2) // 2, 3)
    assert len(m.cells[0].data) == n * n
    assert (_planar_area(m)["triangle"] > 0).all()
    assert abs(_total_area(m) - 3.0) < 1e-12
    assert {r.name for r in m.regions} == {"ymin", "xmin", "hypotenuse"}
    assert all(len(r.entries) == n and r.dim == 1 for r in m.regions)
    assert _euler(m) == 1
    assert m.points[0].tolist() == [1.0, -1.0, 4.0]


def test_reference_triangle_barycentric_order():
    m = mp.reference_triangle(n=2)
    np.testing.assert_array_equal(
        m.points[:, :2],
        [[0, 0], [0.5, 0], [1, 0], [0, 0.5], [0.5, 0.5], [0, 1]],
    )
    # row 0 holds two lattice triangles and the inverted one between them
    assert m.cells[0].data[0].tolist() == [0, 1, 3]
    assert m.cells[0].data[1].tolist() == [1, 4, 3]


@pytest.mark.parametrize("sides", [3, 5, 8])
@pytest.mark.parametrize("rings", [1, 3])
def test_ngon_fans_and_rings(sides, rings):
    m = mp.ngon(sides=sides, radius=2.0, rings=rings, center=(1.0, 1.0, 0.0))
    assert len(m.cells[0].data) == sides * rings * rings
    assert len(m.points) == 1 + sides * rings * (rings + 1) // 2
    assert (_planar_area(m)["triangle"] > 0).all()
    want = 0.5 * sides * 4.0 * np.sin(2 * np.pi / sides)
    assert abs(_total_area(m) - want) < 1e-12
    assert len(_region(m, "rim").entries) == sides * rings
    assert np.unique(m.points.round(10), axis=0).shape[0] == len(m.points)
    assert _euler(m) == 1


@pytest.mark.parametrize("cell_type", ["quad", "triangle"])
def test_annulus_full_ring_is_welded(cell_type):
    m = mp.annulus(0.5, 1.5, 12, 3, center=(1.0, 2.0, 3.0), cell_type=cell_type)
    assert len(m.points) == 12 * 4
    assert len(m.cells[0].data) == (36 if cell_type == "quad" else 72)
    assert np.unique(m.points.round(10), axis=0).shape[0] == len(m.points)
    assert all(a.min() > 0 for a in _planar_area(m).values())
    assert _euler(m) == 0
    assert {r.name for r in m.regions} == {"inner", "outer"}
    assert all(len(r.entries) == 12 for r in m.regions)
    r = np.linalg.norm(m.points[:, :2] - [1.0, 2.0], axis=1)
    assert r.min() > 0.5 - 1e-12 and r.max() < 1.5 + 1e-12
    assert (m.points[:, 2] == 3.0).all()
    rings = 0.5 + np.arange(4) / 3
    area = sum(
        6 * (b * b - a * a) * np.sin(np.pi / 6) for a, b in zip(rings, rings[1:])
    )
    assert abs(_total_area(m) - area) < 1e-12


def test_annulus_sector_keeps_its_seam_apart():
    m = mp.annulus(1.0, 2.0, 6, 2, angle=np.pi / 2)
    assert len(m.points) == 7 * 3 and len(m.cells[0].data) == 12
    assert {r.name for r in m.regions} == {"inner", "outer", "start", "end"}
    assert [len(_region(m, n).entries) for n in ("start", "end")] == [2, 2]
    assert _euler(m) == 1
    np.testing.assert_allclose(m.points[0], [1.0, 0.0, 0.0], atol=1e-15)
    np.testing.assert_allclose(m.points[6], [0.0, 1.0, 0.0], atol=1e-15)


def test_annulus_sector_seams_link_with_periodic():
    m = mp.annulus(1.0, 2.0, 4, 3, angle=np.pi / 2)
    start = {tuple(m.points[i].round(12)) for i in np.unique(_edge_nodes(m, "start"))}
    end = {tuple(m.points[i].round(12)) for i in np.unique(_edge_nodes(m, "end"))}
    # a quarter turn about the origin maps the start ray onto the end ray
    turned = {(-y, x, z) for x, y, z in start}
    assert {tuple(np.round(p, 12)) for p in turned} == end


def _edge_nodes(mesh, name):
    conn = np.asarray(mesh.cells[0].data)
    nodes = []
    for c, f in _region(mesh, name).entries.tolist():
        nodes += [conn[c, f], conn[c, (f + 1) % conn.shape[1]]]
    return np.array(nodes)


@pytest.mark.parametrize("cell_type", ["quad", "triangle"])
def test_lshape_loses_a_quadrant(cell_type):
    m = mp.lshape(n=3, size=2.0, center=(0.0, 0.0, 1.0), cell_type=cell_type)
    assert len(m.cells[0].data) == (27 if cell_type == "quad" else 54)
    assert len(m.points) == 7 * 7 - 9
    assert all(a.min() > 0 for a in _planar_area(m).values())
    assert abs(_total_area(m) - 3.0) < 1e-12
    counts = {r.name: len(r.entries) for r in m.regions}
    assert counts == {
        "xmin": 6,
        "xmax": 3,
        "ymin": 3,
        "ymax": 6,
        "notch_x": 3,
        "notch_y": 3,
    }
    assert _euler(m) == 1
    inside = (m.points[:, 0] > 1e-12) & (m.points[:, 1] < -1e-12)
    assert not inside.any()
    assert ((m.points[:, 0] == 0.0) & (m.points[:, 1] == 0.0)).sum() == 1


@pytest.mark.parametrize("cell_type", ["quad", "triangle"])
def test_plate_with_hole_surrounds_the_hole(cell_type):
    m = mp.plate_with_hole(
        radius=0.25, size=2.0, nt=3, nr=2, center=(1.0, 1.0, 0.0), cell_type=cell_type
    )
    assert len(m.cells[0].data) == (24 if cell_type == "quad" else 48)
    assert len(m.points) == 3 * 12
    assert np.unique(m.points.round(10), axis=0).shape[0] == len(m.points)
    assert all(a.min() > 0 for a in _planar_area(m).values())
    assert _euler(m) == 0
    counts = {r.name: len(r.entries) for r in m.regions}
    assert counts == {"hole": 12, "xmin": 3, "xmax": 3, "ymin": 3, "ymax": 3}
    hole = 0.5 * 12 * 0.0625 * np.sin(2 * np.pi / 12)
    assert abs(_total_area(m) - (4.0 - hole)) < 1e-12
    assert np.abs(m.points[:, :2] - 1.0).max() <= 1.0 + 1e-12


def test_plate_with_hole_quadratic_hole_is_circular():
    m = mp.plate_with_hole(radius=0.5, size=2.0, nt=2, nr=3, order=2)
    assert m.cells[0].type == "quad8"
    conn = np.asarray(m.cells[0].data)
    for c, f in _region(m, "hole").entries.tolist():
        for node in (conn[c, f], conn[c, 4 + f]):
            assert np.isclose(np.hypot(*m.points[node, :2]), 0.5)
    for c, f in _region(m, "xmax").entries.tolist():
        assert m.points[conn[c, 4 + f], 0] == 1.0


@pytest.mark.parametrize("dims", [(1, 1), (3, 2), (4, 5)])
def test_honeycomb_is_a_hexagon_polygon_block(dims):
    m = mp.honeycomb(dims=dims, edge=2.0, origin=(1.0, 2.0, 0.0))
    cells = dims[0] * dims[1]
    assert len(m.cells) == 1 and m.cells[0].type.startswith("polygon")
    conn = np.asarray(m.cells[0].data)
    assert conn.shape == (cells, 6)
    p = m.points[:, :2]
    area = 0.0
    for row in conn:
        ring = p[row]
        nxt = np.roll(ring, -1, axis=0)
        twice = float((ring[:, 0] * nxt[:, 1] - nxt[:, 0] * ring[:, 1]).sum())
        assert twice > 0
        assert np.allclose(np.linalg.norm(nxt - ring, axis=1), 2.0)
        area += twice / 2
    assert np.isclose(area, cells * 1.5 * np.sqrt(3) * 4.0)
    assert np.unique(p.round(10), axis=0).shape[0] == len(p)
    assert p.min(axis=0).tolist() == [1.0, 2.0]
    assert _euler(m) == 1
    boundary = _region(m, "boundary")
    assert boundary.dim == 1
    edges = {
        (min(int(r[i]), int(r[(i + 1) % 6])), max(int(r[i]), int(r[(i + 1) % 6])))
        for r in conn
        for i in range(6)
    }
    assert len(boundary.entries) == 2 * len(edges) - 6 * cells


def test_honeycomb_numbers_nodes_by_first_appearance():
    m = mp.honeycomb(dims=(2, 2))
    assert np.asarray(m.cells[0].data)[0].tolist() == [0, 1, 2, 3, 4, 5]
    assert len(m.points) == 16


def test_planar_shapes_run_through_operations():
    """The domains feed refine, quality, the writers and the ragged paths."""
    for m in (
        mp.reference_triangle(n=3),
        mp.ngon(sides=6, rings=2),
        mp.annulus(),
        mp.lshape(n=2),
        mp.plate_with_hole(nt=4, nr=2),
    ):
        mp.compute_quality(m)
    out = mp.refine(mp.lshape(n=2))
    assert len(out.cells[0].data) > 12


# --------------------------------------------------------------------------- #
# reference cells                                                             #
# --------------------------------------------------------------------------- #


def test_reference_cells_default_covers_every_type():
    m = mp.reference_cells()
    assert [c.type for c in m.cells] == [
        "vertex",
        "line",
        "triangle",
        "quad",
        "tetra",
        "hexahedron",
        "wedge",
        "pyramid",
        "line3",
        "triangle6",
        "quad8",
        "tetra10",
        "hexahedron20",
        "wedge15",
        "pyramid13",
        "polyhedron8",
    ]
    assert all(len(c.data) == 1 for c in m.cells)


def test_reference_cells_are_positively_oriented_and_apart():
    m = mp.reference_cells(["triangle", "quad", "tetra", "hexahedron", "wedge"])
    for kind in ("triangle", "quad"):
        assert (_planar_area(m)[kind] > 0).all()
    vols = _cell_volumes(m)
    assert np.isclose(vols["tetra"][0], 1 / 6)
    assert np.isclose(vols["hexahedron"][0], 1.0)
    assert np.isclose(vols["wedge"][0], 0.5)
    spans = [
        (
            m.points[np.asarray(c.data).reshape(-1)][:, 0].min(),
            m.points[np.asarray(c.data).reshape(-1)][:, 0].max(),
        )
        for c in m.cells
    ]
    assert all(spans[i][1] < spans[i + 1][0] for i in range(len(spans) - 1))


def test_reference_cells_quadratic_midpoints_are_exact():
    m = mp.reference_cells(["tetra10", "hexahedron20", "wedge15", "pyramid13"])
    from meshioplusplus._convert_cells import _ELEVATE

    for block, base in zip(m.cells, ("tetra", "hexahedron", "wedge", "pyramid")):
        row = np.asarray(block.data)[0]
        _, corners, edges = _ELEVATE[base]
        for k, (a, b) in enumerate(edges):
            np.testing.assert_allclose(
                m.points[row[corners + k]], (m.points[row[a]] + m.points[row[b]]) / 2
            )


def test_reference_cells_polyhedron_is_a_closed_cube():
    m = mp.reference_cells(["polyhedron"])
    assert m.cells[0].type == "polyhedron8"
    faces = m.cells[0].data[0]
    assert len(faces) == 6
    edges = np.sort(
        np.array([[f[i], f[(i + 1) % len(f)]] for f in faces for i in range(len(f))]),
        axis=1,
    )
    _, counts = np.unique(edges, axis=0, return_counts=True)
    assert (counts == 2).all()


@pytest.mark.parametrize("name", prim._PAIR_TYPES)
def test_reference_cell_pairs_share_a_facet(name):
    m = mp.reference_cells([name], pairs=True)
    conn = np.asarray(m.cells[0].data)
    assert conn.shape[0] == 2
    shared = set(conn[0].tolist()) & set(conn[1].tolist())
    expect = {
        "line": 1,
        "triangle": 2,
        "quad": 2,
        "tetra": 3,
        "hexahedron": 4,
        "wedge": 3,
        "pyramid": 4,
    }
    assert len(shared) == expect[name]
    for vol in _cell_volumes(m).values():
        assert (vol > 0).all()


# --------------------------------------------------------------------------- #
# twin vs core                                                                #
# --------------------------------------------------------------------------- #

_CASES = [
    ("rectangle", dict(dims=(5, 3), size=(2.0, 0.5), origin=(0.1, 0.2, 0.3)), True),
    ("rectangle", dict(dims=(4, 3), cell_type="triangle", split="zigzag"), True),
    ("rectangle", dict(dims=(4, 3), cell_type="triangle", split="up", order=2), True),
    ("rectangle", dict(dims=(4, 3), cell_type="triangle", split="down"), True),
    (
        "rectangle",
        dict(dims=(4, 3), cell_type="triangle", split="center", order=2),
        True,
    ),
    ("rectangle", dict(dims=(3, 3), order=2), True),
    ("box", dict(dims=(3, 2, 4), size=(1.5, 2.5, 3.5), origin=(0.1, 0.2, 0.3)), True),
    ("box", dict(dims=(3, 2, 2), cell_type="tetra", order=2), True),
    ("box", dict(dims=(3, 2, 2), cell_type="wedge", order=2), True),
    ("box", dict(dims=(2, 2, 2), order=2), True),
    ("icosa_sphere", dict(radius=2.5, subdivisions=3, center=(1.0, 2.0, 3.0)), True),
    ("icosa_sphere", dict(subdivisions=2, order=2), True),
    ("reference_cells", dict(), True),
    ("reference_cells", dict(pairs=True), True),
    ("reference_cells", dict(types=["hexahedron20", "polyhedron", "vertex"]), True),
    ("disk", dict(radius=1.7, n=5, center=(0.5, 0.25, 1.0)), False),
    ("disk", dict(n=4, order=2), False),
    ("disk", dict(n=4, cell_type="triangle"), False),
    ("disk", dict(n=3, cell_type="triangle", order=2), False),
    ("uv_sphere", dict(radius=1.3, nu=9, nv=5), False),
    ("uv_sphere", dict(nu=7, nv=4, cell_type="quad", order=2), False),
    ("uv_sphere", dict(nu=6, nv=2), False),
    (
        "cylinder",
        dict(radius=1.2, height=3.0, nxy=3, nz=4, center=(1.0, 0.0, -1.0)),
        False,
    ),
    ("cylinder", dict(nxy=2, nz=2, order=2), False),
    ("cylinder", dict(nxy=2, nz=3, cell_type="wedge", order=2), False),
    (
        "reference_triangle",
        dict(n=5, size=(2.0, 0.5), origin=(0.1, 0.2, 0.3)),
        True,
    ),
    ("reference_triangle", dict(n=3, order=2), True),
    ("ngon", dict(sides=7, radius=1.4, rings=3, center=(0.5, 0.25, 1.0)), False),
    ("ngon", dict(sides=5, rings=2, order=2), False),
    ("annulus", dict(inner_radius=0.4, outer_radius=1.7, ntheta=9, nr=3), False),
    ("annulus", dict(ntheta=8, nr=2, cell_type="triangle", order=2), False),
    ("annulus", dict(ntheta=5, nr=2, angle=2.5, order=2), False),
    ("lshape", dict(n=3, size=1.3, center=(0.1, 0.2, 0.3)), True),
    ("lshape", dict(n=2, cell_type="triangle", order=2), True),
    ("plate_with_hole", dict(radius=0.3, size=1.7, nt=3, nr=2), False),
    ("plate_with_hole", dict(nt=2, nr=2, cell_type="triangle", order=2), False),
    ("honeycomb", dict(dims=(3, 2), edge=1.3, origin=(0.1, 0.2, 0.3)), True),
]


def _twin(name, kwargs):
    """Run the numpy twin directly, bypassing the core."""
    defaults = {
        "rectangle": dict(
            dims=(1, 1),
            size=(1.0, 1.0),
            origin=(0.0, 0.0, 0.0),
            cell_type="quad",
            split="zigzag",
            order=1,
        ),
        "box": dict(
            dims=(1, 1, 1),
            size=(1.0, 1.0, 1.0),
            origin=(0.0, 0.0, 0.0),
            cell_type="hexahedron",
            order=1,
        ),
        "disk": dict(
            radius=1.0, n=8, center=(0.0, 0.0, 0.0), cell_type="quad", order=1
        ),
        "uv_sphere": dict(
            radius=1.0,
            nu=16,
            nv=8,
            center=(0.0, 0.0, 0.0),
            cell_type="triangle",
            order=1,
        ),
        "icosa_sphere": dict(
            radius=1.0, subdivisions=2, center=(0.0, 0.0, 0.0), order=1
        ),
        "cylinder": dict(
            radius=1.0,
            height=1.0,
            nxy=8,
            nz=8,
            center=(0.0, 0.0, 0.0),
            cell_type="hexahedron",
            order=1,
        ),
        "reference_triangle": dict(
            n=4, size=(1.0, 1.0), origin=(0.0, 0.0, 0.0), order=1
        ),
        "ngon": dict(sides=6, radius=1.0, rings=1, center=(0.0, 0.0, 0.0), order=1),
        "annulus": dict(
            inner=0.5,
            outer=1.0,
            ntheta=16,
            nr=4,
            angle=2.0 * np.pi,
            center=(0.0, 0.0, 0.0),
            cell_type="quad",
            order=1,
        ),
        "lshape": dict(
            n=4, size=2.0, center=(0.0, 0.0, 0.0), cell_type="quad", order=1
        ),
        "plate_with_hole": dict(
            radius=0.25,
            size=1.0,
            nt=8,
            nr=4,
            center=(0.0, 0.0, 0.0),
            cell_type="quad",
            order=1,
        ),
        "honeycomb": dict(dims=(4, 4), edge=1.0, origin=(0.0, 0.0, 0.0)),
    }
    renames = {"annulus": {"inner_radius": "inner", "outer_radius": "outer"}}
    kwargs = {renames.get(name, {}).get(k, k): v for k, v in kwargs.items()}
    if name == "reference_cells":
        types = kwargs.get("types") or list(
            prim._PAIR_TYPES if kwargs.get("pairs") else prim._ALL_TYPES
        )
        return prim._reference_cells_py(types, kwargs.get("pairs", False))
    a = {**defaults[name], **kwargs}
    for key in ("dims", "size", "origin", "center"):
        if key in a and not np.isscalar(a[key]):
            a[key] = np.asarray(a[key], dtype=np.int64 if key == "dims" else np.float64)
    fn = {
        "rectangle": prim._rectangle_py,
        "box": prim._box_py,
        "disk": prim._disk_py,
        "uv_sphere": prim._uv_sphere_py,
        "icosa_sphere": prim._icosa_sphere_py,
        "cylinder": prim._cylinder_py,
        "reference_triangle": prim._reference_triangle_py,
        "ngon": prim._ngon_py,
        "annulus": prim._annulus_py,
        "lshape": prim._lshape_py,
        "plate_with_hole": prim._plate_py,
        "honeycomb": prim._honeycomb_py,
    }[name]
    return fn(**a)


@pytest.mark.parametrize(
    "name, kwargs, exact", _CASES, ids=lambda v: v if isinstance(v, str) else None
)
def test_twin_matches_core(name, kwargs, exact):
    core = getattr(mp, name)(**kwargs)
    twin = _twin(name, kwargs)
    _same_mesh(core, twin, exact)


def test_core_and_dispatch_agree():
    a = mp.box(dims=(2, 2, 2), order=2)
    b = _core.box(
        [2, 2, 2], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0], "hexahedron", 2, 20000000
    )
    assert np.asarray(a.points).tobytes() == np.asarray(b.points).tobytes()


def test_strict_core_never_falls_back(monkeypatch):
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    for name, kwargs, _ in _CASES:
        getattr(mp, name)(**kwargs)
