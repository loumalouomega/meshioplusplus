"""Conforming octree transitions, snapping, refinement criteria and planes.

These need the compiled core: the numpy reference implements the hanging-node
octree only, and says so by name.
"""

import itertools

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import octree, voxelize

pytest.importorskip("meshioplusplus._core")

from .test_octree import BOX, _circle, _sphere  # noqa: E402

_HEX_FACES = (
    (0, 3, 2, 1),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)


def _cell_faces(cb):
    """Every face of every cell of a block, as a tuple of node ids."""
    t = cb.type
    out = []
    if t.startswith("polyhedron"):
        for cell in cb.data:
            out.append([tuple(int(i) for i in f) for f in cell])
    elif t == "hexahedron":
        for row in np.asarray(cb.data):
            out.append([tuple(int(row[i]) for i in f) for f in _HEX_FACES])
    elif t == "pyramid":
        for r in np.asarray(cb.data):
            r = [int(i) for i in r]
            out.append(
                [(r[0], r[1], r[2], r[3])]
                + [(r[i], r[(i + 1) % 4], r[4]) for i in range(4)]
            )
    elif t == "tetra":
        for r in np.asarray(cb.data):
            r = [int(i) for i in r]
            out.append(
                [
                    (r[0], r[1], r[2]),
                    (r[0], r[1], r[3]),
                    (r[0], r[2], r[3]),
                    (r[1], r[2], r[3]),
                ]
            )
    elif t in ("quad", "triangle"):
        for r in np.asarray(cb.data):
            r = [int(i) for i in r]
            out.append([(r[i], r[(i + 1) % len(r)]) for i in range(len(r))])
    elif t.startswith("polygon"):
        for r in cb.data:
            r = [int(i) for i in r]
            out.append([(r[i], r[(i + 1) % len(r)]) for i in range(len(r))])
    else:
        raise AssertionError(f"unexpected block {t}")
    return out


def _cell_volume(cb, pts):
    """Volumes (areas in 2-D) of every cell, signed so that a valid cell is positive."""
    t = cb.type
    vols = []

    def tet(a, b, c, d):
        return (
            float(np.dot(pts[b] - pts[a], np.cross(pts[c] - pts[a], pts[d] - pts[a])))
            / 6
        )

    if t == "hexahedron":
        for r in np.asarray(cb.data):
            vols.append(
                float(
                    np.dot(
                        pts[r[1]] - pts[r[0]],
                        np.cross(pts[r[3]] - pts[r[0]], pts[r[4]] - pts[r[0]]),
                    )
                )
            )
    elif t == "pyramid":
        for r in np.asarray(cb.data):
            vols.append(tet(r[0], r[1], r[2], r[4]) + tet(r[0], r[2], r[3], r[4]))
    elif t == "tetra":
        for r in np.asarray(cb.data):
            vols.append(tet(*r))
    elif t.startswith("polyhedron"):
        for cell in cb.data:
            nodes = sorted({int(i) for f in cell for i in f})
            c = pts[nodes].mean(axis=0)
            v = 0.0
            for f in cell:
                f = [int(i) for i in f]
                for i in range(1, len(f) - 1):
                    v += (
                        tet(f[0], f[i], f[i + 1], -1)
                        if False
                        else float(
                            np.dot(
                                pts[f[0]] - c,
                                np.cross(pts[f[i]] - c, pts[f[i + 1]] - c),
                            )
                        )
                        / 6
                    )
            vols.append(v)
    else:  # 2-D: shoelace in xy
        for r in cb.data:
            xy = pts[[int(i) for i in r]][:, :2]
            x, y = xy[:, 0], xy[:, 1]
            vols.append(0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)))
    return np.array(vols)


def _check_conforming(mesh, domain_lo, domain_hi, expect_volume):
    pts = np.asarray(mesh.points)
    counts = {}
    for cb in mesh.cells:
        for faces in _cell_faces(cb):
            for f in faces:
                counts[tuple(sorted(f))] = counts.get(tuple(sorted(f)), 0) + 1
    assert counts and max(counts.values()) <= 2
    dim = 3 if any(len(k) >= 3 for k in counts) else 2
    for key, n in counts.items():
        if n == 1:  # a boundary face must lie on the domain boundary
            p = pts[list(key)][:, :dim]
            on = [
                np.allclose(p[:, k], domain_lo[k]) or np.allclose(p[:, k], domain_hi[k])
                for k in range(dim)
            ]
            assert any(on), f"face {key} is a dangling interior face"
    total = 0.0
    for cb in mesh.cells:
        v = _cell_volume(cb, pts)
        assert np.all(v > 0.0), f"{cb.type}: non-positive volume {v.min()}"
        total += float(v.sum())
    assert total == pytest.approx(expect_volume, rel=1e-9)


def _tiny_triangle(c, h):
    e = 1e-3 * h
    p = np.array(c, float)
    return [p, p + [e, 0, 0], p + [0, e, 0]]


def _forest_skin(h, offsets):
    """Tiny triangles in the neighbouring roots of a 3x3x3 forest, to force
    refinement towards the centre root from the chosen directions."""
    pts, tris = [], []
    for d in offsets:
        c = []
        for k in range(3):
            if d[k] == -1:
                c.append(1.0 * h - 3e-3 * h)
            elif d[k] == 1:
                c.append(2.0 * h + 2e-3 * h)
            else:
                c.append(1.5 * h)
        base = len(pts)
        pts.extend(_tiny_triangle(c, h))
        tris.append([base, base + 1, base + 2])
    return meshioplusplus.Mesh(
        np.array(pts), [("triangle", np.array(tris, dtype=np.int64))]
    )


_NEIGHBOURS = [d for d in itertools.product((-1, 0, 1), repeat=3) if d != (0, 0, 0)]


def _patterns():
    rng = np.random.default_rng(7)
    out = [(d,) for d in [(1, 0, 0), (0, -1, 0), (0, 0, 1)]]  # a face neighbour
    out += [((1, 1, 0),), ((-1, -1, 1),), ((1, 1, 1),)]  # an edge / a corner neighbour
    for _ in range(24):
        n = int(rng.integers(1, 7))
        idx = rng.choice(len(_NEIGHBOURS), size=n, replace=False)
        out.append(tuple(_NEIGHBOURS[i] for i in idx))
    return out


@pytest.mark.parametrize("transition", ["polyhedra", "templates"])
@pytest.mark.parametrize("balance", ["face", "full"])
@pytest.mark.parametrize("pattern", _patterns())
def test_forest_patterns_are_conforming(transition, balance, pattern):
    h = 1.0
    skin = _forest_skin(h, pattern)
    kw = dict(
        root_cells=(3, 3, 3), max_level=3, bounds=[0, 0, 0, 3, 3, 3], balance=balance
    )
    plain, rep = octree(skin, return_report=True, **kw)
    if rep["num_hanging"] == 0:
        pytest.skip("pattern produced no transition")
    mesh = octree(skin, transition=transition, **kw)
    _check_conforming(mesh, [0, 0, 0], [3, 3, 3], 27.0)
    assert "octree:hanging" not in mesh.point_data
    if transition == "polyhedra":
        assert any(cb.type.startswith("polyhedron") for cb in mesh.cells)
    else:
        assert any(cb.type in ("pyramid", "tetra") for cb in mesh.cells)
    leaf = np.concatenate([np.asarray(a) for a in mesh.cell_data["octree:leaf"]])
    assert leaf.max() < len(plain.cells[0].data)


@pytest.mark.parametrize("transition", ["polyhedra", "templates"])
def test_transition_needs_a_balanced_tree(transition):
    with pytest.raises(ValueError, match="balanced"):
        octree(
            _sphere(), bounds=BOX, max_level=3, balance="none", transition=transition
        )


@pytest.mark.parametrize("transition", ["polyhedra", "templates"])
def test_sphere_transition_conforms_and_conserves_volume(transition):
    mesh = octree(_sphere(), bounds=BOX, max_level=4, transition=transition)
    _check_conforming(mesh, [-1.2] * 3, [1.2] * 3, 2.4**3)
    # data is carried onto every cell of the leaf it was cut from
    for name in ("octree:level", "octree:cut", "octree:leaf"):
        assert len(mesh.cell_data[name]) == len(mesh.cells)


@pytest.mark.parametrize("transition", ["polyhedra", "templates"])
def test_quadtree_transition_conforms(transition):
    mesh, rep = octree(
        _circle(n=128),
        bounds=[-1.2, -1.2, 0, 1.2, 1.2, 0],
        max_level=5,
        transition=transition,
        return_report=True,
    )
    kinds = {cb.type for cb in mesh.cells}
    assert kinds <= (
        {"quad", "polygon"} if transition == "polyhedra" else {"quad", "triangle"}
    )
    assert len(kinds) == 2
    _check_conforming(mesh, [-1.2, -1.2], [1.2, 1.2], 2.4**2)
    assert rep["num_hanging"] == 0


def test_snap_moves_outer_nodes_onto_the_sphere_without_inverting_a_leaf():
    kw = dict(bounds=BOX, max_level=4, fill="inside")
    plain = octree(_sphere(nu=48, nv=24), **kw)
    snapped, rep = octree(_sphere(nu=48, nv=24), snap=True, return_report=True, **kw)
    assert rep["num_snapped"] > 0
    pts = np.asarray(snapped.points)
    ref = np.asarray(plain.points)
    assert pts.shape == ref.shape
    moved = np.any(pts != ref, axis=1)
    assert moved.sum() == rep["num_snapped"]
    r = np.linalg.norm(pts[moved], axis=1)
    r0 = np.linalg.norm(ref[moved], axis=1)
    assert np.abs(r - 1.0).mean() < np.abs(r0 - 1.0).mean() * 0.1
    cb = snapped.cells[0]
    assert np.all(_cell_volume(cb, pts) > 0.0)


def test_snap_needs_a_boundary():
    with pytest.raises(ValueError, match="snap"):
        octree(_sphere(), bounds=BOX, max_level=3, snap=True)


def test_snap_in_2d_lands_on_the_polyline():
    mesh, rep = octree(
        _circle(n=256),
        bounds=[-1.2, -1.2, 0, 1.2, 1.2, 0],
        max_level=6,
        fill="inside",
        snap=True,
        return_report=True,
    )
    pts = np.asarray(mesh.points)
    on = np.isclose(np.linalg.norm(pts[:, :2], axis=1), 1.0, atol=1e-3)
    assert rep["num_snapped"] > 0 and on.sum() >= rep["num_snapped"] // 2
    assert np.all(_cell_volume(mesh.cells[0], pts) > 0.0)


def test_surface_level_and_criteria():
    sph = _sphere(nu=48, nv=24)
    kw = dict(bounds=BOX, max_level=5)
    n_all = octree(sph, **kw).cells[0].data.shape[0]
    n_low = octree(sph, surface_level=3, **kw).cells[0].data.shape[0]
    assert n_low < n_all
    # every skin triangle flagged: the criteria carry the tree to max_level again
    n_curv = (
        octree(sph, surface_level=3, curvature_threshold=0.1, **kw)
        .cells[0]
        .data.shape[0]
    )
    assert n_curv == n_all
    # a threshold nothing reaches changes nothing
    assert (
        octree(sph, surface_level=3, curvature_threshold=1e6, **kw)
        .cells[0]
        .data.shape[0]
        == n_low
    )

    sph.point_data["h"] = np.full(len(sph.points), 0.2)
    sized = octree(sph, surface_level=2, size_array="h", **kw)
    level = np.asarray(sized.cell_data["octree:level"][0])
    cut = np.asarray(sized.cell_data["octree:cut"][0])
    assert np.all(level[cut == 1] == 4)  # 2.4 / 2**4 = 0.15 <= 0.2 < 0.3
    with pytest.raises(ValueError, match="size_array"):
        octree(sph, size_array="nope", **kw)


def test_feature_angle_refines_the_edges_of_a_cube():
    from .test_voxelize import _box

    cube = _box(1.0)
    kw = dict(bounds=[-1.5] * 3 + [1.5] * 3, max_level=4, surface_level=2)
    n_flat = octree(cube, **kw).cells[0].data.shape[0]
    n_edge = octree(cube, feature_angle=30.0, **kw).cells[0].data.shape[0]
    n_full = octree(cube, bounds=kw["bounds"], max_level=4).cells[0].data.shape[0]
    assert n_flat < n_edge < n_full


@pytest.mark.parametrize("plane,axes", [("xz", (0, 2, 1)), ("yz", (1, 2, 0))])
def test_planes(plane, axes):
    xy = _circle(n=128)
    pts = np.asarray(xy.points)
    moved = np.zeros_like(pts)
    moved[:, list(axes)] = pts  # frame -> world
    skin = meshioplusplus.Mesh(moved, xy.cells)
    ref = octree(xy, max_level=4, fill="inside")
    got = octree(skin, max_level=4, fill="inside")  # auto picks the zero-extent axis
    forced = octree(skin, max_level=4, fill="inside", plane=plane)
    assert got.cells[0].type == "quad"
    flat = [k for k in range(3) if k not in axes[:2]][0]
    assert np.all(np.asarray(got.points)[:, flat] == 0.0)
    assert np.asarray(got.points).tobytes() == np.asarray(forced.points).tobytes()
    assert len(got.cells[0].data) == len(ref.cells[0].data)
    # pixels follow the same rule
    pix = voxelize(skin, resolution=(8, 8), fill="inside")
    assert np.all(np.asarray(pix.points)[:, flat] == 0.0)


def test_planes_must_be_known():
    with pytest.raises(ValueError, match="plane"):
        octree(_circle(), max_level=2, plane="zz")


def test_numpy_reference_declines_core_only_options():
    from meshioplusplus._octree import _octree_py

    with pytest.raises(NotImplementedError, match="compiled core"):
        _octree_py(
            _circle(),
            [1, 1, 1],
            0,
            2,
            None,
            0.0,
            0.0,
            0,
            0.0,
            "face",
            "all",
            False,
            False,
            1000,
            "pseudonormal",
            transition="templates",
        )
