"""tetrahedralize (roadmap 6.7): tetrahedra filling a closed triangle surface.

Through the Python wrapper: with Netgen, the sphere's CAD skin filled with its
triangles kept and its enclosed volume matched and a box in a box giving two
regions; with the native octree backend, always built, the same skin filled
with its boundary fitted within the reported distance; the report, the
refusals by name, determinism, and the fallback of a core built without
Netgen. The CI ``netgen`` job runs the Netgen tests; the algorithm itself is
tested in the core (``tests/cpp/test_tetrahedralize.cpp``).
"""

import pathlib

import numpy as np
import pytest

import meshioplusplus as mio

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - the numpy-only install
    _core = None

HAS_NETGEN = bool(getattr(_core, "__has_netgen__", False))

needs_core = pytest.mark.skipif(
    not hasattr(_core, "tetrahedralize"), reason="core predates tetrahedralize"
)
needs_netgen = pytest.mark.skipif(
    not HAS_NETGEN, reason="core built without MESHIOPLUSPLUS_WITH_NETGEN"
)

CAD = pathlib.Path(__file__).parent / "meshes" / "cad_json"

_REPORT_KEYS = {
    "num_tetrahedra",
    "num_points",
    "num_new_points",
    "num_regions",
    "num_boundary_triangles",
    "boundary_preserved",
    "hausdorff_to_surface",
    "min_dihedral_deg",
    "max_dihedral_deg",
    "volume",
    "backend",
    "warnings",
}

_BOX_TRIANGLES = np.array(
    [
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],
        [3, 0, 4],
        [3, 4, 7],
    ]
)


def _box(lo, hi, inward=False):
    corners = np.array(
        [
            [lo, lo, lo],
            [hi, lo, lo],
            [hi, hi, lo],
            [lo, hi, lo],
            [lo, lo, hi],
            [hi, lo, hi],
            [hi, hi, hi],
            [lo, hi, hi],
        ],
        dtype=float,
    )
    tris = _BOX_TRIANGLES[:, [0, 2, 1]] if inward else _BOX_TRIANGLES
    return corners, tris


def _boxes(*specs):
    points, tris = [], []
    for spec in specs:
        p, t = _box(*spec)
        tris.append(t + sum(len(q) for q in points))
        points.append(p)
    return mio.Mesh(np.vstack(points), [("triangle", np.vstack(tris))])


def _sphere_skin(h=0.3):
    return mio.mesh_cad_surface(mio.read(CAD / "sphere.cad.json"), h)


def _triangles(mesh):
    return np.vstack([c.data for c in mesh.cells if c.type == "triangle"])


def _enclosed_volume(mesh):
    """The volume a closed, outward triangle surface encloses (divergence
    theorem)."""
    p = np.asarray(mesh.points, dtype=float)
    t = _triangles(mesh)
    a, b, c = p[t[:, 0]], p[t[:, 1]], p[t[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


@needs_core
@needs_netgen
def test_the_sphere_is_filled_keeping_its_skin():
    skin = _sphere_skin()
    out, report = mio.tetrahedralize(skin, max_size=0.3, return_report=True)
    assert set(report) == _REPORT_KEYS
    assert report["boundary_preserved"] and report["hausdorff_to_surface"] == 0.0
    assert report["backend"] == "netgen" and report["num_regions"] == 1
    assert report["num_tetrahedra"] > 0 and report["num_new_points"] > 0
    assert 0.0 < report["min_dihedral_deg"] < report["max_dihedral_deg"] < 180.0
    # The tetrahedra fill exactly what the skin encloses, which is within the
    # skin's chordal error of the sphere.
    assert report["volume"] == pytest.approx(_enclosed_volume(skin), rel=1e-9)
    assert report["volume"] == pytest.approx(4 / 3 * np.pi, rel=0.05)
    assert report["volume"] < 4 / 3 * np.pi
    # The input's points first and unchanged, its triangles as given.
    n = len(skin.points)
    assert report["num_points"] == len(out.points) == n + report["num_new_points"]
    np.testing.assert_array_equal(out.points[:n], skin.points)
    np.testing.assert_array_equal(_triangles(out), _triangles(skin))
    assert report["num_boundary_triangles"] == len(_triangles(skin))
    tets = [c for c in out.cells if c.type == "tetra"]
    assert len(tets) == 1 and len(tets[0].data) == report["num_tetrahedra"]
    assert "tetrahedralize:region" in out.cell_data
    # The CAD faces ride along on the triangles.
    assert "cad:face" in out.cell_data


@needs_core
@needs_netgen
@pytest.mark.parametrize("inward", [False, True])
def test_a_box_in_a_box_gives_two_regions(inward):
    surface = _boxes((0.0, 3.0), (1.0, 2.0, inward))
    out, report = mio.tetrahedralize(surface, max_size=0.5, return_report=True)
    assert report["num_regions"] == 2 and report["boundary_preserved"]
    assert report["volume"] == pytest.approx(27.0, rel=1e-12)
    regions = np.concatenate(
        [
            np.asarray(out.cell_data["tetrahedralize:region"][b])
            for b, c in enumerate(out.cells)
            if c.type == "tetra"
        ]
    )
    assert set(regions.tolist()) == {0, 1}


@needs_core
@needs_netgen
def test_two_runs_agree():
    skin = _sphere_skin(0.4)
    a = mio.tetrahedralize(skin, max_size=0.3)
    b = mio.tetrahedralize(skin, max_size=0.3)
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for ca, cb in zip(a.cells, b.cells):
        np.testing.assert_array_equal(ca.data, cb.data)


@needs_core
@needs_netgen
def test_a_smaller_size_gives_more_tetrahedra():
    skin = _sphere_skin(0.4)
    _, coarse = mio.tetrahedralize(skin, return_report=True)
    _, fine = mio.tetrahedralize(skin, max_size=0.15, return_report=True)
    assert fine["num_tetrahedra"] > coarse["num_tetrahedra"]
    with pytest.raises(ValueError, match="more than 10 tetrahedra .max_cells."):
        mio.tetrahedralize(skin, max_size=0.15, max_cells=10)


@needs_core
def test_refusals_are_by_name():
    box = _boxes((0.0, 1.0))
    for kwargs, name in [
        ({"backend": "tetgen"}, "backend"),
        ({"max_size": -1.0}, "max_size"),
        ({"max_volume": float("nan")}, "max_volume"),
        ({"grading": 2.0}, "grading"),
        ({"optimize_steps": 101}, "optimize_steps"),
        ({"optimize_steps": 1.5}, "optimize_steps"),
        ({"max_cells": 0}, "max_cells"),
    ]:
        with pytest.raises(ValueError, match=name):
            mio.tetrahedralize(box, **kwargs)
    # An open surface: a triangle left out.
    open_box = mio.Mesh(box.points, [("triangle", box.cells[0].data[:-1])])
    with pytest.raises(ValueError, match="open edges"):
        mio.tetrahedralize(open_box)
    quads = mio.Mesh(np.eye(4, 3), [("quad", [[0, 1, 2, 3]])])
    with pytest.raises(ValueError, match="quad"):
        mio.tetrahedralize(quads)


@needs_core
@pytest.mark.skipif(HAS_NETGEN, reason="Netgen is built in")
def test_without_netgen_auto_falls_back_and_netgen_is_named():
    box = _boxes((0.0, 1.0))
    with pytest.warns(RuntimeWarning, match="MESHIOPLUSPLUS_WITH_NETGEN"):
        _, report = mio.tetrahedralize(box, return_report=True)
    assert report["backend"] == "octree" and not report["boundary_preserved"]
    with pytest.raises(NotImplementedError, match="MESHIOPLUSPLUS_WITH_NETGEN"):
        mio.tetrahedralize(box, backend="netgen")


@needs_core
def test_auto_takes_netgen_when_built_in():
    _, report = mio.tetrahedralize(_boxes((0.0, 1.0)), return_report=True)
    assert report["backend"] == ("netgen" if HAS_NETGEN else "octree")


@needs_core
def test_the_octree_backend_fits_the_sphere():
    h = 0.3
    skin = _sphere_skin(h)
    out, report = mio.tetrahedralize(skin, backend="octree", return_report=True)
    assert set(report) == _REPORT_KEYS
    assert report["backend"] == "octree" and report["warnings"] == []
    assert not report["boundary_preserved"]
    assert 0.0 < report["hausdorff_to_surface"] < 0.1 * h
    assert report["num_new_points"] == report["num_points"] == len(out.points)
    assert 0.0 < report["min_dihedral_deg"] < report["max_dihedral_deg"] < 180.0
    # The volume differs from the skin's by at most the boundary's distance
    # to it over the skin's area.
    p = np.asarray(skin.points, dtype=float)
    t = _triangles(skin)
    area = (
        0.5
        * np.linalg.norm(
            np.cross(p[t[:, 1]] - p[t[:, 0]], p[t[:, 2]] - p[t[:, 0]]), axis=1
        ).sum()
    )
    gap = abs(report["volume"] - _enclosed_volume(skin))
    assert gap <= report["hausdorff_to_surface"] * area
    # The triangle block is the tetrahedra's boundary, closed and outward.
    boundary = mio.Mesh(out.points, [("triangle", _triangles(out))])
    assert report["num_boundary_triangles"] == len(_triangles(out))
    assert _enclosed_volume(boundary) == pytest.approx(report["volume"], rel=1e-9)
    # The CAD faces reach the fitted boundary, one Side region per face.
    faces = np.concatenate(
        [
            np.asarray(out.cell_data["cad:face"][b])
            for b, c in enumerate(out.cells)
            if c.type == "triangle"
        ]
    )
    assert (faces >= 0).all()
    sides = [r.name for r in out.regions if r.name.startswith("face:")]
    assert sorted(sides) == [f"face:{v}" for v in sorted(set(faces.tolist()))]
    again = mio.tetrahedralize(skin, backend="octree")
    np.testing.assert_array_equal(again.points, out.points)


@needs_core
def test_the_octree_backend_keeps_two_regions():
    surface = _boxes((0.0, 3.0), (1.0, 2.0))
    with pytest.warns(RuntimeWarning, match="regions"):
        out, report = mio.tetrahedralize(
            surface, backend="octree", max_size=0.25, return_report=True
        )
    assert report["num_regions"] == 2
    assert report["volume"] == pytest.approx(27.0, rel=0.02)
    regions = np.concatenate(
        [
            np.asarray(out.cell_data["tetrahedralize:region"][b])
            for b, c in enumerate(out.cells)
            if c.type == "tetra"
        ]
    )
    assert set(regions.tolist()) == {0, 1}
    with pytest.raises(ValueError, match="max_cells"):
        mio.tetrahedralize(surface, backend="octree", max_cells=10)
