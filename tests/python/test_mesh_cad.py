"""mesh_cad_surface: conforming surface meshes of CAD B-reps, and
triangulate_polygon_2d.

The CAD JSON fixtures through the Python wrapper: the closed solids
(``unit_cube``, ``capped_cylinder``, ``sphere``) come out watertight with
their area and enclosed volume converging at second order, the trimmed plate
keeps its hole, the Scordelis-Lo roof (an untrimmed patch whose parameter
scales differ fortyfold) and the two-patch beam (a coupling edge shared by
two faces, and an embedded edge meshed as mesh edges) mesh conformingly;
the report, the refusals before the core, strict mode, a bare B-rep dict,
round trips through VTU and Gmsh, both CLI verbs and the pipeline step that
meshes the model its Input reads. With OpenCASCADE, a STEP solid.
The algorithm itself is tested in the core (``tests/cpp/test_mesh_cad.cpp``,
``test_cdt2d.cpp``, ``test_edge_sampler.cpp``, ``test_predicates.cpp``).
"""

import copy
import json
import math
import pathlib

import numpy as np
import pytest

import meshioplusplus as mio

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - the numpy-only install
    _core = None

needs_core = pytest.mark.skipif(
    not hasattr(_core, "mesh_cad_surface"), reason="core predates mesh_cad_surface"
)

CAD = pathlib.Path(__file__).parent / "meshes" / "cad_json"


def _read(name):
    return mio.read(CAD / name)


def _measure(mesh):
    """Area and enclosed volume (divergence theorem) of the triangle block."""
    tri = [c for c in mesh.cells if c.type == "triangle"][0].data
    p = mesh.points
    a, b, c = p[tri[:, 0]], p[tri[:, 1]], p[tri[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()
    volume = np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0
    return area, volume


def _watertight(report):
    assert report["num_failed"] == 0, report["faces"]
    assert report["num_open_edges"] == 0
    assert report["num_non_manifold_edges"] == 0
    assert report["num_misoriented_edges"] == 0


@needs_core
@pytest.mark.parametrize(
    "name, area, volume",
    [
        ("unit_cube.cad.json", 6.0, 1.0),
        ("capped_cylinder.cad.json", 2 * math.pi * 2 + 2 * math.pi, 2 * math.pi),
        ("sphere.cad.json", 4 * math.pi, 4 * math.pi / 3),
    ],
)
def test_closed_solids_are_watertight_and_converge(name, area, volume):
    model = _read(name)
    errors = []
    for h in (0.2, 0.1):
        out, report = mio.mesh_cad_surface(model, h, return_report=True)
        _watertight(report)
        a, v = _measure(out)
        assert v > 0.0  # outward
        errors.append((abs(a - area), abs(v - volume)))
        # Every point on the exact surface.
        if name.startswith("sphere"):
            np.testing.assert_allclose(
                np.linalg.norm(out.points, axis=1), 1.0, atol=1e-12
            )
    if name.startswith("unit_cube"):
        for ea, ev in errors:
            assert ea < 1e-12 and ev < 1e-12
    else:
        (a0, v0), (a1, v1) = errors
        assert a1 < a0 / 2.5 and v1 < v0 / 2.5  # about 4: second order


@needs_core
def test_trimmed_plate_keeps_its_hole():
    model = _read("plate_with_hole.cad.json")
    exact = 16.0 - math.pi
    errors = []
    for h in (0.2, 0.1):
        out, report = mio.mesh_cad_surface(model, h, return_report=True)
        assert report["num_failed"] == 0
        assert report["min_angle_deg"] > 10.0
        assert report["max_edge"] <= 1.6 * h
        # The hole is an inscribed polygon: the plate is slightly larger.
        a, _ = _measure(out)
        assert a > exact
        errors.append(a - exact)
    assert errors[1] < errors[0] / 2.5


@needs_core
def test_scordelis_lo_roof_untrimmed_and_anisotropic():
    model = _read("cylindrical_roof.cad.json")
    exact = 25.0 * math.radians(80.0) * 50.0
    errors = []
    for h in (2.0, 1.0):
        out, report = mio.mesh_cad_surface(model, h, return_report=True)
        assert report["num_failed"] == 0
        assert report["min_angle_deg"] > 10.0
        assert report["max_edge"] <= 1.5 * h
        np.testing.assert_allclose(
            np.hypot(out.points[:, 0], out.points[:, 2]), 25.0, rtol=1e-12
        )
        errors.append(abs(_measure(out)[0] - exact))
    assert errors[1] < errors[0] / 2.5


@needs_core
def test_two_patch_beam_shares_its_coupling_edge():
    out, report = mio.mesh_cad_surface(
        _read("two_patch_beam.cad.json"), 0.2, return_report=True
    )
    assert report["num_faces"] == 2 and report["num_failed"] == 0
    assert report["num_non_manifold_edges"] == 0
    assert report["num_misoriented_edges"] == 0
    # Some points belong to triangles of both faces: the shared edge.
    tri = out.cells[0].data
    face = out.cell_data["cad:face"][0]
    on0 = set(tri[face == 0].ravel())
    on1 = set(tri[face == 1].ravel())
    assert len(on0 & on1) > 2
    # Roadmap 6.1.1: the first face's embedded edge is meshed, not warned
    # about: its line cells are edges of two triangles of that face.
    assert not any("embedded" in w for w in report["warnings"])
    uses = {}
    for t in tri[face == 0]:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            key = (min(a, b), max(a, b))
            uses[key] = uses.get(key, 0) + 1
    lines = out.cells[1].data
    interior = [(a, b) for a, b in lines if uses.get((min(a, b), max(a, b)), 0) == 2]
    assert len(interior) >= 5  # the curve is about 3 long; h = 0.2


@needs_core
def test_closed_edge_started_elsewhere_is_rotated():
    # Roadmap 6.1.2: the top cap's circle starts a quarter turn after the
    # band's rim; the cap is a loop of its own, so it is rotated onto the
    # rim's corner and the solid closes.
    out, report = mio.mesh_cad_surface(
        _read("capped_cylinder_offset.cad.json"), 0.2, return_report=True
    )
    _watertight(report)
    assert report["max_edge_gap"] < 1e-12
    area, volume = _measure(out)
    assert abs(volume - 2 * math.pi) < 0.02 * 2 * math.pi
    assert abs(area - 6 * math.pi) < 0.02 * 6 * math.pi


@needs_core
def test_loops_crossing_at_a_coarse_size_are_split():
    # Roadmap 6.1.3: the hole comes within 0.02 of the rim, and at h = 1.1 a
    # rim chord passes inside the hole's sample. The pieces involved are
    # halved until the loops do not cross, and the plate is meshed.
    model = _read("plate_with_near_hole.cad.json")
    out, report = mio.mesh_cad_surface(model, 1.1, return_report=True)
    assert report["num_failed"] == 0, report["faces"]
    assert report["num_non_manifold_edges"] == 0
    area, _ = _measure(out)
    assert 0.9 * 3 * math.pi < area < 3 * math.pi
    # The same at a size that never crossed: nothing to split.
    fine = mio.mesh_cad_surface(model, 0.1)
    assert 0.99 * 3 * math.pi < _measure(fine)[0] < 3 * math.pi


@needs_core
def test_report_layout_and_regions():
    out, report = mio.mesh_cad_surface(
        _read("unit_cube.cad.json"), 0.25, return_report=True
    )
    for key in (
        "num_points",
        "num_triangles",
        "num_edge_segments",
        "num_faces",
        "num_failed",
        "num_edges",
        "tolerance",
        "min_angle_deg",
        "max_edge",
        "max_edge_gap",
        "faces",
        "warnings",
    ):
        assert key in report
    assert report["num_points"] == len(out.points)
    assert [c.type for c in out.cells] == ["triangle", "line"]
    assert report["num_edges"] == 12
    assert len(report["faces"]) == 6
    assert {f["status"] for f in report["faces"]} == {"ok"}
    assert set(out.cell_data) == {"cad:face", "cad:edge"}
    face = out.cell_data["cad:face"]
    assert set(face[0]) == set(range(6)) and set(face[1]) == {-1}
    edge = out.cell_data["cad:edge"]
    assert set(edge[0]) == {-1} and min(edge[1]) >= 0
    names = {r.name for r in out.regions}
    assert names == {f"face:{k}" for k in range(6)}
    bare = mio.mesh_cad_surface(
        _read("unit_cube.cad.json"), 0.25, edges=False, record_face_ids=False
    )
    assert [c.type for c in bare.cells] == ["triangle"]
    assert not bare.cell_data and not bare.regions


@needs_core
def test_bare_dict_and_determinism():
    model = _read("capped_cylinder.cad.json")
    a = mio.mesh_cad_surface(model, 0.2)
    b = mio.mesh_cad_surface(model.info["brep"], 0.2)
    assert np.array_equal(a.points, b.points)
    assert np.array_equal(a.cells[0].data, b.cells[0].data)


@needs_core
def test_a_failed_face_is_reported_and_strict_raises():
    brep = copy.deepcopy(_read("unit_cube.cad.json").info["brep"])
    brep["faces"][2]["patch"] = 99
    out, report = mio.mesh_cad_surface(brep, 0.25, return_report=True)
    assert report["num_failed"] == 1
    bad = report["faces"][2]
    assert bad["status"] == "unsupported" and bad["reason"]
    assert report["num_open_edges"] > 0  # the hole is visible
    # The refusal says why; the edges naming the face's trims add nothing.
    assert not any("does not exist" in w for w in report["warnings"])
    assert 2 not in set(out.cell_data["cad:face"][0])
    with pytest.raises(RuntimeError, match="face 2"):
        mio.mesh_cad_surface(brep, 0.25, strict=True)


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"target_size": 0.0}, "target_size"),
        ({"target_size": -1.0}, "target_size"),
        ({"target_size": float("nan")}, "target_size"),
        ({"target_size": 0.1, "max_chordal_deviation": -1.0}, "max_chordal"),
        ({"target_size": 0.1, "min_points_per_edge": -1}, "min_points_per_edge"),
        ({"target_size": 0.1, "tolerance": -1.0}, "tolerance"),
        ({"target_size": 0.1, "max_points": 0}, "max_points"),
        ({"target_size": 0.1, "grading": 0.5}, "grading"),
        ({"target_size": 0.1, "grading": 1.0}, "grading"),
        ({"target_size": 0.1, "grading": float("nan")}, "grading"),
        ({"target_size": 0.1, "curvature_deviation": -1.0}, "curvature_deviation"),
        ({"target_size": 0.1, "cells_across_gaps": -1.0}, "cells_across_gaps"),
        ({"target_size": 0.1, "pole_rings": -1}, "pole_rings"),
        ({"target_size": 0.1, "pole_rings": 9}, "pole_rings"),
    ],
)
def test_refusals_before_the_core(kwargs, match):
    model = _read("unit_cube.cad.json")
    with pytest.raises(ValueError, match=match):
        mio.mesh_cad_surface(model, **kwargs)


def test_a_mesh_without_a_brep_is_refused():
    with pytest.raises(ValueError, match="no B-rep"):
        mio.mesh_cad_surface(mio.rectangle(), 0.1)
    with pytest.raises(ValueError, match="not a B-rep"):
        mio.mesh_cad_surface({"patches": []}, 0.1)


@needs_core
@pytest.mark.parametrize("ext, fmt", [(".vtu", "vtu"), (".msh", "gmsh")])
def test_round_trip(tmp_path, ext, fmt):
    out = mio.mesh_cad_surface(_read("sphere.cad.json"), 0.3)
    path = tmp_path / f"sphere{ext}"
    mio.write(path, out, file_format=fmt)
    back = mio.read(path, file_format=fmt)
    if fmt == "vtu":
        np.testing.assert_array_equal(back.points, out.points)
        tri = np.concatenate([c.data for c in back.cells if c.type == "triangle"])
        np.testing.assert_array_equal(
            np.sort(tri, axis=0), np.sort(out.cells[0].data, axis=0)
        )
        return

    # Gmsh 4.1 writes one entity per face (roadmap 6.1.2) and groups $Nodes by
    # entity, so the points come back grouped by face: the same triangles,
    # point for point, and every face its own physical surface.
    def triangles(mesh):
        tri = np.concatenate([c.data for c in mesh.cells if c.type == "triangle"])
        xyz = mesh.points[tri]
        xyz = np.array([t[np.lexsort(t.T[::-1])] for t in xyz]).reshape(len(tri), 9)
        return xyz[np.lexsort(xyz.T[::-1])]

    np.testing.assert_array_equal(triangles(back), triangles(out))
    faces = sorted(r.name for r in back.regions if r.name.startswith("face:"))
    assert faces == sorted(r.name for r in out.regions)


def _tet_volume(mesh):
    v = 0.0
    for c in mesh.cells:
        if c.type == "tetra":
            p = mesh.points[c.data]
            a, b, d = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]
            v += np.abs(np.einsum("ij,ij->i", a, np.cross(b, d))).sum() / 6.0
    return v


@needs_core
@pytest.mark.parametrize(
    "name, volume, dihedrals",
    [
        ("sphere.cad.json", 4 * math.pi / 3, (10.7, 164.8)),
        ("unit_cube.cad.json", 1.0, (35.0, 126.0)),
    ],
)
def test_octree_mesh_fills_the_skin(name, volume, dihedrals):
    """Roadmap 6.7.1: a closed mesh_cad_surface result is the skin octree_mesh
    fills with tetrahedra; under uniform skin grading at the surface's size the
    volume follows the skin's, the boundary stays within a few hundredths of h
    of it and the dihedral angles keep inside the stuffing bound."""
    h = 0.2
    skin = mio.mesh_cad_surface(_read(name), h, min_angle=25)
    vol, report = mio.octree_mesh(
        skin, skin_grading="uniform", surface_size=h, return_report=True
    )
    _, skin_volume = _measure(skin)
    assert abs(_tet_volume(vol) - skin_volume) < 0.005 * volume  # 0.30% on the sphere
    assert abs(_tet_volume(vol) - volume) < 0.02 * volume
    assert 0.0 <= report["hausdorff_to_skin"] < 0.06 * h
    assert report["min_dihedral"] > dihedrals[0]
    assert report["max_dihedral"] < dihedrals[1]


@needs_core
def test_triangulate_polygon_2d():
    outer = [[0, 0], [4, 0], [4, 4], [0, 4]]
    hole = [[1, 1], [1, 3], [3, 3], [3, 1]]
    m = mio.triangulate_polygon_2d(outer, [hole])
    assert len(m.points) == 8 and len(m.cells[0].data) == 8
    assert _measure(m)[0] == pytest.approx(12.0, rel=1e-12)
    fine = mio.triangulate_polygon_2d(outer, [hole], target_size=0.25)
    assert _measure(fine)[0] == pytest.approx(12.0, rel=1e-12)
    assert np.all(fine.points[:, 2] == 0.0)
    with pytest.raises(ValueError, match="self-intersecting"):
        mio.triangulate_polygon_2d([[0, 0], [1, 1], [1, 0], [0, 1]])
    with pytest.raises(ValueError, match="target_size"):
        mio.triangulate_polygon_2d(outer, target_size=-1.0)
    with pytest.raises(ValueError, match=r"\(n, 2\)"):
        mio.triangulate_polygon_2d([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    # Beyond the range the exact predicates hold in (a fuzz finding).
    with pytest.raises(ValueError, match="exact predicates"):
        mio.triangulate_polygon_2d([[5e-299, -8.9e304], [-8.9e304, 1e-323], [0, 0]])


def _min_angles(mesh):
    """The smallest angle (degrees) of every triangle."""
    tri = [c for c in mesh.cells if c.type == "triangle"][0].data
    p = mesh.points[tri]
    out = np.full(len(tri), 180.0)
    for k in range(3):
        u = p[:, (k + 1) % 3] - p[:, k]
        v = p[:, (k + 2) % 3] - p[:, k]
        ang = np.degrees(
            np.arctan2(
                np.linalg.norm(np.cross(u, v), axis=1), np.einsum("ij,ij->i", u, v)
            )
        )
        out = np.minimum(out, ang)
    return out


@needs_core
@pytest.mark.parametrize(
    "name", ["plate_with_hole.cad.json", "capped_cylinder.cad.json", "sphere.cad.json"]
)
def test_refinement_meets_the_bound_or_reports_the_corner(name):
    """Roadmap 6.2.7 on the CAD JSON fixtures: every triangle meets 25 degrees
    except within the reach of a reported pole; 0 keeps the unrefined mesh."""
    model = _read(name)
    plain = mio.mesh_cad_surface(model, 0.2)
    zero = mio.mesh_cad_surface(model, 0.2, min_angle=0.0)
    np.testing.assert_array_equal(plain.points, zero.points)
    out, report = mio.mesh_cad_surface(model, 0.2, min_angle=25.0, return_report=True)
    assert report["num_failed"] == 0
    poles = [c for f in report["faces"] for c in f["corners"] if c["pole"]]
    assert all(not c["pole"] or c["reach"] >= 0.2 for c in poles)
    tri = [c for c in out.cells if c.type == "triangle"][0].data
    below = _min_angles(out) < 25.0 - 1e-9
    for t in np.nonzero(below)[0]:
        x = out.points[tri[t]]
        assert any(
            np.min(np.linalg.norm(x - np.asarray(c["x"]), axis=1)) <= c["reach"] + 1e-12
            for c in poles
        ), "a triangle under the bound touches no reported corner"
    assert sum(report["angle_histogram"]) == 3 * report["num_triangles"]
    assert report["num_faces_below_bound"] == sum(
        f["below_bound"] for f in report["faces"]
    )
    if name.startswith("sphere"):
        assert len(poles) == 2
    else:
        assert not poles and report["num_faces_below_bound"] == 0


@needs_core
@pytest.mark.parametrize("order", [2, 3])
def test_curved_elements_put_every_node_on_the_sphere(order):
    """Roadmap 6.6: order 2 and 3 give triangle6/line3 and triangle10/line4
    whose corners are the linear mesh and whose every node is on the sphere."""
    model = _read("sphere.cad.json")
    linear = mio.mesh_cad_surface(model, 0.3)
    out, report = mio.mesh_cad_surface(model, 0.3, order=order, return_report=True)
    _watertight(report)
    assert [c.type for c in out.cells] == (
        ["triangle6", "line3"] if order == 2 else ["triangle10", "line4"]
    )
    assert report["num_curved_straightened"] == 0
    np.testing.assert_array_equal(out.points[: len(linear.points)], linear.points)
    np.testing.assert_array_equal(out.cells[0].data[:, :3], linear.cells[0].data)
    np.testing.assert_allclose(np.linalg.norm(out.points, axis=1), 1.0, atol=1e-12)
    for bad in (0, 4, True):
        with pytest.raises(ValueError, match="order"):
            mio.mesh_cad_surface(model, 0.3, order=bad)


@needs_core
@pytest.mark.parametrize("order", [2, 3])
@pytest.mark.parametrize("ext", [".mesh", ".meshb", ".msh"])
def test_curved_elements_survive_a_file(tmp_path, order, ext):
    """The curved elements round-trip through Medit (its P2/P3 keywords) and
    Gmsh with their node order."""
    out = mio.mesh_cad_surface(_read("sphere.cad.json"), 0.5, order=order)
    path = tmp_path / ("curved" + ext)
    mio.write(path, out, file_format="gmsh" if ext == ".msh" else None)
    back = mio.read(path, file_format="gmsh" if ext == ".msh" else None)
    if ext == ".msh":
        # One entity per face (roadmap 6.1.2): Gmsh 4.1 groups $Nodes by
        # entity, so the points come back grouped by face. Each element keeps
        # its nodes, in their order, point for point.
        for kind in {c.type for c in out.cells}:
            got = np.concatenate(
                [back.points[c.data] for c in back.cells if c.type == kind]
            )
            want = np.concatenate(
                [out.points[c.data] for c in out.cells if c.type == kind]
            )
            got, want = got.reshape(len(got), -1), want.reshape(len(want), -1)
            np.testing.assert_array_equal(
                got[np.lexsort(got.T[::-1])], want[np.lexsort(want.T[::-1])]
            )
        return
    by_type = {c.type: c.data for c in back.cells}
    for c in out.cells:
        np.testing.assert_array_equal(by_type[c.type], c.data)
    np.testing.assert_allclose(back.points, out.points, rtol=0, atol=1e-15)


@needs_core
@needs_core
@pytest.mark.parametrize("deviation", [0.002, 0.0005])
def test_curvature_sizing_keeps_the_chordal_error_under_the_bound(deviation):
    """Roadmap 6.3.2: curvature_deviation sizes the cylinder's side so no
    triangle sits farther than it from the surface; the flat caps stay near
    target_size, and the result is still watertight."""
    model = _read("capped_cylinder.cad.json")
    out, report = mio.mesh_cad_surface(
        model, 0.2, curvature_deviation=deviation, return_report=True
    )
    _watertight(report)
    tri, p = out.cells[0].data, out.points
    radius = np.hypot(p[:, 0], p[:, 1])
    side = np.all(np.abs(radius[tri] - 1.0) < 1e-9, axis=1)
    centroid = p[tri[side]].mean(axis=1)
    assert (1.0 - np.hypot(centroid[:, 0], centroid[:, 1])).max() <= deviation
    assert report["max_edge"] <= 1.5 * 0.2
    uniform = mio.mesh_cad_surface(model, 0.2, return_report=True)[1]
    assert report["num_triangles"] > uniform["num_triangles"]


@needs_core
def test_sizing_defaults_and_determinism():
    """Every sizing option at its default is the unsized mesh byte for byte; a
    sized run is the same twice, and a stronger grading limit spends more
    triangles."""
    model = _read("capped_cylinder.cad.json")
    base = mio.mesh_cad_surface(model, 0.2)
    same = mio.mesh_cad_surface(
        model, 0.2, curvature_deviation=0.0, cells_across_gaps=0.0, grading=0.0
    )
    np.testing.assert_array_equal(base.points, same.points)
    np.testing.assert_array_equal(base.cells[0].data, same.cells[0].data)
    a = mio.mesh_cad_surface(model, 0.2, curvature_deviation=0.002, grading=1.2)
    b = mio.mesh_cad_surface(model, 0.2, curvature_deviation=0.002, grading=1.2)
    np.testing.assert_array_equal(a.points, b.points)
    np.testing.assert_array_equal(a.cells[0].data, b.cells[0].data)
    loose = mio.mesh_cad_surface(model, 0.2, curvature_deviation=0.002, grading=2.0)
    assert len(a.cells[0].data) > len(loose.cells[0].data)


@needs_core
def test_pole_rings_lift_the_sphere_to_the_bound():
    """Roadmap 6.4.2: without rings the sphere's pole fans stay below 25
    degrees and are reported; two rings grade them and the bound is met."""
    model = _read("sphere.cad.json")
    fan = mio.mesh_cad_surface(model, 0.2, min_angle=25, return_report=True)[1]
    assert fan["min_angle_deg"] < 25 and fan["num_faces_below_bound"] > 0
    out, ringed = mio.mesh_cad_surface(
        model, 0.2, min_angle=25, pole_rings=2, return_report=True
    )
    _watertight(ringed)
    assert ringed["min_angle_deg"] >= 25 - 1e-9
    assert ringed["num_faces_below_bound"] == 0


@needs_core
def test_a_varying_metric_meets_the_bound_at_the_poles():
    """Roadmap 6.4.1: near its poles the sphere's parameters stretch the
    surface far from each hemisphere's mean metric; meshed in the surface's
    own metric (Lawson's flips and encroachment measured in it) the pole fans
    meet 25 degrees without rings, and the mesh stays watertight."""
    model = _read("sphere.cad.json")
    mean = mio.mesh_cad_surface(model, 0.2, min_angle=25, return_report=True)[1]
    out, varying = mio.mesh_cad_surface(
        model, 0.2, min_angle=25, metric="varying", return_report=True
    )
    _watertight(varying)
    assert mean["min_angle_deg"] < 25 <= varying["min_angle_deg"] + 1e-9
    assert varying["num_faces_below_bound"] == 0
    assert varying["num_points"] <= 1.5 * mean["num_points"]
    with pytest.raises(ValueError, match="metric"):
        mio.mesh_cad_surface(model, 0.2, metric="riemann")


@needs_core
def test_given_sizes_refine_and_report():
    """Roadmap 6.3.1 / 6.3.8: a face given its size meshes at it, an edge given
    its size is sampled at it, and the report names exactly the edge a finer
    neighbour limits a face by. Sizes only refine."""
    cube = _read("unit_cube.cad.json")
    edges = cube.info["brep"]["edges"]
    ids = {f["id"]: k for k, f in enumerate(cube.info["brep"]["faces"])}
    faces_of = [sorted(ids[r["id"]] for r in e["topology"]) for e in edges]
    shared = next(k for k, fs in enumerate(faces_of) if fs == [1, 2])
    out, report = mio.mesh_cad_surface(
        cube, 0.4, face_sizes={1: 0.1, 2: 0.4}, return_report=True
    )
    _watertight(report)
    tri, p = out.cells[0].data, out.points
    face = out.cell_data["cad:face"][0]
    e = np.concatenate(
        [
            np.linalg.norm(p[tri[face == 1][:, a]] - p[tri[face == 1][:, b]], axis=1)
            for a, b in ((0, 1), (1, 2), (2, 0))
        ]
    )
    assert np.median(e) == pytest.approx(0.1, rel=0.15)
    assert [f["limited_by_edges"] for f in report["faces"]] == [
        [shared] if k == 2 else [] for k in range(6)
    ]
    assert any(f"limited by edge {shared}" in w for w in report["warnings"])
    # The same as pairs; the edge size.
    again = mio.mesh_cad_surface(cube, 0.4, face_sizes=[(2, 0.4), (1, 0.1)])
    np.testing.assert_array_equal(again.points, out.points)
    out, report = mio.mesh_cad_surface(
        cube, 0.4, edge_sizes={shared: 0.05}, return_report=True
    )
    lines = [c for c in out.cells if c.type == "line"][0].data
    on = lines[out.cell_data["cad:edge"][1] == shared]
    gaps = out.points[on[:, 0]] - out.points[on[:, 1]]
    assert np.linalg.norm(gaps, axis=1).max() <= 0.05 + 1e-12
    # Empty sources are the plain mesh.
    plain = mio.mesh_cad_surface(cube, 0.4)
    same = mio.mesh_cad_surface(cube, 0.4, face_sizes={}, edge_sizes=[])
    np.testing.assert_array_equal(plain.points, same.points)


@needs_core
def test_a_background_mesh_sizes_the_model():
    cube = _read("unit_cube.cad.json")
    bg = mio.rectangle(dims=(20, 20), size=(1.0, 1.0), cell_type="triangle")
    bg.point_data["size"] = 0.05 + 0.3 * bg.points[:, 0]
    out, report = mio.mesh_cad_surface(cube, 0.4, background=bg, return_report=True)
    _watertight(report)
    tri, p = out.cells[0].data, out.points
    bottom = np.all(np.abs(p[tri][:, :, 2]) < 1e-9, axis=1)
    c = p[tri].mean(axis=1)
    for lo in (0.0, 0.25, 0.5, 0.75):
        sel = tri[bottom & (c[:, 0] >= lo) & (c[:, 0] < lo + 0.25)]
        e = np.concatenate(
            [
                np.linalg.norm(p[sel[:, a]] - p[sel[:, b]], axis=1)
                for a, b in ((0, 1), (1, 2), (2, 0))
            ]
        )
        target = 0.05 + 0.3 * (lo + 0.125)
        assert np.median(e) == pytest.approx(target, rel=0.25)


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"face_sizes": {6: 0.1}}, "face_sizes"),
        ({"face_sizes": {0: 0.0}}, "face_sizes"),
        ({"face_sizes": [(0.5, 0.1)]}, "face_sizes"),
        ({"edge_sizes": {-1: 0.1}}, "edge_sizes"),
        ({"edge_sizes": [0.1]}, "edge_sizes"),
        ({"metric": "anisotropic"}, "metric"),
    ],
)
def test_given_size_refusals(kwargs, match):
    with pytest.raises(ValueError, match=match):
        mio.mesh_cad_surface(_read("unit_cube.cad.json"), 0.4, **kwargs)


def test_a_background_without_its_array_is_refused():
    bg = mio.rectangle(dims=(2, 2), size=(1.0, 1.0), cell_type="triangle")
    with pytest.raises(ValueError, match="background"):
        mio.mesh_cad_surface(_read("unit_cube.cad.json"), 0.4, background=bg)


def _sphere_patch(radius=1.0):
    """One untrimmed rational patch: a full turn in u (four quarter arcs), a
    half meridian in v (two quarter arcs), so its loop is a seam twice and
    two poles."""
    brep = copy.deepcopy(_read("cylindrical_roof.cad.json").info["brep"])
    s = np.sqrt(0.5)
    ring = [(1, 0, 1), (1, 1, s), (0, 1, 1), (-1, 1, s), (-1, 0, 1)]
    ring += [(-1, -1, s), (0, -1, 1), (1, -1, s), (1, 0, 1)]
    meridian = [(0, -1, 1), (1, -1, s), (1, 0, 1), (1, 1, s), (0, 1, 1)]
    cps, weights = [], []
    for rr, z, wv in meridian:
        for x, y, wu in ring:
            cps.append([radius * rr * x, radius * rr * y, radius * z])
            weights.append(wu * wv)
    patch = brep["patches"][0]
    patch["degrees"] = [2, 2]
    patch["knots"] = [
        np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4], float),
        np.array([0, 0, 0, 1, 1, 2, 2, 2], float),
    ]
    patch["control_points"] = np.array(cps)
    patch["weights"] = np.array(weights)
    brep["faces"][0]["control_point_ids"] = np.arange(1, len(cps) + 1)
    return brep


@needs_core
@pytest.mark.parametrize("min_angle", [0, 25])
def test_a_surface_of_revolution_smaller_than_h_is_meshed(min_angle):
    """Roadmap 6.5.7: a face whose loop is a seam and two poles, smaller than
    h, has no lattice point, and every triangle of its boundary alone touches
    a pole or the seam twice (118 faces of the corpus). It is meshed again
    with an interior grid: watertight, outward, its volume positive."""
    out, report = mio.mesh_cad_surface(
        _sphere_patch(), 10.0, min_angle=min_angle, return_report=True
    )
    _watertight(report)
    area, volume = _measure(out)
    assert 0.5 * 4 * math.pi < area < 4 * math.pi
    assert 0.2 * 4 * math.pi / 3 < volume < 4 * math.pi / 3
    fine = mio.mesh_cad_surface(_sphere_patch(), 0.3)
    assert len(fine.cells[0].data) > len(out.cells[0].data)


@needs_core
@pytest.mark.parametrize("min_angle", [0, 25])
def test_virtual_edges_skip_a_sliver_of_an_edge(min_angle):
    """Roadmap 6.5.2: a rim split by an edge of h / 50 forces a point h / 50
    from its neighbour on both plates. With ``virtual_edges`` one corner of
    the sliver is suppressed: the plates still conform across the rim, the
    chain is reported, and the triangles near it are not slivers."""
    model = _read("split_rim.cad.json")
    off, r_off = mio.mesh_cad_surface(
        model, 0.25, min_angle=min_angle, return_report=True
    )
    on, r_on = mio.mesh_cad_surface(
        model, 0.25, min_angle=min_angle, virtual_edges=True, return_report=True
    )
    assert r_off["virtual_edges"] == [] and r_on["virtual_edges"] == [[3, 5]]
    assert r_on["num_failed"] == 0
    assert r_on["num_open_edges"] == r_off["num_open_edges"]
    assert np.isclose(_measure(on)[0], _measure(off)[0], rtol=0, atol=1e-12)
    if min_angle == 0:
        assert r_off["min_angle_deg"] < 5 < 20 < r_on["min_angle_deg"]
    else:
        assert r_on["min_angle_deg"] >= 25 - 1e-9
        assert len(on.cells[0].data) < 0.7 * len(off.cells[0].data)
    # The suppressed corner stays a point, used by no cell (as subdivide's
    # orphans do): clean(remove_orphans=True) drops it.
    used = np.unique(np.concatenate([c.data.ravel() for c in on.cells]))
    assert len(on.points) - len(used) == 1
    with pytest.raises(ValueError, match="virtual_edges"):
        mio.mesh_cad_surface(model, 0.25, order=2, virtual_edges=True)


def test_smoothing_moves_the_inside_only():
    model = _read("sphere.cad.json")
    plain, r0 = mio.mesh_cad_surface(model, 0.2, return_report=True)
    for method in ("laplacian", "odt"):
        out, report = mio.mesh_cad_surface(
            model, 0.2, smooth_iterations=2, smooth_method=method, return_report=True
        )
        _watertight(report)
        assert report["num_smoothed"] > 0
        boundary = r0["num_points"] - sum(f["num_interior_points"] for f in r0["faces"])
        np.testing.assert_array_equal(out.points[:boundary], plain.points[:boundary])
        np.testing.assert_allclose(np.linalg.norm(out.points, axis=1), 1.0, atol=1e-12)
    with pytest.raises(ValueError, match="smooth_method"):
        mio.mesh_cad_surface(model, 0.2, smooth_method="none")
    with pytest.raises(ValueError, match="min_angle"):
        mio.mesh_cad_surface(model, 0.2, min_angle=60.0)


@needs_core
def test_triangulate_polygon_2d_refines_to_a_min_angle():
    outer = [[0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]]
    m = mio.triangulate_polygon_2d(outer, target_size=0.3, min_angle=30.0)
    assert _measure(m)[0] == pytest.approx(3.0, rel=1e-12)
    assert _min_angles(m).min() >= 30.0 - 1e-9
    plain = mio.triangulate_polygon_2d(outer, target_size=0.3)
    np.testing.assert_array_equal(
        plain.points,
        mio.triangulate_polygon_2d(outer, target_size=0.3, min_angle=0).points,
    )
    with pytest.raises(ValueError, match="min_angle"):
        mio.triangulate_polygon_2d(outer, min_angle=75.0)


@needs_core
def test_cli_verb(tmp_path, capsys):
    out = tmp_path / "cube.vtu"
    code = mio._cli.main(
        [
            "mesh-cad-surface",
            str(CAD / "unit_cube.cad.json"),
            str(out),
            "--size",
            "0.25",
            "--json",
        ]
    )
    assert code == 0
    report = json.loads(capsys.readouterr().out)
    assert report["num_faces"] == 6 and report["num_open_edges"] == 0
    assert out.exists()
    code = mio._cli.main(
        [
            "mesh-cad-surface",
            str(CAD / "sphere.cad.json"),
            str(out),
            "--size",
            "0.2",
            "--min-angle",
            "25",
            "--curvature-deviation",
            "0.01",
            "--cells-across-gaps",
            "2",
            "--grading",
            "1.3",
            "--pole-rings",
            "2",
            "--json",
        ]
    )
    assert code == 0
    report = json.loads(capsys.readouterr().out)
    assert report["num_failed"] == 0 and report["min_angle_deg"] >= 25 - 1e-9
    code = mio._cli.main(
        [
            "mesh-cad-surface",
            str(CAD / "unit_cube.cad.json"),
            str(out),
            "--size",
            "0.4",
            "--face-sizes",
            "1=0.1,2=0.4",
            "--edge-sizes",
            "0=0.05",
            "--metric",
            "mean",
            "--json",
        ]
    )
    assert code == 0
    report = json.loads(capsys.readouterr().out)
    assert report["num_failed"] == 0 and report["virtual_edges"] == []
    assert any(f["limited_by_edges"] for f in report["faces"])
    code = mio._cli.main(
        [
            "mesh-cad-surface",
            str(CAD / "split_rim.cad.json"),
            str(out),
            "--size",
            "0.25",
            "--virtual-edges",
            "--json",
        ]
    )
    assert code == 0
    assert json.loads(capsys.readouterr().out)["virtual_edges"] == [[3, 5]]
    # A metric from a Hessian (roadmap 6.10.2), read from a file.
    bg_path = _layer_hessian(tmp_path)
    code = mio._cli.main(
        [
            "mesh-cad-surface",
            str(CAD / "sphere.cad.json"),
            str(out),
            "--size",
            "0.3",
            "--metric-background",
            str(bg_path),
            "--metric-array",
            "f:hessian",
            "--metric-error",
            "0.01",
            "--json",
        ]
    )
    assert code == 0
    assert json.loads(capsys.readouterr().out)["num_failed"] == 0


def _layer_hessian(tmp_path):
    """A coarse sphere carrying the Hessian of a layer at z = 0.3, written."""
    bg = mio.mesh_cad_surface(_read("sphere.cad.json"), 0.2, edges=False)
    x = np.asarray(bg.points)
    bg.point_data["f"] = np.tanh((x[:, 2] - 0.3) / 0.1)
    path = tmp_path / "hessian.vtu"
    # Uncompressed, so the core reads it in builds without zlib (the Windows wheel).
    mio.write(path, mio.hessian(bg, "f", location="point"), compression=None)
    return path


@needs_core
def test_pipeline_takes_a_metric_like_the_operation(tmp_path):
    """The metric keys (roadmap 6.10.2) reach mesh_cad_surface in both
    engines."""
    bg_path = _layer_hessian(tmp_path)
    doc = {
        "Version": 2,
        "Input": {"Path": str(CAD / "sphere.cad.json")},
        "Operations": [
            {
                "Op": "MeshCadSurface",
                "TargetSize": 0.3,
                "MetricBackground": str(bg_path),
                "MetricArray": "f:hessian",
                "MetricError": 0.01,
            }
        ],
    }
    direct = mio.mesh_cad_surface(
        _read("sphere.cad.json"),
        0.3,
        metric_background=mio.read(bg_path),
        metric_array="f:hessian",
        metric_error=0.01,
    )
    for report, mesh in _run_both(doc, tmp_path):
        np.testing.assert_array_equal(mesh.points, direct.points)
        assert report["steps"][0]["NumFailed"] == 0.0


def _run_both(doc, tmp_path):
    """The pipeline document through the Python engine and, with the core,
    the C++ one: their reports and their meshes."""
    out = []
    engines = [("py", mio.run_pipeline)]
    if hasattr(_core, "run_pipeline_json"):
        engines.append(("cpp", lambda d: _core.run_pipeline_json(json.dumps(d))))
    for name, run in engines:
        path = tmp_path / f"{name}.vtu"
        report = run(dict(doc, Output={"Path": str(path)}))
        out.append((report, mio.read(path)))
    return out


@needs_core
@pytest.mark.parametrize("version", [1, 2])
def test_pipeline_meshes_the_model_its_input_reads(tmp_path, version):
    """Roadmap 6.1.1: the B-rep a CAD read returns travels from Input to the
    MeshCadSurface step; both engines mesh it as mesh_cad_surface does."""
    doc = {
        "Version": version,
        "Input": {"Path": str(CAD / "unit_cube.cad.json")},
        "Operations": [{"Op": "MeshCadSurface", "TargetSize": 0.25, "MinAngle": 25}],
    }
    direct = mio.mesh_cad_surface(_read("unit_cube.cad.json"), 0.25, min_angle=25)
    results = _run_both(doc, tmp_path)
    for report, mesh in results:
        np.testing.assert_array_equal(mesh.points, direct.points)
        np.testing.assert_array_equal(mesh.cells[0].data, direct.cells[0].data)
        assert report["steps"][0]["op"] == "MeshCadSurface"
        assert report["steps"][0]["NumFaces"] == 6.0
        assert report["steps"][0]["NumOpenEdges"] == 0.0
    assert all(r == results[0][0] for r, _ in results)


@needs_core
def test_pipeline_sizes_and_rings_like_the_operation(tmp_path):
    """The sizing and pole keys reach mesh_cad_surface in both engines."""
    doc = {
        "Version": 2,
        "Input": {"Path": str(CAD / "sphere.cad.json")},
        "Operations": [
            {
                "Op": "MeshCadSurface",
                "TargetSize": 0.3,
                "MinAngle": 25,
                "CurvatureDeviation": 0.01,
                "Grading": 1.3,
                "CellsAcrossGaps": 2,
                "PoleRings": 2,
            }
        ],
    }
    direct = mio.mesh_cad_surface(
        _read("sphere.cad.json"),
        0.3,
        min_angle=25,
        curvature_deviation=0.01,
        grading=1.3,
        cells_across_gaps=2,
        pole_rings=2,
    )
    results = _run_both(doc, tmp_path)
    for _, mesh in results:
        np.testing.assert_array_equal(mesh.points, direct.points)
        np.testing.assert_array_equal(mesh.cells[0].data, direct.cells[0].data)
    assert all(r == results[0][0] for r, _ in results)


@needs_core
def test_pipeline_takes_given_sizes_like_the_operation(tmp_path):
    """FaceSizes / EdgeSizes (flat index, size lists), Background and Metric
    reach mesh_cad_surface in both engines."""
    bg = mio.rectangle(dims=(10, 10), size=(1.0, 1.0), cell_type="triangle")
    bg.point_data["size"] = 0.1 + 0.2 * bg.points[:, 0]
    bg_path = tmp_path / "bg.vtu"
    mio.write(bg_path, bg, compression=None)
    doc = {
        "Version": 2,
        "Input": {"Path": str(CAD / "unit_cube.cad.json")},
        "Operations": [
            {
                "Op": "MeshCadSurface",
                "TargetSize": 0.4,
                "FaceSizes": [1, 0.1],
                "EdgeSizes": [0, 0.05],
                "Background": str(bg_path),
                "Metric": "mean",
            }
        ],
    }
    direct = mio.mesh_cad_surface(
        _read("unit_cube.cad.json"),
        0.4,
        face_sizes={1: 0.1},
        edge_sizes={0: 0.05},
        background=mio.read(bg_path),
    )
    results = _run_both(doc, tmp_path)
    for _, mesh in results:
        np.testing.assert_array_equal(mesh.points, direct.points)
    assert all(r == results[0][0] for r, _ in results)
    doc = {
        "Version": 2,
        "Input": {"Path": str(CAD / "split_rim.cad.json")},
        "Operations": [
            {"Op": "MeshCadSurface", "TargetSize": 0.25, "VirtualEdges": True}
        ],
    }
    direct = mio.mesh_cad_surface(_read("split_rim.cad.json"), 0.25, virtual_edges=True)
    for _, mesh in _run_both(doc, tmp_path):
        np.testing.assert_array_equal(mesh.points, direct.points)
        np.testing.assert_array_equal(mesh.cells[0].data, direct.cells[0].data)


@needs_core
def test_pipeline_refines_then_meshes_the_refined_model(tmp_path):
    doc = {
        "Version": 2,
        "Input": {"Path": str(CAD / "unit_cube.cad.json")},
        "Operations": [
            {"Op": "IgaRefine", "Subdivisions": [2, 3]},
            {"Op": "MeshCadSurface", "TargetSize": 0.3, "Order": 2},
        ],
    }
    results = _run_both(doc, tmp_path)
    for report, mesh in results:
        assert [s["op"] for s in report["steps"]] == ["IgaRefine", "MeshCadSurface"]
        assert mesh.cells[0].type == "triangle6"
    assert all(r == results[0][0] for r, _ in results)


@pytest.mark.parametrize(
    "ops, message",
    [
        (
            [{"Op": "Quality"}, {"Op": "MeshCadSurface", "TargetSize": 0.25}],
            "step 1 'Quality' changed the mesh, which drops the model",
        ),
        ([{"Op": "MeshCadSurface"}], "'TargetSize' is required"),
    ],
)
def test_pipeline_model_steps_fail_by_name(tmp_path, ops, message):
    doc = {
        "Input": {"Path": str(CAD / "unit_cube.cad.json")},
        "Output": {"Path": str(tmp_path / "out.vtu")},
        "Operations": ops,
    }
    with pytest.raises(ValueError, match=message):
        mio.run_pipeline(doc)
    if hasattr(_core, "run_pipeline_json"):
        with pytest.raises(ValueError, match=message):
            _core.run_pipeline_json(json.dumps(doc))


@needs_core
def test_pipeline_meshes_the_solid_fit_cad_fits(tmp_path):
    """FitCad hands its fitted B-rep on (outside the native legacy freeform
    mode, which fits none), so a skin is fitted and re-meshed in one run."""
    src = tmp_path / "box.vtu"
    mio.write(src, mio.extract_surface(mio.box()), compression=None)
    fit = {"Op": "FitCad", "Backend": "native", "FreeformMode": "adaptive"}
    doc = {
        "Input": {"Path": str(src)},
        "Operations": [fit, {"Op": "MeshCadSurface", "TargetSize": 0.25}],
    }
    results = _run_both(doc, tmp_path)
    for report, mesh in results:
        assert report["steps"][1]["NumFaces"] == 6.0
        assert report["steps"][1]["NumOpenEdges"] == 0.0
    # FitCad's own report holds NaN counters (no volume natively), which never
    # compare equal; the meshing step's must agree.
    assert all(r["steps"][1] == results[0][0]["steps"][1] for r, _ in results)
    doc["Operations"][0] = dict(fit, FreeformMode="legacy")
    doc["Output"] = {"Path": str(tmp_path / "legacy.vtu")}
    with pytest.raises(ValueError, match="'FitCad' gave no B-rep"):
        mio.run_pipeline(doc)


def test_pipeline_model_step_needs_a_cad_input(tmp_path):
    src = tmp_path / "plain.vtu"
    mio.write(src, mio.extract_surface(mio.box()), compression=None)
    with pytest.raises(ValueError, match="Input format 'vtu' carries none"):
        mio.run_pipeline(
            {
                "Input": {"Path": str(src)},
                "Output": {"Path": str(tmp_path / "out.vtu")},
                "Operations": [{"Op": "IgaSample"}],
            }
        )


@needs_core
@pytest.mark.skipif(
    not getattr(_core, "__has_occt__", False), reason="needs OpenCASCADE"
)
def test_step_solid_through_opencascade(tmp_path):
    # A box skin written as a STEP solid, read back as a B-rep by
    # OpenCASCADE: watertight, area 6 and volume 1, against BRepMesh's own
    # triangulation (the read's mesh).
    path = tmp_path / "box.step"
    mio.write(path, mio.extract_surface(mio.box()))
    model = mio.read(path)
    assert len(model.info["brep"]["faces"]) == 6
    out, report = mio.mesh_cad_surface(model, 0.1, return_report=True)
    _watertight(report)
    a, v = _measure(out)
    assert a == pytest.approx(6.0, rel=1e-9)
    assert v == pytest.approx(1.0, rel=1e-9)
    assert a == pytest.approx(_measure(model)[0], rel=1e-9)


@needs_core
@pytest.mark.skipif(
    not getattr(_core, "__has_occt__", False), reason="needs OpenCASCADE"
)
def test_a_step_read_declares_its_tolerances(tmp_path):
    """Roadmap 6.5.1: an OpenCASCADE read writes each edge's tolerance and its
    vertices' into the edge's extra JSON, which model_tolerances reads; a
    sound box declares nothing the global tolerance would refuse."""
    path = tmp_path / "box.step"
    mio.write(path, mio.extract_surface(mio.box()))
    model = mio.read(path)
    for e in model.info["brep"]["edges"]:
        extra = json.loads(e["extra"])
        assert extra["tolerance"] >= 0.0
        assert len(extra["end_tolerances"]) == 2
    # Its eight vertices, each starting a trim of each of its three faces.
    vertices = model.info["brep"]["vertices"]
    assert len(vertices) == 8
    for v in vertices:
        assert json.loads(v["extra"])["tolerance"] >= 0.0
        assert len(v["topology"]) == 3
    out, report = mio.mesh_cad_surface(
        model, 0.1, model_tolerances=True, return_report=True
    )
    _watertight(report)
    assert report["tolerance_uses"] == []


def _collapsed_patch(path):
    """The planar bilinear patch with corners (0, 0), (1, 0), (0, 2) and
    (0, 2): its side v = 1 collapses, so it is the triangle with a 26.6-degree
    corner there (roadmap 6.4.4's probe)."""
    doc = {
        "breps": [
            {
                "brep_id": 1,
                "faces": [
                    {
                        "brep_id": 2,
                        "swapped_surface_normal": False,
                        "surface": {
                            "is_trimmed": False,
                            "is_rational": False,
                            "degrees": [1, 1],
                            "knot_vectors": [[0, 1], [0, 1]],
                            "control_points": [
                                [1, [0, 0, 0, 1]],
                                [2, [1, 0, 0, 1]],
                                [3, [0, 2, 0, 1]],
                                [4, [0, 2, 0, 1]],
                            ],
                        },
                        "boundary_loops": [],
                        "embedded_edges": [],
                    }
                ],
                "edges": [],
                "vertices": [],
            }
        ]
    }
    path.write_text(json.dumps(doc), encoding="utf-8")
    return mio.read(path, file_format="cad_json")


@needs_core
def test_a_collapsed_patch_is_flattened_and_meets_the_bound(tmp_path):
    """Roadmap 6.4.4: in its own parameters the collapsed corner is a pole
    whose fan neither metric rounds; meshed in the flattening of its coarse
    mesh (the conformal map), every triangle meets 25 degrees, and the face
    is reported with its condition numbers."""
    model = _collapsed_patch(tmp_path / "collapsed.cad.json")
    _, mean = mio.mesh_cad_surface(model, 0.1, min_angle=25, return_report=True)
    _, varying = mio.mesh_cad_surface(
        model, 0.1, min_angle=25, metric="varying", return_report=True
    )
    out, flat = mio.mesh_cad_surface(
        model, 0.1, min_angle=25, reparam_condition_number=4, return_report=True
    )
    assert mean["min_angle_deg"] < 25 and varying["min_angle_deg"] < 25
    assert flat["num_failed"] == 0 and flat["warnings"] == []
    assert flat["min_angle_deg"] >= 25 - 1e-9
    assert flat["num_points"] < mean["num_points"]
    (rec,) = flat["reparametrised"]
    assert rec["taken"] and rec["method"] == "conformal" and rec["face"] == 0
    assert rec["condition_number"] > 4 > rec["condition_number_disc"]
    # Every point on the patch, which is the plane z = 0.
    assert np.all(out.points[:, 2] == 0.0)
    a, _ = _measure(out)
    assert a == pytest.approx(1.0, rel=1e-12)
    # Below the threshold nothing changes; the refusals.
    plain = mio.mesh_cad_surface(model, 0.1, min_angle=25)
    high = mio.mesh_cad_surface(
        model, 0.1, min_angle=25, reparam_condition_number=1e300
    )
    assert np.array_equal(plain.points, high.points)
    for bad in (-1.0, 0.5, 1.0, float("nan")):
        with pytest.raises(ValueError, match="reparam_condition_number"):
            mio.mesh_cad_surface(model, 0.1, reparam_condition_number=bad)
    # Roadmap 6.4.4: curved on the flattened patch too, none straightened.
    _, curved = mio.mesh_cad_surface(
        model, 0.1, reparam_condition_number=4, order=2, return_report=True
    )
    assert curved["num_failed"] == 0 and curved["num_curved_straightened"] == 0


@needs_core
def test_the_sphere_is_flattened_through_its_poles():
    """Each hemisphere is flattened, its pole a point of the seam's cut: the
    sphere meets 25 degrees without pole rings, watertight, every point on
    it."""
    model = _read("sphere.cad.json")
    out, report = mio.mesh_cad_surface(
        model, 0.2, min_angle=25, reparam_condition_number=4, return_report=True
    )
    _watertight(report)
    assert report["min_angle_deg"] >= 25 - 1e-9
    assert [r["taken"] for r in report["reparametrised"]] == [True, True]
    np.testing.assert_allclose(np.linalg.norm(out.points, axis=1), 1.0, atol=1e-12)


@needs_core
def test_virtual_faces_mesh_a_sliver_with_its_neighbour():
    """Roadmap 6.5.2: three plates in a row, the middle one 0.01 wide. Alone,
    the sliver is a strip of needles; meshed with the neighbour it meets
    along the most edge, the pair meets 25 degrees with a fraction of the
    triangles, each triangle keeps its own face in cad:face, the shared edge
    leaves the line block and the merge is reported."""
    model = _read("sliver_strip.cad.json")
    _, alone = mio.mesh_cad_surface(model, 0.2, min_angle=25, return_report=True)
    out, merged = mio.mesh_cad_surface(
        model, 0.2, min_angle=25, virtual_faces=True, return_report=True
    )
    assert merged["num_failed"] == 0 and merged["num_misoriented_edges"] == 0
    assert merged["min_angle_deg"] >= 25 - 1e-9
    assert 3 * merged["num_triangles"] < alone["num_triangles"]
    (rec,) = merged["merged_faces"]
    assert rec["sliver"] == 1 and rec["merged"] and rec["host"] in (0, 2)
    assert len(rec["edges"]) == 1
    tri = [k for k, c in enumerate(out.cells) if c.type == "triangle"][0]
    faces = out.cell_data["cad:face"][tri]
    assert set(np.unique(faces)) == {0, 1, 2}
    line = [k for k, c in enumerate(out.cells) if c.type == "line"][0]
    assert rec["edges"][0] not in set(out.cell_data["cad:edge"][line])
    a, _ = _measure(out)
    assert a == pytest.approx(2.01, rel=1e-12)
    # Roadmap 6.5.2: curved on the merged face too.
    curved = mio.mesh_cad_surface(model, 0.2, virtual_faces=True, order=2)
    assert "triangle6" in {c.type for c in curved.cells}


def _planar_lens(path, length=4.0, theta_deg=2.5):
    """A planar face on [0, 5]^2 bounded by two circular trims (rational
    quadratic arcs) from (0.5, 1) to (0.5 + length, 1), each leaving the
    chord at ``theta_deg``: both corners are twice that."""
    th = math.radians(theta_deg)
    w, dy, x0, y0 = math.cos(th), 0.5 * length * math.tan(th), 0.5, 1.0

    def arc(index, ax, bx, my):
        return {
            "trim_index": index,
            "curve_direction": True,
            "parameter_curve": {
                "is_rational": True,
                "degree": 2,
                "knot_vector": [0, 0, 0, 1, 1, 1],
                "active_range": [0, 1],
                "control_points": [
                    [[ax, y0, 0, 1]],
                    [[x0 + 0.5 * length, my, 0, w]],
                    [[bx, y0, 0, 1]],
                ],
            },
        }

    face = {
        "brep_id": 2,
        "swapped_surface_normal": False,
        "surface": {
            "is_trimmed": True,
            "is_rational": False,
            "degrees": [1, 1],
            "knot_vectors": [[0, 5], [0, 5]],
            "control_points": [
                [1, [0, 0, 0, 1]],
                [2, [5, 0, 0, 1]],
                [3, [0, 5, 0, 1]],
                [4, [5, 5, 0, 1]],
            ],
        },
        "boundary_loops": [
            {
                "loop_type": "outer",
                "trimming_curves": [
                    arc(0, x0, x0 + length, y0 - dy),
                    arc(1, x0 + length, x0, y0 + dy),
                ],
            }
        ],
        "embedded_edges": [],
    }
    doc = {"breps": [{"brep_id": 1, "faces": [face], "edges": [], "vertices": []}]}
    path.write_text(json.dumps(doc))
    return mio.read(path)


@needs_core
def test_corner_lopping_states_its_bound_on_a_lens(tmp_path):
    """Roadmap 6.2.3, Pav and Walkington (2005): a planar face with two
    5-degree corners between circular trims, refined at 20 degrees with its
    corners lopped, ends under the stated bound with no budget warning, and
    every triangle under 20 degrees lies within twice its corner's radius
    (Lemma 7). Off, the mesh is the concentric shells'."""
    model = _planar_lens(tmp_path / "lens.cad.json")
    out, report = mio.mesh_cad_surface(
        model, 0.05, min_angle=20, lop_corners=True, return_report=True
    )
    assert report["num_failed"] == 0
    assert not any("did not finish" in w for w in report["warnings"])
    corners = report["lopped_corners"]
    assert len(corners) == 2
    for c in corners:
        assert c["face"] == 0 and 4.9 < c["angle_deg"] <= 5.0
        assert c["lfs0"] == pytest.approx(4.0, rel=1e-9)
        assert 0 < c["lop_radius"] <= (math.sqrt(2) - 1) * 4.0
        assert 20.0 < c["stated_bound_deg"] <= 20.705
        assert c["guaranteed"] and c["reason"] == "" and c["within_two_radii"]
    tri = [c for c in out.cells if c.type == "triangle"][0].data
    angles = _min_angles(out)
    centres = np.array([out.points[c["point"]] for c in corners])
    radii = np.array([c["lop_radius"] for c in corners])
    for t in np.nonzero(angles < 20 - 1e-9)[0]:
        d = np.linalg.norm(out.points[tri[t]][:, None, :] - centres[None], axis=2)
        assert np.any(np.all(d <= 2 * radii * (1 + 1e-9), axis=0))
    _, off = mio.mesh_cad_surface(model, 0.05, min_angle=20, return_report=True)
    assert off["lopped_corners"] == []


@needs_core
def test_model_tolerances_join_what_the_model_declares():
    """Roadmap 6.5.1: the capped cylinder whose top cap starts a quarter turn
    round meshes either way (the cap is rotated); a tolerance declared on
    every edge is read but accepts nothing the global one refuses, and a
    malformed one is warned."""
    model = _read("capped_cylinder_offset.cad.json")
    brep = copy.deepcopy(model.info["brep"])
    for e in brep["edges"]:
        e["extra"] = json.dumps({"tolerance": 1e-6})
    _, report = mio.mesh_cad_surface(
        brep, 0.3, model_tolerances=True, return_report=True
    )
    _watertight(report)
    assert report["tolerance_uses"] == []
    brep["edges"][0]["extra"] = json.dumps({"tolerance": -1})
    _, bad = mio.mesh_cad_surface(brep, 0.3, model_tolerances=True, return_report=True)
    assert any("declared tolerance" in w for w in bad["warnings"])


@needs_core
def test_the_models_tolerance_pairs_the_trims_no_edge_names():
    """Roadmap 6.5.1: the unit cube without edge topology, one face moved
    1e-5: its trims pair with their neighbours' only under the file's own
    ``tolerances.model_tolerance``, each pair reported with no edge."""
    model = _read("unit_cube.cad.json")
    brep = copy.deepcopy(model.info["brep"])
    brep["edges"] = []
    for b in brep["breps"]:
        b["edges"] = []
    patch = brep["patches"][brep["faces"][0]["patch"]]
    patch["control_points"] = np.asarray(patch["control_points"]) + 1e-5
    extra = json.loads(brep["extra"]) if brep["extra"] else {}
    extra["tolerances"] = {"model_tolerance": 1e-4}
    brep["extra"] = json.dumps(extra)
    _, off = mio.mesh_cad_surface(brep, 0.2, return_report=True)
    assert off["num_edges"] == 16
    _, on = mio.mesh_cad_surface(brep, 0.2, model_tolerances=True, return_report=True)
    _watertight(on)
    assert on["num_edges"] == 12
    uses = on["tolerance_uses"]
    assert len(uses) == 4
    assert all(u["edge"] == -1 and u["where"] == "pairing" for u in uses)
    assert all(u["declared"] == 1e-4 for u in uses)


def _stacked_plates(path, ys):
    """Untrimmed planar bilinear patches stacked in y, ``[0, 1] x [ys[k],
    ys[k + 1]]`` each, with no edge topology (their trims pair by image)."""
    faces = []
    cid = 1
    for k in range(len(ys) - 1):
        y0, y1 = ys[k], ys[k + 1]
        corners = [[0, y0, 0], [1, y0, 0], [0, y1, 0], [1, y1, 0]]
        cps = []
        for c in corners:
            cps.append([cid, [*c, 1]])
            cid += 1
        faces.append(
            {
                "brep_id": 100 + k,
                "swapped_surface_normal": False,
                "surface": {
                    "is_trimmed": False,
                    "is_rational": False,
                    "degrees": [1, 1],
                    "knot_vectors": [[0, 1], [0, 1]],
                    "control_points": cps,
                },
                "boundary_loops": [],
                "embedded_edges": [],
            }
        )
    doc = {"breps": [{"brep_id": 1, "faces": faces, "edges": [], "vertices": []}]}
    path.write_text(json.dumps(doc), encoding="utf-8")
    return mio.read(path, file_format="cad_json")


@needs_core
def test_virtual_faces_join_a_chain_of_slivers_into_one_host(tmp_path):
    """Roadmap 6.5.2: three strips h / 20 wide beside a plate. The first
    merges with the plate; the others border only slivers, and join the group
    it forms in the rounds that follow, each triangle keeping its own face."""
    h, w = 0.2, 0.01
    ys = [0.0, 1.0, 1.0 + w, 1.0 + 2 * w, 1.0 + 3 * w]
    model = _stacked_plates(tmp_path / "chain.cad.json", ys)
    _, alone = mio.mesh_cad_surface(model, h, min_angle=25, return_report=True)
    out, merged = mio.mesh_cad_surface(
        model, h, min_angle=25, virtual_faces=True, return_report=True
    )
    assert merged["num_failed"] == 0
    assert merged["min_angle_deg"] >= 25 - 1e-9
    assert merged["num_triangles"] < alone["num_triangles"]
    recs = merged["merged_faces"]
    assert [r["sliver"] for r in recs] == [1, 2, 3]
    assert all(r["merged"] and r["host"] == 0 for r in recs)
    faces = out.cell_data["cad:face"][0]
    assert sorted(set(int(f) for f in faces)) == [0, 1, 2, 3]
    a, _ = _measure(out)
    assert a == pytest.approx(1.0 + 3 * w, abs=1e-12)


@needs_core
def test_a_metric_from_the_hessian_stretches_along_the_layer():
    """Roadmap 6.10.2: the sphere sized by the metric of a layer's Hessian
    (``R |L| R^T / (8 eps)``): the triangles in the layer stretch along it,
    with far fewer points than a uniform mesh fine enough for the layer."""
    sphere = _read("sphere.cad.json")
    bg = mio.mesh_cad_surface(sphere, 0.1, edges=False, pole_rings=2, min_angle=25)
    x = np.asarray(bg.points)
    bg.point_data["f"] = np.tanh((x[:, 2] - 0.3) / 0.05)
    bg = mio.hessian(bg, "f", location="point")
    out, rep = mio.mesh_cad_surface(
        sphere,
        0.2,
        min_angle=25,
        pole_rings=2,
        edges=False,
        metric_background=bg,
        metric_array="f:hessian",
        metric_error=8e-3,
        return_report=True,
    )
    _, uniform = mio.mesh_cad_surface(
        sphere, 0.05, min_angle=25, pole_rings=2, edges=False, return_report=True
    )
    assert rep["num_failed"] == 0 and rep["num_open_edges"] == 0
    assert rep["num_points"] < 0.15 * uniform["num_points"]
    p = np.asarray(out.points)
    tri = [c for c in out.cells if c.type == "triangle"][0].data
    z = p[tri][:, :, 2].mean(axis=1)
    t = p[tri[np.abs(z - 0.3) < 0.05]]
    edge = np.linalg.norm(t - np.roll(t, 1, axis=1), axis=2)
    area2 = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1)
    assert np.mean(edge.max(axis=1) ** 2 / area2) > 3.0
    with pytest.raises(ValueError, match="metric_error"):
        mio.mesh_cad_surface(
            sphere, 0.2, metric_background=bg, metric_array="f:hessian"
        )
    bg.point_data["two"] = np.zeros((len(x), 2))
    with pytest.raises(ValueError, match="nine components"):
        mio.mesh_cad_surface(
            sphere, 0.2, metric_background=bg, metric_array="two", metric_error=1
        )
    with pytest.raises(ValueError, match="need a metric_background"):
        mio.mesh_cad_surface(sphere, 0.2, metric_error=1.0)


@needs_core
def test_a_hessian_sized_mesh_reaches_the_layer_with_fewer_points():
    """Roadmap 6.11: a background size from the layer's tangential Hessian
    (sqrt(8 eps / |lambda|max), Loseille and Alauzet 2011) reaches the
    uniform mesh's L2 interpolation error with well under half its points."""
    sphere = _read("sphere.cad.json")

    def layer(x):
        return np.tanh((x[..., 2] - 0.3) / 0.05)

    a1, b1 = 0.059715871789770, 0.470142064105115
    a2, b2 = 0.797426985353087, 0.101286507323456
    bary = np.array(
        [
            [1 / 3, 1 / 3, 1 / 3],
            [a1, b1, b1],
            [b1, a1, b1],
            [b1, b1, a1],
            [a2, b2, b2],
            [b2, a2, b2],
            [b2, b2, a2],
        ]
    )
    weights = np.array([0.225] + [0.132394152788506] * 3 + [0.125939180544827] * 3)

    def l2(mesh):
        p = np.asarray(mesh.points)
        t = [c.data for c in mesh.cells if c.type == "triangle"][0]
        c = p[t]
        area = 0.5 * np.linalg.norm(
            np.cross(c[:, 1] - c[:, 0], c[:, 2] - c[:, 0]), axis=1
        )
        q = np.einsum("qk,nkd->nqd", bary, c)
        on = q / np.linalg.norm(q, axis=2, keepdims=True)
        e = np.einsum("qk,nk->nq", bary, layer(p)[t]) - layer(on)
        return np.sqrt((area[:, None] * weights * e**2).sum())

    opts = dict(min_angle=25, pole_rings=2, edges=False)
    uniform = mio.mesh_cad_surface(sphere, 0.05, **opts)
    bg = mio.mesh_cad_surface(sphere, 0.1, **opts)
    bg.point_data["f"] = layer(np.asarray(bg.points))
    hess = mio.hessian(bg, "f", location="point").point_data["f:hessian"]
    H = np.asarray(hess).reshape(-1, 3, 3)
    H = 0.5 * (H + H.transpose(0, 2, 1))
    n = np.asarray(bg.points) / np.linalg.norm(bg.points, axis=1, keepdims=True)
    P = np.eye(3) - n[:, :, None] * n[:, None, :]
    lam = np.abs(np.linalg.eigvalsh(P @ H @ P)).max(axis=1)
    size = np.sqrt(8 * 8e-3 / np.maximum(lam, 1e-300))
    bg.point_data["size"] = np.clip(size, 0.05 / 4, 0.2)
    sized = mio.mesh_cad_surface(sphere, 0.2, background=bg, **opts)
    assert l2(sized) < l2(uniform)
    assert len(sized.points) < 0.4 * len(uniform.points)


@needs_core
def test_a_chain_of_slivers_alone_is_hosted_by_its_widest(tmp_path):
    """Roadmap 6.5.2: three strips narrower than h and nothing wider beside
    them. The widest (the middle one) hosts the chain; its own record names
    itself as the host."""
    h, w = 0.2, 0.01
    ys = [0.0, w, 3 * w, 4 * w]
    model = _stacked_plates(tmp_path / "strips.cad.json", ys)
    _, alone = mio.mesh_cad_surface(model, h, min_angle=25, return_report=True)
    out, merged = mio.mesh_cad_surface(
        model, h, min_angle=25, virtual_faces=True, return_report=True
    )
    assert merged["num_failed"] == 0
    assert merged["num_triangles"] < alone["num_triangles"]
    recs = merged["merged_faces"]
    assert [r["sliver"] for r in recs] == [0, 1, 2]
    assert all(r["merged"] and r["host"] == 1 for r in recs)
    a, _ = _measure(out)
    assert a == pytest.approx(4 * w, abs=1e-12)


@needs_core
def test_untangling_keeps_folded_elements_curved():
    """Roadmap 6.6.3: the plate whose hole nearly touches its rim, curved at
    h = 0.6. Without untangling some elements fold and are straightened; with
    it the plate's own nodes move first and none is."""
    model = _read("plate_with_near_hole.cad.json")
    for order in (2, 3):
        _, plain = mio.mesh_cad_surface(model, 0.6, order=order, return_report=True)
        out, report = mio.mesh_cad_surface(
            model, 0.6, order=order, untangle=True, return_report=True
        )
        assert plain["num_curved_straightened"] > 0
        assert report["num_curved_straightened"] == 0
        assert report["num_untangled"] == plain["num_curved_straightened"]
        assert out.cells[0].type == ("triangle6" if order == 2 else "triangle10")
        # Every node stays on the plane.
        assert np.allclose(np.asarray(out.points)[:, 2], 0.0)
    with pytest.raises(ValueError, match="untangle"):
        mio.mesh_cad_surface(model, 0.6, untangle=True)


@needs_core
def test_untangling_slides_nodes_along_the_edges():
    """Roadmap 6.6.3: at h = 1.75 the quadratic plate's two folds are set by
    the rim and the hole, which the plate's own nodes cannot undo; the nodes
    on those edges slide along them, so none is straightened and every edge
    node stays on its circle."""
    model = _read("plate_with_near_hole.cad.json")
    _, plain = mio.mesh_cad_surface(model, 1.75, order=2, return_report=True)
    out, report = mio.mesh_cad_surface(
        model, 1.75, order=2, untangle=True, return_report=True
    )
    assert plain["num_curved_straightened"] == 2
    assert report["num_curved_straightened"] == 0
    assert report["num_untangled"] == 2
    p = np.asarray(out.points)
    ids = np.unique(
        np.concatenate([c.data.ravel() for c in out.cells if c.type == "line3"])
    )
    q = p[ids, :2]
    rim = np.isclose(np.hypot(q[:, 0], q[:, 1]), 2.0, atol=1e-9)
    hole = q[~rim]
    # The hole's nodes lie on one circle: fit its centre and radius.
    a = np.c_[2 * hole, np.ones(len(hole))]
    c = np.linalg.lstsq(a, (hole**2).sum(axis=1), rcond=None)[0]
    r = np.sqrt(c[2] + c[0] ** 2 + c[1] ** 2)
    assert np.allclose(np.hypot(*(hole - c[:2]).T), r, atol=1e-9)
    again, _ = mio.mesh_cad_surface(
        model, 1.75, order=2, untangle=True, return_report=True
    )
    assert np.array_equal(np.asarray(again.points), p)


@needs_core
def test_quads_dominant_and_all(tmp_path):
    """Roadmap 6.8: the unit cube all-quad (no triangle left, convex quads,
    each keeping its face, written to Gmsh with a physical surface per face),
    and the sphere quad-dominant."""
    cube = _read("unit_cube.cad.json")
    out, report = mio.mesh_cad_surface(
        cube, 0.2, min_angle=25, quads="all", return_report=True
    )
    assert {c.type for c in out.cells} == {"quad", "line"}
    assert report["num_triangles"] == 0 and report["num_leftover_triangles"] == 0
    assert report["num_open_edges"] == 0 and report["num_misoriented_edges"] == 0
    assert report["min_quad_scaled_jacobian"] > 0.3
    q = mio.compute_quality(out)["metrics"]["quality:scaled_jacobian"]
    assert q["min"] == pytest.approx(report["min_quad_scaled_jacobian"])
    faces = out.cell_data["cad:face"][0]
    assert sorted(set(int(f) for f in faces)) == list(range(6))
    a = 0.0
    p = np.asarray(out.points)
    for quad in out.cells[0].data:
        a += 0.5 * np.linalg.norm(
            np.cross(p[quad[2]] - p[quad[0]], p[quad[3]] - p[quad[1]])
        )
    assert a == pytest.approx(6.0, abs=1e-12)
    path = tmp_path / "cube.msh"
    mio.write(path, out, file_format="gmsh")
    back = mio.read(path)
    assert (
        sum(len(c.data) for c in back.cells if c.type == "quad") == report["num_quads"]
    )
    sphere = _read("sphere.cad.json")
    _, plain = mio.mesh_cad_surface(sphere, 0.2, min_angle=25, return_report=True)
    dom, d = mio.mesh_cad_surface(
        sphere, 0.2, min_angle=25, quads="dominant", return_report=True
    )
    assert {c.type for c in dom.cells} == {"triangle", "quad", "line"}
    assert 2 * d["num_quads"] + d["num_leftover_triangles"] == plain["num_triangles"]
    assert d["num_open_edges"] == 0
    with pytest.raises(ValueError, match="quads"):
        mio.mesh_cad_surface(sphere, 0.2, quads="hex")
    # Roadmap 6.8.3: quad16, untangling with quads, and curved quad-dominant.
    cubic, c = mio.mesh_cad_surface(
        sphere, 0.3, quads="all", order=3, return_report=True
    )
    assert {b.type for b in cubic.cells} == {"quad16", "line4"}
    assert c["num_curved_straightened"] == 0
    _, u = mio.mesh_cad_surface(
        sphere, 0.3, quads="all", order=2, untangle=True, return_report=True
    )
    assert u["num_curved_straightened"] == 0
    curved_dom = mio.mesh_cad_surface(sphere, 0.3, quads="dominant", order=2)
    assert {"triangle6", "quad9"} <= {b.type for b in curved_dom.cells}
    # Roadmap 6.8.1: the matching (the default) leaves fewer triangles than
    # greedy pairing.
    _, greedy = mio.mesh_cad_surface(
        sphere,
        0.2,
        min_angle=25,
        quads="dominant",
        quad_pairs="greedy",
        return_report=True,
    )
    assert d["num_leftover_triangles"] < greedy["num_leftover_triangles"]
    with pytest.raises(ValueError, match="quad_pairs"):
        mio.mesh_cad_surface(sphere, 0.2, quads="dominant", quad_pairs="best")


@needs_core
def test_quads_all_curved_and_sized():
    """Roadmap 6.8.3: quads="all" with order 2 writes quad9 whose every node
    lies on the sphere, closed and none straightened; and quads="all" takes
    a size source, asking it for twice its size."""
    sphere = _read("sphere.cad.json")
    out, report = mio.mesh_cad_surface(
        sphere, 0.2, quads="all", order=2, return_report=True
    )
    assert {c.type for c in out.cells} == {"quad9", "line3"}
    assert report["num_open_edges"] == 0
    assert report["num_curved_straightened"] == 0
    p = np.asarray(out.points)
    assert np.allclose(np.linalg.norm(p, axis=1), 1.0, atol=1e-12)
    cyl = _read("capped_cylinder.cad.json")
    _, uniform = mio.mesh_cad_surface(cyl, 0.2, quads="all", return_report=True)
    sized, r = mio.mesh_cad_surface(
        cyl, 0.2, quads="all", curvature_deviation=0.002, return_report=True
    )
    assert {c.type for c in sized.cells} == {"quad", "line"}
    assert r["num_open_edges"] == 0
    assert r["num_quads"] > 2 * uniform["num_quads"]


@needs_core
def test_pipeline_meshes_quads_and_untangles(tmp_path):
    """Roadmap 6.8 and 6.6.3 through both engines: Quads and Untangle reach
    mesh_cad_surface, and the step reports NumQuads and NumUntangled."""
    doc = {
        "Version": 2,
        "Input": {"Path": str(CAD / "unit_cube.cad.json")},
        "Operations": [{"Op": "MeshCadSurface", "TargetSize": 0.25, "Quads": "all"}],
    }
    direct = mio.mesh_cad_surface(_read("unit_cube.cad.json"), 0.25, quads="all")
    results = _run_both(doc, tmp_path)
    for report, mesh in results:
        assert {c.type for c in mesh.cells} == {"quad", "line"}
        np.testing.assert_array_equal(mesh.points, direct.points)
        assert report["steps"][0]["NumQuads"] == float(len(direct.cells[0].data))
        assert report["steps"][0]["NumUntangled"] == 0.0
    assert all(r == results[0][0] for r, _ in results)
    doc["Operations"][0].update(Quads="none", Order=2, Untangle=True)
    for report, _ in _run_both(doc, tmp_path):
        assert report["steps"][0]["NumUntangled"] == 0.0
