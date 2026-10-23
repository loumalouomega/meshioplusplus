"""mesh_cad_surface: conforming surface meshes of CAD B-reps, and
triangulate_polygon_2d.

The CAD JSON fixtures through the Python wrapper: the closed solids
(``unit_cube``, ``capped_cylinder``, ``sphere``) come out watertight with
their area and enclosed volume converging at second order, the trimmed plate
keeps its hole, the Scordelis-Lo roof (an untrimmed patch whose parameter
scales differ fortyfold) and the two-patch beam (a coupling edge shared by
two faces, and an embedded edge meshed as mesh edges) mesh conformingly;
the report, the refusals before the core, strict mode, a bare B-rep dict,
round trips through VTU and Gmsh, both CLI verbs and the pipeline's refusal
by name. With OpenCASCADE, a STEP solid.
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
    np.testing.assert_array_equal(back.points, out.points)
    tri = np.concatenate([c.data for c in back.cells if c.type == "triangle"])
    np.testing.assert_array_equal(
        np.sort(tri, axis=0), np.sort(out.cells[0].data, axis=0)
    )


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


def test_pipeline_refuses_the_step_by_name():
    with pytest.raises(ValueError, match="mesh-cad-surface"):
        mio.run_pipeline(
            {
                "Input": {"Path": "a"},
                "Output": {"Path": "b"},
                "Operations": [{"Op": "MeshCadSurface"}],
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
