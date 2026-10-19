"""Tests for the ``step`` writer (roadmap §6.3.5.1): a closed surface as an
AP214 B-rep, one face per facet (``faces="facet"``), per coplanar group
(``"planar"``, the default without fitted patches) or per fitted patch
(``"fitted"``). The engines write the same faceted bytes; an OpenCASCADE check
(skipped without ``OCP``, from ``pip install cadquery-ocp``) reads each file
back as a valid solid of the right volume. The C++ behaviour is pinned in
``tests/cpp/test_step.cpp``."""

import re

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus.step import _step

from .analytic_surfaces import box_skin, cylinder, torus

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")


@pytest.fixture(autouse=True)
def _fixed_provenance(monkeypatch):
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "0")


def _two_boxes():
    b = box_skin(2)
    pts = np.vstack([b.points, b.points + 3.0])
    conn = np.vstack([b.cells[0].data, b.cells[0].data + len(b.points)])
    return mio.Mesh(pts, [("quad", conn)])


def _inward():
    b = box_skin(2)
    return mio.Mesh(b.points, [("quad", b.cells[0].data[:, ::-1].copy())])


def _tet():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1.0]])
    return mio.Mesh(pts, [("tetra", np.array([[0, 1, 2, 3]]))])


def _warped():
    b = box_skin(1)
    pts = b.points.copy()
    pts[0] += [0.0, 0.0, -0.1]
    return mio.Mesh(pts, b.cells)


def _polygons():
    b = box_skin(1)
    return mio.Mesh(b.points, [("polygon", [row for row in b.cells[0].data])])


# mesh, solids, faces with faces="facet", faces with "planar", volume
CASES = {
    "box": (box_skin(3), 1, 54, 6, 1.0),
    "torus": (torus(n=16, m=12), 1, 192, 192, None),
    "inward": (_inward(), 1, 24, 6, 1.0),
    "two": (_two_boxes(), 2, 48, 12, 2.0),
    "tet": (_tet(), 1, 4, 4, 1.0 / 6.0),
    "warped": (_warped(), 1, 7, 7, None),
    "polygons": (_polygons(), 1, 6, 6, 1.0),
}


def _count(text, entity):
    return len(re.findall(r"=" + entity + r"\(", text))


def test_the_format_reads_through_opencascade_only():
    formats = mio.formats()
    assert "step" in formats["writable"] and "step" in formats["readable"]
    assert "iges" in formats["readable"] and "iges" not in formats["writable"]
    assert formats["extensions"][".step"] == ["step"]
    assert formats["extensions"][".stp"] == ["step"]
    assert formats["extensions"][".igs"] == ["iges"]
    with pytest.raises(Exception):
        mio.read("nothing.step")


HAS_OCCT = bool(getattr(_core, "__has_occt__", False)) if _core is not None else False


@pytest.mark.skipif(HAS_OCCT, reason="checks the build without OpenCASCADE")
@pytest.mark.parametrize("fmt", ["step", "iges"])
def test_reading_without_opencascade_is_named(tmp_path, fmt):
    path = tmp_path / "box.step"
    _step.write(path, box_skin(2))
    with pytest.raises(NotImplementedError, match="MESHIOPLUSPLUS_WITH_OCCT"):
        mio.read(path, file_format=fmt)


@pytest.mark.skipif(
    not HAS_OCCT, reason="core built without -DMESHIOPLUSPLUS_WITH_OCCT=ON"
)
def test_opencascade_reads_back_the_faceted_solid(tmp_path):
    path = tmp_path / "box.step"
    mio.write(path, box_skin(2), faces="facet")
    mesh = mio.read(path)
    assert mesh.cells[0].type == "triangle"
    brep = mesh.info["brep"]
    assert len(brep["faces"]) == 24 and len(brep["patches"]) == 24
    assert all(f["trimmed"] and len(f["loops"]) == 1 for f in brep["faces"])
    assert np.allclose(mesh.points.min(axis=0), 0.0) and np.allclose(
        mesh.points.max(axis=0), 1.0
    )
    face_ids = {f["id"] for f in brep["faces"]}
    assert set(np.unique(mesh.cell_data["cad:id"][0]).tolist()) == face_ids
    # Every inner mesh edge of the box bounds two faces, every one is a 3-D curve.
    assert all(
        e["curve"] is not None and len(e["topology"]) == 2 for e in brep["edges"]
    )
    # A STEP read writes as Kratos CAD JSON, which samples back onto the box.
    mio.write(tmp_path / "box.cad.json", mesh)
    back = mio.read(tmp_path / "box.cad.json")
    assert len(back.info["brep"]["faces"]) == 24
    assert np.allclose(back.points.min(axis=0), 0.0, atol=1e-12)
    assert np.allclose(back.points.max(axis=0), 1.0, atol=1e-12)


@pytest.mark.skipif(
    not HAS_OCCT, reason="core built without -DMESHIOPLUSPLUS_WITH_OCCT=ON"
)
def test_deflection_controls_the_triangulation(tmp_path):
    # A curved face: the faceted writer's planar ones triangulate the same at
    # any deflection, so the torus goes through fit_cad as one toroidal face.
    path = tmp_path / "torus.step"
    mio.fit_cad(mio.fit_primitives(mio.segment(torus()), tolerance=1e-3), path)
    coarse = mio.step.read(path, deflection=0.2)
    fine = mio.step.read(path, deflection=0.01)
    assert len(fine.cells[0].data) > len(coarse.cells[0].data)
    with pytest.raises(ValueError):
        mio.step.read(path, angular_deflection=0.0)
    with pytest.raises(Exception):
        bad = tmp_path / "bad.step"
        bad.write_text("not a step file\n")
        mio.read(bad)


@pytest.mark.skipif(
    not HAS_OCCT, reason="core built without -DMESHIOPLUSPLUS_WITH_OCCT=ON"
)
def test_coplanar_facets_read_back_as_one_face(tmp_path):
    path = tmp_path / "box.step"
    mio.write(path, box_skin(2))
    brep = mio.read(path).info["brep"]
    assert len(brep["faces"]) == 6
    assert all(f["trimmed"] and len(f["loops"]) == 1 for f in brep["faces"])


@pytest.mark.parametrize("faces", ["facet", "planar"])
@pytest.mark.parametrize("name", sorted(CASES))
def test_topology_is_closed_and_counted(tmp_path, name, faces):
    if faces != "facet" and _core is None:
        pytest.skip("needs the compiled core")
    mesh, solids, n_facet, n_planar, _volume = CASES[name]
    path = tmp_path / "part.stp"
    mio.write(path, mesh, faces=faces)
    faces = n_facet if faces == "facet" else n_planar
    text = path.read_text()
    assert text.startswith("ISO-10303-21;\n") and text.endswith("END-ISO-10303-21;\n")
    assert _count(text, "MANIFOLD_SOLID_BREP") == solids
    assert _count(text, "CLOSED_SHELL") == solids
    assert _count(text, "ADVANCED_FACE") == faces
    assert _count(text, "PLANE") == faces
    # every edge curve is used by exactly two oriented edges, once each way
    edges = re.findall(r"#(\d+)=EDGE_CURVE\(", text)
    uses = re.findall(r"ORIENTED_EDGE\('',\*,\*,#(\d+),\.([TF])\.\)", text)
    for e in edges:
        senses = sorted(s for i, s in uses if i == e)
        assert senses == ["F", "T"], e


def test_a_cube_has_six_faces_twelve_edges_eight_vertices(tmp_path):
    path = tmp_path / "cube.step"
    mio.write(path, box_skin(1))
    text = path.read_text()
    assert (_count(text, "ADVANCED_FACE"), _count(text, "EDGE_CURVE")) == (6, 12)
    assert _count(text, "VERTEX_POINT") == 8
    assert "FILE_NAME('cube.step'," in text
    assert "SI_UNIT(.MILLI.,.METRE.)" in text


@needs_core
@pytest.mark.parametrize("name", sorted(CASES))
def test_the_core_and_the_reference_write_the_same_bytes(tmp_path, name):
    mesh = CASES[name][0]
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    _core.step_write(str(tmp_path / "a" / "x.step"), mesh, "facet")
    _step.write(str(tmp_path / "b" / "x.step"), mesh)
    assert (tmp_path / "a" / "x.step").read_bytes() == (
        tmp_path / "b" / "x.step"
    ).read_bytes()


def test_reals_carry_a_decimal_point():
    assert _step._real(3.0) == "3."
    assert _step._real(1e-5) == "1.0000000000000001E-05"
    assert _step._real(1e20) == "1.E+20"
    assert _step._real(-0.25) == "-0.25"
    assert _step._real(2.0**-20) == "9.5367431640625E-07"


@pytest.mark.parametrize("writer", ["core", "reference"])
def test_an_open_surface_is_refused_naming_its_edges(tmp_path, writer):
    if writer == "core" and _core is None:
        pytest.skip("needs the compiled core")
    write = (
        (lambda p, m: _core.step_write(str(p), m)) if writer == "core" else _step.write
    )
    with pytest.raises(Exception, match=r"64 open.*first \(0, 1\)"):
        write(tmp_path / "x.step", cylinder())


@needs_core
def test_both_engines_refuse_with_the_same_message(tmp_path):
    messages = []
    for write in (lambda p, m: _core.step_write(str(p), m), _step.write):
        with pytest.raises(Exception) as info:
            write(tmp_path / "x.step", cylinder())
        messages.append(str(info.value))
    assert messages[0] == messages[1]


def test_lines_are_dropped_and_no_facets_is_refused(tmp_path):
    b = box_skin(1)
    mesh = mio.Mesh(b.points, [b.cells[0], ("line", np.array([[0, 1]]))])
    mio.write(tmp_path / "a.step", mesh)
    line = mio.Mesh(b.points, [("line", np.array([[0, 1]]))])
    with pytest.raises(Exception, match="no surface facets"):
        mio.write(tmp_path / "b.step", line)


@pytest.mark.parametrize("faces", ["facet", "planar"])
@pytest.mark.parametrize("name", sorted(CASES))
def test_opencascade_reads_a_valid_solid(tmp_path, name, faces):
    if faces != "facet" and _core is None:
        pytest.skip("needs the compiled core")
    pytest.importorskip("OCP")
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    mesh, solids, _n_facet, _n_planar, volume = CASES[name]
    path = tmp_path / "x.step"
    mio.write(path, mesh, faces=faces)
    reader = STEPControl_Reader()
    assert reader.ReadFile(str(path)) == IFSelect_RetDone
    reader.TransferRoots()
    shape = reader.OneShape()
    assert BRepCheck_Analyzer(shape).IsValid()
    found, total = 0, 0.0
    explorer = TopExp_Explorer(shape, TopAbs_SOLID)
    while explorer.More():
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(explorer.Current(), props)
        assert props.Mass() > 0.0
        total += props.Mass()
        found += 1
        explorer.Next()
    assert found == solids
    if volume is not None:
        assert total == pytest.approx(volume, rel=1e-9)
