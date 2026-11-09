"""STEP curves and element solids (roadmap 7.6.5)."""

import math

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "step_write_elements"), reason="needs the C++ core"
)


def _fitted(surface):
    return mp.fit_primitives(mp.segment(surface, method="primitive"))


def _occ(path):
    """The shape OpenCASCADE reads from ``path`` and whether it is valid, or
    ``None`` without OpenCASCADE's Python bindings (``OCP``)."""
    try:
        import OCP  # noqa: F401
    except ImportError:
        return None
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Reader

    reader = STEPControl_Reader()
    assert reader.ReadFile(str(path)) == IFSelect_RetDone
    reader.TransferRoots()
    shape = reader.OneShape()
    return shape, BRepCheck_Analyzer(shape).IsValid()


def _solid_volumes(shape):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    out = []
    e = TopExp_Explorer(shape, TopAbs_SOLID)
    while e.More():
        g = GProp_GProps()
        BRepGProp.VolumeProperties_s(e.Current(), g)
        out.append(g.Mass())
        e.Next()
    return out


def test_a_fitted_cylinder_and_its_planes(tmp_path):
    fitted = _fitted(mp.extract_skin(mp.cylinder(2.0, 5.0, 12, 10)))
    path = tmp_path / "cyl.step"
    rep = mp.write_step_elements(path, fitted)
    assert (rep["num_cylinders"], rep["num_planes"], rep["num_skipped"]) == (1, 2, 0)
    text = path.read_text()
    assert text.count("CYLINDRICAL_SURFACE(") == 1
    assert text.count("SHELL_BASED_SURFACE_MODEL(") == 2
    occ = _occ(path)
    if occ is not None:
        assert occ[1]
        # The measured radius and length: pi r^2 L to the writer's precision.
        assert _solid_volumes(occ[0]) == [pytest.approx(math.pi * 4.0 * 5.0, rel=1e-9)]


def test_a_sphere_is_an_exact_ball(tmp_path):
    path = tmp_path / "ball.step"
    rep = mp.write_step_elements(path, _fitted(mp.icosa_sphere(1.5, 3)))
    assert rep["num_spheres"] == 1
    assert "SPHERICAL_SURFACE(" in path.read_text()
    occ = _occ(path)
    if occ is not None:
        assert occ[1]
        assert _solid_volumes(occ[0]) == [
            pytest.approx(4 / 3 * math.pi * 1.5**3, rel=1e-6)
        ]


def test_selection_and_refusals(tmp_path):
    fitted = _fitted(mp.extract_skin(mp.box((2, 2, 2))))
    path = tmp_path / "planes.step"
    assert mp.write_step_elements(path, fitted)["num_planes"] == 6
    assert mp.write_step_elements(path, fitted, patches=[0])["num_planes"] == 1
    with pytest.raises(mp.WriteError, match="no fitted patch 99"):
        mp.write_step_elements(path, fitted, patches=[99])
    with pytest.raises(mp.WriteError, match="no fitted patches"):
        mp.write_step_elements(path, mp.icosa_sphere())
    occ = _occ(path)
    assert occ is None or occ[1]


def test_sketches_and_sections_are_curves(tmp_path):
    outline = mp.silhouette(mp.annulus(0.5, 1.0, 64, 2, cell_type="triangle"))
    sketch = mp.fit_curves(outline, tolerance=1e-3)
    path = tmp_path / "sketch.step"
    mp.write(path, sketch)
    text = path.read_text()
    assert "GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION(" in text
    assert text.count("CIRCLE(") == 2
    mp.write(tmp_path / "poly.step", outline)
    assert (tmp_path / "poly.step").read_text().count("POLYLINE(") == 2
    occ = _occ(path)
    assert occ is None or occ[1]


def test_a_slot_sketch_has_lines_and_arcs(tmp_path):
    t = np.linspace(-np.pi / 2, np.pi / 2, 25)
    ring = [(2 + 0.75 * np.cos(a), 0.75 * np.sin(a)) for a in t] + [
        (-2 - 0.75 * np.cos(a), -0.75 * np.sin(a)) for a in t
    ]
    n = len(ring)
    pts = [(x, y, 0.0) for x, y in ring] + [(x, y, 1.0) for x, y in ring]
    tris = [
        t
        for i in range(n)
        for t in ((i, (i + 1) % n, (i + 1) % n + n), (i, (i + 1) % n + n, i + n))
    ]
    prism = mp.Mesh(np.array(pts), [("triangle", np.array(tris))])
    section = mp.slice(prism, (0, 0, 0.5), (0, 0, 1), chain=True)
    sketch = mp.fit_curves(section, tolerance=1e-6)
    path = tmp_path / "slot.step"
    mp.write(path, sketch)
    text = path.read_text()
    assert text.count("TRIMMED_CURVE(") == 4  # two lines, two arcs
    assert text.count("LINE(") - text.count("POLYLINE(") == 2
