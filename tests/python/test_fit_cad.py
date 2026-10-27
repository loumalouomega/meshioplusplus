"""fit_cad: analytic CAD faces through the optional OpenCASCADE backend.

The OpenCASCADE cases run only against a core built with
``-DMESHIOPLUSPLUS_WITH_OCCT=ON``; elsewhere the named ``NotImplementedError``
is the contract. The STEP read-back check also needs ``OCP`` (``pip install
cadquery-ocp``), an independent OpenCASCADE.
"""

import json

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _fit_cad
from meshioplusplus._cli._main import main as cli_main

from .analytic_surfaces import (
    box_skin,
    filleted_plate,
    inclined_boss,
    rounded_box,
    rounded_thread_rod,
    section_blend_block,
    tapered_rounded_thread_rod,
    tapered_thread_rod,
    threaded_rod,
    torus,
    turned_part,
    variable_fillet_block,
    variable_pitch_thread_rod,
    wavy_bar,
    wavy_capped,
)

HAS_OCCT = _fit_cad.has_occt()
needs_occt = pytest.mark.skipif(
    not HAS_OCCT, reason="core built without -DMESHIOPLUSPLUS_WITH_OCCT=ON"
)
try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python install
    _core = None
needs_refine = pytest.mark.skipif(
    not hasattr(_core, "refine_patches"), reason="core predates refine_patches"
)


def _fitted(skin, curvature_tolerance=0.0, tolerance=1e-3, min_patch_cells=1):
    seg = mio.segment(
        skin,
        curvature_tolerance=curvature_tolerance,
        min_patch_cells=min_patch_cells,
    )
    return mio.fit_primitives(seg, tolerance=tolerance)


@pytest.mark.skipif(HAS_OCCT, reason="checks the build without OpenCASCADE")
def test_absent_backend_is_named():
    # The automatic backend falls back to the native one; asking for
    # OpenCASCADE by name is what raises.
    _, report = mio.fit_cad(box_skin(2), return_report=True)
    assert report["backend"] == "native"
    with pytest.raises(NotImplementedError, match="MESHIOPLUSPLUS_WITH_OCCT"):
        mio.fit_cad(box_skin(2), backend="occt")


def test_bad_arguments_are_refused_before_the_core():
    with pytest.raises(ValueError, match="unknown primitive type"):
        mio.fit_cad(box_skin(2), types=["blob"])
    with pytest.raises(ValueError, match="tolerance"):
        mio.fit_cad(box_skin(2), tolerance=-1.0)
    with pytest.raises(ValueError, match="max_demotion_rounds"):
        mio.fit_cad(box_skin(2), max_demotion_rounds=-1)


@needs_occt
def test_box_is_six_exact_planes(tmp_path):
    path = tmp_path / "box.step"
    out, report = mio.fit_cad(_fitted(box_skin(3)), path, return_report=True)
    assert report["valid"]
    assert report["num_faces"] == report["num_plane"] == 6
    assert report["num_exact_edges"] == 12
    assert report["num_demoted"] == 0
    assert report["volume_ratio"] == pytest.approx(1.0, abs=1e-9)
    assert report["max_deviation"] < 1e-9
    assert path.read_text().startswith("ISO-10303-21;")
    assert out.cells[0].type == "triangle"
    for name in ("fit_cad:face", "fit_cad:type", "fit_cad:patch"):
        assert name in out.cell_data
    assert set(out.cell_data["fit_cad:type"][0].tolist()) == {1}


@needs_occt
@pytest.mark.parametrize("noise", [0.0, 1e-3])
def test_turned_part_keeps_every_primitive(noise):
    skin = turned_part(noise=noise)
    fit = _fitted(skin, curvature_tolerance=1.0, tolerance=2e-3, min_patch_cells=20)
    _out, report = mio.fit_cad(fit, return_report=True)
    assert report["valid"]
    assert report["num_demoted"] == 0
    assert report["num_faceted_faces"] == 0
    assert [
        report["num_" + t] for t in ("plane", "sphere", "cylinder", "cone", "torus")
    ] == [
        2,
        1,
        1,
        1,
        1,
    ]
    assert report["max_deviation"] <= report["tolerance"]
    assert report["volume_ratio"] == pytest.approx(1.0, abs=0.01)


@pytest.mark.parametrize("backend", ["native", "occt"])
def test_constrained_primitives_share_one_axis(backend):
    if backend == "occt" and not HAS_OCCT:
        pytest.skip("needs a core built with OpenCASCADE")
    # Seed 1: the constrained solid keeps every face on both backends.
    skin = turned_part(noise=1e-3, seed=1)
    fit = _fitted(skin, curvature_tolerance=1.0, tolerance=2e-3, min_patch_cells=20)
    _out, plain = mio.fit_cad(fit, backend=backend, return_report=True)
    assert "num_constrained_primitives" not in plain
    _out, report = mio.fit_cad(
        fit, backend=backend, constrain_primitives=True, return_report=True
    )
    # Two planes, a sphere, a cylinder, a cone and a torus about one axis.
    assert report["num_constrained_primitives"] == 6
    assert report["num_constraint_groups"] >= 2
    assert report["num_demoted"] == 0
    assert report["num_faces"] == plain["num_faces"]
    assert report["max_deviation"] <= report["tolerance"]


@needs_occt
def test_constrained_primitives_fall_back_when_the_solid_loses_a_junction():
    # Seed 7: the constrained torus fillet no longer builds its tangent junction
    # with the wall in OpenCASCADE, which demotes two faces. The unconstrained
    # solid, which demotes none, is returned instead, with the counts at zero.
    skin = turned_part(noise=1e-3)
    fit = _fitted(skin, curvature_tolerance=1.0, tolerance=2e-3, min_patch_cells=20)
    _out, plain = mio.fit_cad(fit, backend="occt", return_report=True)
    _out, report = mio.fit_cad(
        fit, backend="occt", constrain_primitives=True, return_report=True
    )
    assert report["valid"] is True
    assert report["num_constrained_primitives"] == 0
    assert report["num_constraint_groups"] == 0
    assert report["num_demoted"] == plain["num_demoted"] == 0
    assert report["num_faces"] == plain["num_faces"] == 6


@needs_occt
def test_closed_torus_is_one_natural_face():
    _out, report = mio.fit_cad(_fitted(torus()), return_report=True)
    assert report["valid"]
    assert report["num_faces"] == report["num_torus"] == 1


@needs_occt
def test_types_mask_and_empty_types_facet_the_rest():
    fit = _fitted(box_skin(2))
    _out, report = mio.fit_cad(fit, types=[], return_report=True)
    assert report["valid"]
    assert report["num_plane"] == 0
    # The facets left out merge into one planar face per coplanar group.
    assert report["num_faceted_faces"] == report["num_faces"] == 6
    assert report["num_demoted"] == 0


def _bumped_box(n=32, amp=0.15):
    """The unit cube whose top carries a smooth bump no primitive explains."""
    skin = box_skin(n)
    p = skin.points.copy()
    top = np.isclose(p[:, 2], 1.0)
    x, y = p[top, 0], p[top, 1]
    p[top, 2] += (
        amp
        * np.exp(-((x - 0.5) ** 2 + (y - 0.5) ** 2) / 0.03)
        * 16
        * x
        * (1 - x)
        * y
        * (1 - y)
    )
    return mio.Mesh(p, skin.cells)


@needs_occt
def test_unexplained_patch_becomes_a_freeform_face(tmp_path):
    fit = _fitted(_bumped_box())
    assert 0 in fit.field_data["fit:patch_type"].tolist()  # the bump is no primitive
    out, report = mio.fit_cad(fit, tmp_path / "bump.step", return_report=True)
    assert report["valid"]
    assert report["num_bspline"] == 1 and report["num_plane"] == 5
    assert report["num_faceted_faces"] == 0 and report["num_demoted"] == 0
    assert report["max_deviation"] <= report["tolerance"]
    assert report["volume_ratio"] == pytest.approx(1.0, abs=1e-4)
    assert 6 in out.cell_data["fit_cad:type"][0].tolist()
    # The solid as NURBS, ready for Kratos: six trimmed faces.
    brep = out.info["brep"]
    assert len(brep["faces"]) == 6 and all(f["trimmed"] for f in brep["faces"])
    mio.write(tmp_path / "bump.cad.json", out)
    assert len(mio.read(tmp_path / "bump.cad.json").info["brep"]["faces"]) == 6


@needs_occt
def test_freeform_off_or_a_coarse_skin_stays_faceted():
    _out, report = mio.fit_cad(
        _fitted(_bumped_box()), freeform=False, return_report=True
    )
    assert report["num_bspline"] == 0 and report["num_faceted_faces"] > 0
    # A skin coarser than the tolerance cannot be a smooth face within it: the
    # fit refuses it up front rather than building a face to demote.
    _out, report = mio.fit_cad(
        _fitted(_bumped_box(16)), tolerance=2e-3, return_report=True
    )
    assert report["num_bspline"] == 0 and report["num_demoted"] == 0


@needs_occt
def test_a_curled_patch_is_split_into_freeform_faces():
    """The wavy wall folds over its plane and, closed around the axis, over
    its unrolled cylinder too: it is split into pieces that each unroll
    (roadmap 6.3.2), each a B-spline face under a new patch label."""
    fit = _fitted(wavy_capped())
    out, report = mio.fit_cad(fit, return_report=True)
    assert report["valid"]
    assert report["num_split_patches"] >= 2
    assert report["num_bspline"] == report["num_split_patches"]
    assert report["num_faceted_faces"] == 0
    assert report["max_deviation"] <= report["tolerance"]
    largest = max(int(np.max(b)) for b in fit.cell_data["segment:patch"])
    labels = set(out.cell_data["fit_cad:patch"][0].tolist())
    assert max(labels) == largest + report["num_split_patches"]
    _out, report = mio.fit_cad(fit, freeform_split=False, return_report=True)
    assert report["num_split_patches"] == 0 and report["num_bspline"] == 0


@needs_occt
def test_a_rolling_ball_fillet_is_its_own_face():
    """The foot of a leaning boss rounded with a constant radius: a canal
    surface along an ellipse, no primitive (roadmap 6.3.1)."""
    fit = _fitted(inclined_boss(), curvature_tolerance=3.0)
    out, report = mio.fit_cad(fit, return_report=True)
    assert report["valid"]
    assert report["num_fillet"] == 1 and report["num_faceted_faces"] == 0
    assert report["num_variable_fillet"] == 0 and report["num_blend"] == 0
    assert report["max_deviation"] <= report["tolerance"]
    assert 7 in out.cell_data["fit_cad:type"][0].tolist()
    # Every face reaches the NURBS side channel, trimmed (OpenCASCADE 7.9's
    # conversion throws on this solid; the faces are then approximated).
    faces = out.info["brep"]["faces"]
    assert len(faces) == report["num_faces"] and all(f["trimmed"] for f in faces)
    _out, report = mio.fit_cad(fit, fillets=False, return_report=True)
    assert report["num_fillet"] == 0


@needs_occt
def test_a_fillet_of_varying_radius_is_its_own_face():
    """A block's edge rounded by a ball whose radius goes from 0.3 to 0.2
    along it: no one radius explains it, a radius law along the rail does
    (roadmap 6.3.1)."""
    fit = _fitted(variable_fillet_block(), curvature_tolerance=3.0)
    out, report = mio.fit_cad(fit, return_report=True)
    assert report["valid"]
    assert report["num_fillet"] == report["num_variable_fillet"] == 1
    assert report["num_plane"] == 6 and report["num_faceted_faces"] == 0
    assert report["max_deviation"] <= report["tolerance"]
    assert report["volume_ratio"] == pytest.approx(1.0, abs=0.01)
    assert set(out.cell_data["fit_cad:type"][0].tolist()) == {1, 7}
    _out, report = mio.fit_cad(fit, fillets=False, return_report=True)
    assert report["num_fillet"] == report["num_variable_fillet"] == 0


@needs_occt
def test_a_blend_between_two_freeform_faces_is_its_own_face():
    """A bar whose wavy top and leaning front are both freeform B-spline
    faces, its edge rounded by a rolling ball: the blend is fitted after the
    freeform faces, against them (roadmap 6.3.1)."""
    fit = _fitted(wavy_bar(), curvature_tolerance=3.0)
    out, report = mio.fit_cad(fit, return_report=True)
    assert report["valid"]
    assert report["num_fillet"] == report["num_blend"] == 1
    assert report["num_bspline"] == 2 and report["num_faceted_faces"] == 0
    assert report["max_deviation"] <= report["tolerance"]
    assert {6, 7} <= set(out.cell_data["fit_cad:type"][0].tolist())
    _out, report = mio.fit_cad(fit, fillets=False, return_report=True)
    assert report["num_fillet"] == report["num_blend"] == 0


@needs_occt
def test_thread_flanks_are_helical_faces():
    """A trapezoidal thread: both flanks are ruled helicoids about the axis
    of the root and crest cylinders (roadmap 6.3.1)."""
    fit = _fitted(threaded_rod(turns=3))
    out, report = mio.fit_cad(fit, return_report=True)
    assert report["valid"]
    assert report["num_helical"] == 2 and report["num_cylinder"] == 2
    assert report["num_rounded_thread"] == 0
    assert report["max_deviation"] <= report["tolerance"]
    assert 8 in out.cell_data["fit_cad:type"][0].tolist()
    faces = out.info["brep"]["faces"]
    assert len(faces) == report["num_faces"] and all(f["trimmed"] for f in faces)


@needs_occt
def test_a_rounded_thread_with_no_cylinder_is_helical():
    """Root and crest rounded into arcs and no cylinder anywhere: every face
    of the thread finds the screw from its own facets, the flanks straight
    and the root and crest rounded (roadmap 6.3.1)."""
    fit = _fitted(rounded_thread_rod(turns=3))
    out, report = mio.fit_cad(fit, return_report=True)
    assert report["valid"]
    assert report["num_helical"] == 4 and report["num_rounded_thread"] == 2
    assert report["num_cylinder"] == 0
    assert report["max_deviation"] <= report["tolerance"]
    assert report["volume_ratio"] == pytest.approx(1.0, abs=0.01)
    assert 8 in out.cell_data["fit_cad:type"][0].tolist()


@needs_occt
def test_tapered_thread_flanks_are_helical_faces():
    """threaded_rod cut on a cone: root and crest are cone strips, the flanks
    ruled helicoids reported tapered (roadmap 6.3.3)."""
    out, report = mio.fit_cad(_fitted(tapered_thread_rod()), return_report=True)
    assert report["valid"]
    assert report["num_helical"] == report["num_tapered_thread"] == 2
    assert report["num_cone"] == 2 and report["num_variable_pitch_thread"] == 0
    assert report["max_deviation"] <= report["tolerance"]


@needs_occt
def test_a_thread_of_varying_pitch_is_helical():
    """The pitch grows by 0.05 per turn: the straight flanks are refitted with
    a pitch varying along them (roadmap 6.3.3)."""
    fit = _fitted(variable_pitch_thread_rod())
    out, report = mio.fit_cad(fit, return_report=True)
    assert report["valid"]
    assert report["num_helical"] == report["num_variable_pitch_thread"] == 2
    assert report["num_cylinder"] == 2 and report["num_demoted"] == 0
    assert report["max_deviation"] <= report["tolerance"]
    _out, report = mio.fit_cad(fit, threads=False, return_report=True)
    assert report["num_helical"] == 0


@needs_occt
def test_a_tapered_rounded_thread_finds_its_screw_and_taper():
    out, report = mio.fit_cad(_fitted(tapered_rounded_thread_rod()), return_report=True)
    assert report["valid"]
    assert report["num_helical"] == 4 and report["num_rounded_thread"] == 2
    assert report["num_tapered_thread"] == 2
    assert report["max_deviation"] <= report["tolerance"]


@needs_occt
def test_a_chamfer_of_varying_width_is_its_own_face():
    """A chamfer whose width grows along the edge: a twisted ruled surface no
    plane explains (roadmap 6.3.2)."""
    skin = section_blend_block(0.0)
    for fit in (_fitted(skin), mio.fit_primitives(skin, tolerance=1e-3)):
        out, report = mio.fit_cad(fit, return_report=True)
        assert report["valid"]
        assert report["num_chamfer"] == 1 and report["num_plane"] == 6
        assert report["num_faceted_faces"] == 0 and report["num_faces"] == 7
        assert report["max_deviation"] <= report["tolerance"]
        assert 10 in out.cell_data["fit_cad:type"][0].tolist()
    _out, report = mio.fit_cad(_fitted(skin), fillets=False, return_report=True)
    assert report["num_chamfer"] == 0


@needs_occt
def test_a_conic_blend_is_a_fillet_face():
    """Elliptic sections tangent to both faces with legs of unequal, varying
    length: neither a rolling ball nor a chamfer (roadmap 6.3.2)."""
    fit = mio.fit_primitives(section_blend_block(0.5), tolerance=1e-3)
    out, report = mio.fit_cad(fit, return_report=True)
    assert report["valid"]
    assert report["num_fillet"] == report["num_conic_blend"] == 1
    assert report["num_chamfer"] == 0 and report["num_faceted_faces"] == 0
    assert report["max_deviation"] <= report["tolerance"]


@needs_occt
def test_rounded_box_corners_are_vertex_blends():
    """A rounded box scaled unequally: twelve conic edge blends, and at each
    corner a vertex blend tangent to the three meeting there (roadmap 6.3.1)."""
    fit = mio.fit_primitives(rounded_box(), tolerance=1e-3)
    out, report = mio.fit_cad(fit, return_report=True)
    assert report["valid"]
    assert report["num_plane"] == 6 and report["num_conic_blend"] == 12
    assert report["num_vertex_blend"] == 8 and report["num_faces"] == 26
    assert report["num_faceted_faces"] == 0
    assert report["max_deviation"] <= report["tolerance"]
    assert {1, 7, 9} == set(out.cell_data["fit_cad:type"][0].tolist())
    _out, report = mio.fit_cad(fit, fillets=False, return_report=True)
    assert report["num_vertex_blend"] == report["num_conic_blend"] == 0


def test_freeform_options_are_checked_before_the_core():
    with pytest.raises(ValueError, match="freeform_degree"):
        mio.fit_cad(box_skin(2), freeform_degree=9)
    with pytest.raises(ValueError, match="freeform_max_control"):
        mio.fit_cad(box_skin(2), freeform_max_control=2)


@needs_occt
def test_unlabelled_skin_is_segmented_first():
    _out, report = mio.fit_cad(box_skin(2), return_report=True)
    assert report["num_plane"] == 6


@needs_occt
def test_open_skin_is_refused():
    open_skin = mio.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2]]))],
    )
    with pytest.raises(ValueError, match="closed"):
        mio.fit_cad(open_skin)


@needs_occt
def test_pipeline_step_in_both_engines(tmp_path):
    from meshioplusplus import _core

    src = tmp_path / "box.vtu"
    mio.write(src, box_skin(2), compression=None)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {"Op": "Segment"},
            {"Op": "FitPrimitives"},
            {"Op": "FitCad", "Path": str(tmp_path / "py.step")},
        ],
        "Output": {"Path": str(tmp_path / "out.vtu")},
    }
    reports = [mio.run_pipeline(settings)]
    if hasattr(_core, "run_pipeline_json"):
        settings["Operations"][-1]["Path"] = str(tmp_path / "cpp.step")
        settings["Output"]["Path"] = str(tmp_path / "out_cpp.vtu")
        reports.append(_core.run_pipeline_json(json.dumps(settings)))
        assert (tmp_path / "py.step").read_bytes().replace(b"py.step", b"x") == (
            tmp_path / "cpp.step"
        ).read_bytes().replace(b"cpp.step", b"x")
    for report in reports:
        entry = report["steps"][-1]
        assert entry["Valid"] == 1.0
        assert entry["NumFaces"] == 6
        assert entry["NumDemoted"] == 0
    assert "fit_cad:face" in mio.read(str(tmp_path / "out.vtu")).cell_data
    bad = dict(settings, Operations=[{"Op": "FitCad", "Types": ["blob"]}])
    with pytest.raises(Exception, match="unknown primitive type 'blob'"):
        mio.run_pipeline(bad)


@needs_refine
def test_refine_segmentation_explains_the_fillets():
    # Segmented by angles, the plate's corner fillet and the walls it is
    # tangent to share a patch; refining separates them.
    seg = mio.segment(filleted_plate())
    _out, report = mio.fit_cad(
        seg, refine_segmentation=True, backend="native", return_report=True
    )
    assert report["refined"] is True
    # The native backend builds only from the refined patches.
    assert report["refinement_kept"] is True
    assert report["explained_area_fraction"] == pytest.approx(1.0, abs=1e-9)
    assert report["explained_area_fraction"] > report["initial_explained_area_fraction"]
    assert report["num_refine_splits"] > 0
    assert report["refine_iterations"] >= 0 and report["num_refine_merges"] >= 0
    assert "refine:patch_proxy" in _out.field_data


@needs_refine
def test_default_report_has_no_refine_keys():
    _out, report = mio.fit_cad(
        mio.segment(filleted_plate()), backend="native", return_report=True
    )
    for key in (
        "refined",
        "refinement_kept",
        "refine_iterations",
        "num_refine_splits",
        "num_refine_merges",
        "initial_explained_area_fraction",
        "explained_area_fraction",
    ):
        assert key not in report


@needs_refine
def test_refine_segmentation_pipeline_step_in_both_engines(tmp_path):
    src = tmp_path / "plate.vtu"
    mio.write(src, filleted_plate(), compression=None)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {"Op": "Segment"},
            {
                "Op": "FitCad",
                "Path": str(tmp_path / "py.step"),
                "Backend": "native",
                "RefineSegmentation": True,
                "RansacSeed": 7,
            },
        ],
        "Output": {"Path": str(tmp_path / "out.vtu")},
    }
    reports = [mio.run_pipeline(settings)]
    if hasattr(_core, "run_pipeline_json"):
        settings["Operations"][-1]["Path"] = str(tmp_path / "cpp.step")
        settings["Output"]["Path"] = str(tmp_path / "out_cpp.vtu")
        reports.append(_core.run_pipeline_json(json.dumps(settings)))
        assert (tmp_path / "py.step").read_bytes().replace(b"py.step", b"x") == (
            tmp_path / "cpp.step"
        ).read_bytes().replace(b"cpp.step", b"x")
    for report in reports:
        entry = report["steps"][-1]
        assert entry["Refined"] == 1.0
        assert entry["ExplainedAreaFraction"] == pytest.approx(1.0, abs=1e-9)
        for key in ("NumRefineSplits", "NumRefineMerges", "RefinementKept"):
            assert key in entry
    if len(reports) == 2:
        assert list(reports[0]["steps"][-1]) == list(reports[1]["steps"][-1])
    # Without the option the step reports what it always did.
    settings["Operations"][-1].pop("RefineSegmentation")
    settings["Operations"][-1]["Path"] = str(tmp_path / "plain.step")
    entry = mio.run_pipeline(settings)["steps"][-1]
    assert "Refined" not in entry and "ExplainedAreaFraction" not in entry
    assert "RefinementKept" not in entry


_PHASE3_KEYS = {
    "freeform_general": (
        "num_freeform_general",
        "num_freeform_refinements",
        "avg_freeform_poles",
    ),
    "extend_surfaces": ("num_extended_edges",),
    "snap_junctions": ("num_snapped_junctions",),
    "face_deviation": ("hausdorff_max", "hausdorff_rms"),
}


def test_phase_three_report_keys_appear_only_when_asked():
    skin = box_skin(2)
    _out, plain = mio.fit_cad(skin, backend="native", return_report=True)
    # `closed` is always there; every other Phase 3 key only with its option.
    assert plain["closed"] is True
    for keys in _PHASE3_KEYS.values():
        for key in keys:
            assert key not in plain
    for option, keys in _PHASE3_KEYS.items():
        _out, report = mio.fit_cad(
            skin, backend="native", return_report=True, **{option: True}
        )
        assert set(report) - set(plain) == set(keys), option
    # The native backend builds no freeform faces and measures nothing.
    _out, report = mio.fit_cad(
        skin,
        backend="native",
        return_report=True,
        freeform_general=True,
        face_deviation=True,
    )
    assert report["num_freeform_general"] == 0
    assert report["avg_freeform_poles"] == 0.0
    assert np.isnan(report["hausdorff_max"]) and np.isnan(report["hausdorff_rms"])
    assert "fit_cad:face_hausdorff" not in _out.field_data


@pytest.mark.parametrize(
    "key, value",
    [
        ("freeform_tolerance", -1.0),
        ("freeform_tolerance", float("inf")),
        ("freeform_smoothing", -1.0),
        ("max_freeform_poles", 3),
        ("max_freeform_poles", 1025),
        ("max_freeform_iterations", 0),
        ("max_freeform_splits", -1),
        ("sew_tolerance", -1.0),
        ("sew_tolerance", float("nan")),
    ],
)
def test_phase_three_bad_options_are_refused(key, value):
    with pytest.raises(ValueError, match=key):
        mio.fit_cad(box_skin(2), backend="native", **{key: value})


def test_phase_three_options_reach_the_core():
    # Values in range but off their defaults pass through; the core's own
    # checks would refuse what the wrapper let by.
    _out, report = mio.fit_cad(
        box_skin(2),
        backend="native",
        return_report=True,
        freeform_general=True,
        freeform_tolerance=1e-4,
        freeform_smoothing=0.0,
        max_freeform_poles=4,
        max_freeform_iterations=1,
        max_freeform_splits=0,
        sew_tolerance=1e-6,
    )
    assert report["num_plane"] == 6 and report["closed"] is True


@needs_occt
def test_face_deviation_and_snapped_junctions():
    out, report = mio.fit_cad(
        _fitted(box_skin(3)),
        snap_junctions=True,
        face_deviation=True,
        return_report=True,
    )
    assert report["valid"] is True and report["closed"] is True
    assert report["num_snapped_junctions"] == 8  # the box's corners
    assert np.isfinite(report["hausdorff_max"])
    assert report["hausdorff_rms"] <= report["hausdorff_max"]
    # The per-face field data survives the conversion to a Python mesh.
    hausdorff = out.field_data["fit_cad:face_hausdorff"]
    rms = out.field_data["fit_cad:face_rms"]
    assert hausdorff.dtype == np.float64
    assert hausdorff.shape == (report["num_faces"], 2)
    assert rms.shape[0] == report["num_faces"]
    assert hausdorff.max() == report["hausdorff_max"]
    assert np.all(rms <= hausdorff.max(axis=1) + 1e-15)


def test_phase_three_pipeline_keys_in_both_engines(tmp_path):
    src = tmp_path / "box.vtu"
    mio.write(src, box_skin(2), compression=None)
    step = {
        "Op": "FitCad",
        "Backend": "native",
        "FreeformGeneral": True,
        "FreeformTolerance": 0.0,
        "FreeformSmoothing": 1e-10,
        "MaxFreeformPoles": 20,
        "MaxFreeformIterations": 5,
        "MaxFreeformSplits": 1,
        "SewTolerance": 0.0,
        "ExtendSurfaces": True,
        "SnapJunctions": True,
        "FaceDeviation": True,
        "ConstrainPrimitives": True,
    }
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [step],
        "Output": {"Path": str(tmp_path / "out.vtu")},
    }
    plain_settings = dict(settings, Operations=[{"Op": "FitCad", "Backend": "native"}])
    engines = [mio.run_pipeline]
    if hasattr(_core, "run_pipeline_json"):
        engines.append(lambda s: _core.run_pipeline_json(json.dumps(s)))
    entries = []
    for run in engines:
        plain = run(plain_settings)["steps"][-1]
        full = run(settings)["steps"][-1]
        # `Closed` is always reported, last of the unconditional counters.
        assert plain["Closed"] == 1.0
        assert list(plain)[-1] == "Closed"
        added = [k for k in full if k not in plain]
        assert added == [
            "NumFreeformGeneral",
            "NumFreeformRefinements",
            "AvgFreeformPoles",
            "NumExtendedEdges",
            "NumSnappedJunctions",
            "HausdorffMax",
            "HausdorffRms",
            "NumConstrainedPrimitives",
            "NumConstraintGroups",
        ]
        entries.append((list(plain), list(full)))
    # The two engines agree key for key, in order.
    assert all(e == entries[0] for e in entries)
    bad = dict(settings, Operations=[dict(step, MaxFreeformPoles=3)])
    for run in engines:
        with pytest.raises(Exception, match="max_freeform_poles"):
            run(bad)


@needs_occt
def test_to_cad_cli_json(tmp_path, capsys):
    src = tmp_path / "box.vtu"
    mio.write(src, box_skin(2))
    dst = tmp_path / "box.step"
    solid = tmp_path / "solid.vtu"
    assert (
        cli_main(["to-cad", str(src), str(dst), "--mesh-out", str(solid), "--json"])
        == 0
    )
    report = json.loads(capsys.readouterr().out)
    assert report["valid"] is True
    assert report["closed"] is True
    assert report["num_plane"] == 6
    assert dst.exists() and solid.exists()
    assert "num_snapped_junctions" not in report and "hausdorff_max" not in report
    args = ["to-cad", str(src), str(dst), "--json", "--snap-junctions"]
    assert cli_main(args + ["--face-deviation", "--sew-tolerance", "0"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["num_snapped_junctions"] == 8
    assert report["hausdorff_max"] >= 0.0
    assert "num_freeform_general" not in report


@needs_occt
def test_opencascade_reads_the_solid_back(tmp_path):
    pytest.importorskip("OCP")
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Reader

    skin = turned_part()
    fit = _fitted(skin, curvature_tolerance=1.0, tolerance=2e-3)
    path = tmp_path / "part.step"
    _out, report = mio.fit_cad(fit, path, return_report=True)
    reader = STEPControl_Reader()
    assert reader.ReadFile(str(path)) == IFSelect_RetDone
    reader.TransferRoots()
    shape = reader.OneShape()
    assert BRepCheck_Analyzer(shape).IsValid()
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    assert props.Mass() == pytest.approx(report["volume"], rel=1e-6)
