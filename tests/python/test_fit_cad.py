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

from .analytic_surfaces import box_skin, torus, turned_part

HAS_OCCT = _fit_cad.has_occt()
needs_occt = pytest.mark.skipif(
    not HAS_OCCT, reason="core built without -DMESHIOPLUSPLUS_WITH_OCCT=ON"
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
    with pytest.raises(NotImplementedError, match="MESHIOPLUSPLUS_WITH_OCCT"):
        mio.fit_cad(box_skin(2))


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
    assert report["num_faceted_faces"] == report["num_faces"] > 6
    assert report["num_demoted"] == 0


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
    assert report["num_plane"] == 6
    assert dst.exists() and solid.exists()


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
