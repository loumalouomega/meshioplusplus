"""Length units as metadata (roadmap 7.1.4): ``field_data["meshio:length_unit"]``
(metres per unit) filled by the UNV, Radioss and USD readers, honoured by their
writers, and updated by ``transform``'s ``scale_units`` -- in both engines."""

import importlib.util
import json
import pathlib
import subprocess

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core, _units
from meshioplusplus.radioss import _radioss
from meshioplusplus.unv import _unv

from . import helpers

KEY = "meshio:length_unit"
UNV_MM = pathlib.Path(__file__).parent / "meshes" / "unv" / "units_cs_legacy.unv"


def _tets(unit=None):
    mesh = mp.Mesh(helpers.tet_mesh.points.copy(), helpers.tet_mesh.cells)
    if unit is not None:
        _units.set_length_unit(mesh, unit)
    return mesh


def _unit(mesh):
    return _units.length_unit(mesh)


def test_helpers():
    mesh = _tets()
    assert _unit(mesh) is None
    _units.set_length_unit(mesh, 0.001)
    assert mesh.field_data[KEY].shape == (1,) and _unit(mesh) == 0.001
    _units.set_length_unit(mesh, None)
    assert KEY not in mesh.field_data
    with pytest.raises(ValueError, match="positive"):
        _units.set_length_unit(mesh, -1.0)
    mesh.field_data[KEY] = np.array([np.nan])
    assert _unit(mesh) is None


# --------------------------------------------------------------------------- #
# UNV                                                                         #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("engine", ["core", "python"])
def test_unv_reads_its_length_unit(engine):
    if UNV_MM.read_bytes().startswith(b"version https://git-lfs"):
        pytest.skip("LFS fixture not fetched")
    mesh = (
        _core.unv_read(str(UNV_MM), False, None, 0)
        if engine == "core"
        else _unv.read(str(UNV_MM))
    )
    assert _unit(mesh) == pytest.approx(1e-3, rel=1e-15)


@pytest.mark.parametrize("unit,code", [(1.0, 1), (1e-3, 5), (0.0254, 9), (None, None)])
def test_unv_writes_the_length_unit_in_both_engines(tmp_path, unit, code):
    mesh = _tets(unit)
    core, py = tmp_path / "core.unv", tmp_path / "py.unv"
    _core.unv_write(str(core), mesh, code_aster=False, node_dataset=2411)
    _unv.write(str(py), mesh)
    assert core.read_bytes() == py.read_bytes()
    text = py.read_text()
    if unit is None:
        assert "\n   164\n" not in text
        assert KEY not in mp.read(str(py)).field_data
        return
    assert f"\n   164\n{code:10d}" in text
    for back in (mp.read(str(py)), _unv.read(str(py))):
        assert _unit(back) == pytest.approx(unit, rel=1e-15)


def test_unv_keeps_its_own_units_when_they_agree(tmp_path):
    mesh = _tets(1e-3)
    mesh.field_data["unv:units"] = np.array([5], dtype=np.int64)
    mesh.field_data["unv:unit_factors"] = np.array([1e3, 1e3, 1.0, 273.15])
    path = tmp_path / "a.unv"
    mp.write(str(path), mesh)
    assert "\n   164\n         5" in path.read_text()
    # ... and are rewritten when the shared unit says otherwise (scaled to m).
    scaled = mp.transform(mesh, scale_units=1e-3)
    assert _unit(scaled) == pytest.approx(1.0)
    mp.write(str(path), scaled)
    assert "\n   164\n         1" in path.read_text()


# --------------------------------------------------------------------------- #
# Radioss                                                                     #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "unit,name",
    [(None, "m"), (1.0, "m"), (1e-3, "mm"), (0.0254, "in"), (5e-4, "0.0005")],
)
def test_radioss_states_and_reads_the_unit_in_both_engines(tmp_path, unit, name):
    mesh = _tets(unit)
    (tmp_path / "core").mkdir()
    (tmp_path / "py").mkdir()
    # One name: the run title is the file's.
    core, py = tmp_path / "core" / "x_0000.rad", tmp_path / "py" / "x_0000.rad"
    _core.radioss_write(str(core), mesh, False)
    _radioss.write(str(py), mesh)
    assert core.read_bytes() == py.read_bytes()
    begin = py.read_text().splitlines()
    units_line = begin[begin.index("/BEGIN") + 3]
    assert units_line == f"{'kg':>20}{name:>20}{'s':>20}"
    expected = 1.0 if unit is None else unit
    for back in (_core.radioss_read(str(py)), _radioss.read(str(py))):
        assert _unit(back) == pytest.approx(expected, rel=1e-12)
        np.testing.assert_allclose(back.points, mesh.points, atol=1e-12)


def test_radioss_converted_input_reads_in_the_work_unit(tmp_path):
    """Input in mm, work in m: coordinates are converted, the unit is metres."""
    mesh = _tets(1e-3)
    path = tmp_path / "x_0000.rad"
    _radioss.write(str(path), mesh)
    lines = path.read_text().splitlines()
    k = lines.index("/BEGIN")
    lines[k + 4] = f"{'kg':>20}{'m':>20}{'s':>20}"
    path.write_text("\n".join(lines) + "\n")
    for back in (_core.radioss_read(str(path)), _radioss.read(str(path))):
        assert _unit(back) == pytest.approx(1.0)
        np.testing.assert_allclose(back.points, mesh.points * 1e-3, atol=1e-15)


# --------------------------------------------------------------------------- #
# transform and the pipeline                                                  #
# --------------------------------------------------------------------------- #
def test_transform_scale_units_updates_the_unit():
    mm = _tets(1e-3)
    m = mp.transform(mm, scale_units=1e-3)
    assert _unit(m) == pytest.approx(1.0)
    np.testing.assert_allclose(m.points, mm.points * 1e-3)
    assert _unit(mm) == 1e-3  # the input is untouched
    # Other transforms, and meshes without a unit, state nothing new.
    assert _unit(mp.transform(mm, scale=2.0)) == 1e-3
    assert _unit(mp.transform(_tets(), scale_units=1e-3)) is None


def test_pipeline_scale_units_in_both_engines(tmp_path):
    src = tmp_path / "in.vtu"
    mp.write(str(src), _tets(1e-3), compression=None)
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "Transform", "ScaleUnits": 0.001}],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    mp.run_pipeline(doc)
    assert _unit(mp.read(str(tmp_path / "py.vtu"))) == pytest.approx(1.0)
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    _core.run_pipeline_json(json.dumps(doc))
    assert _unit(mp.read(str(tmp_path / "cpp.vtu"))) == pytest.approx(1.0)
    out, _ = mp.apply_operations(_tets(1e-3), [{"Op": "Transform", "ScaleUnits": 1e-3}])
    assert _unit(out) == pytest.approx(1.0)


def test_native_cli_scale_units(tmp_path):
    from .test_cli_json import NATIVE

    if NATIVE is None:
        pytest.skip("no native CLI binary built")
    src, out = tmp_path / "in.vtu", tmp_path / "out.vtu"
    mp.write(str(src), _tets(1e-3), compression=None)
    r = subprocess.run(
        [NATIVE, "transform", str(src), str(out), "--scale-units", "0.001"],
        capture_output=True,
    )
    assert r.returncode == 0, r.stderr
    assert _unit(mp.read(str(out))) == pytest.approx(1.0)


# --------------------------------------------------------------------------- #
# USD: the roadmap's probe                                                    #
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(importlib.util.find_spec("pxr") is None, reason="needs usd-core")
def test_a_millimetre_unv_converts_to_usd_in_millimetres(tmp_path):
    from pxr import Usd, UsdGeom

    path = tmp_path / "mm.unv"
    mp.write(str(path), _tets(1e-3))
    usd = tmp_path / "mm.usda"
    mp.write(str(usd), mp.read(str(path)))
    stage = Usd.Stage.Open(str(usd))
    assert UsdGeom.GetStageMetersPerUnit(stage) == pytest.approx(1e-3)
    assert _unit(mp.read(str(usd))) == pytest.approx(1e-3)
    # An explicit argument still wins, and a unit-less mesh is written in metres.
    explicit, plain = tmp_path / "explicit.usda", tmp_path / "plain.usda"
    mp.write(str(explicit), _tets(1e-3), meters_per_unit=1.0)
    assert UsdGeom.GetStageMetersPerUnit(Usd.Stage.Open(str(explicit))) == 1.0
    mp.write(str(plain), _tets())
    assert UsdGeom.GetStageMetersPerUnit(Usd.Stage.Open(str(plain))) == 1.0
