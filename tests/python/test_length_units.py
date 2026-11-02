"""Length units as metadata (roadmap 7.1.4): ``field_data["meshio:length_unit"]``
(metres per unit) filled by the UNV, Radioss, USD, CGNS, MED, STEP and IGES readers,
honoured by their writers and the STEP writer, and updated by ``transform``'s
``scale_units`` -- in both engines."""

import importlib.util
import json
import pathlib
import subprocess

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core, _units
from meshioplusplus.cgns import _cgns
from meshioplusplus.med import _med
from meshioplusplus.radioss import _radioss
from meshioplusplus.step import _step
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


# --------------------------------------------------------------------------- #
# CGNS: a DimensionalUnits_t node under the base                             #
# --------------------------------------------------------------------------- #
HAS_HDF5 = hasattr(_core, "cgns_write")
needs_hdf5 = pytest.mark.skipif(not HAS_HDF5, reason="needs the HDF5-enabled core")


def _cgns_units(path):
    """The five SIDS unit names of the base, or None without the node."""
    import h5py

    with h5py.File(path, "r") as f:
        base = f["Base"]
        if "DimensionalUnits" not in base:
            return None
        node = base["DimensionalUnits"]
        assert bytes(node.attrs["label"]).rstrip(b"\0") == b"DimensionalUnits_t"
        assert bytes(node.attrs["type"]).rstrip(b"\0") == b"C1"
        data = node[" data"][()]
        assert data.shape == (5, 32)
        return [bytes(row).decode().rstrip() for row in data]


@needs_hdf5
@pytest.mark.parametrize(
    "unit,name",
    [(1.0, "Meter"), (0.01, "Centimeter"), (1e-3, "Millimeter"), (0.0254, "Inch")],
)
def test_cgns_writes_and_reads_dimensional_units_in_both_engines(tmp_path, unit, name):
    mesh = _tets(unit)
    core, py = tmp_path / "core.cgns", tmp_path / "py.cgns"
    _core.cgns_write(str(core), mesh, -1)
    _cgns.write(py, mesh, compression=None)
    for path in (core, py):
        assert _cgns_units(path) == ["Null", name, "Null", "Null", "Null"]
        for back in (_core.cgns_read(str(path)), _cgns.read(path)):
            assert _unit(back) == pytest.approx(unit, rel=1e-15)
        # The unit lives in the node, not again in the zone's FieldData.
        import h5py

        with h5py.File(path, "r") as f:
            assert KEY not in f["Base/Zone1"].get("FieldData", {})


@needs_hdf5
@pytest.mark.parametrize("engine", ["core", "python"])
def test_cgns_unnamed_unit_stays_field_data(tmp_path, engine):
    path = tmp_path / "x.cgns"
    if engine == "core":
        _core.cgns_write(str(path), _tets(0.5), -1)
    else:
        with pytest.warns(UserWarning, match="DimensionalUnits name"):
            _cgns.write(path, _tets(0.5), compression=None)
    assert _cgns_units(path) is None
    assert _unit(mp.read(str(path))) == pytest.approx(0.5)


@needs_hdf5
def test_cgns_zone_units_override_the_base(tmp_path):
    import h5py

    path = tmp_path / "zone.cgns"
    _cgns.write(path, _tets(1.0), compression=None)
    with h5py.File(path, "r+") as f:
        f.copy(f["Base/DimensionalUnits"], f["Base/Zone1"], "DimensionalUnits")
        row = np.frombuffer(b"Millimeter".ljust(32), dtype="i1")
        f["Base/Zone1/DimensionalUnits/ data"][1] = row
    for back in (_core.cgns_read(str(path)), _cgns.read(path)):
        assert _unit(back) == pytest.approx(1e-3)


# --------------------------------------------------------------------------- #
# MED: the coordinate units (UNI)                                             #
# --------------------------------------------------------------------------- #
def _med_write(engine, path, mesh):
    if engine == "core":
        _core.med_write(
            str(path),
            mesh,
            {},
            {},
            [],
            "mesh",
            "",
            "",
            getattr(mesh, "unit_coords", ""),
            {},
            {},
            "4.1.0",
        )
    else:
        _med.write(str(path), mesh)


@needs_hdf5
@pytest.mark.parametrize("engine", ["core", "python"])
@pytest.mark.parametrize(
    "unit,name", [(1.0, "m"), (1e-3, "mm"), (1e-6, "um"), (0.3048, "ft")]
)
def test_med_writes_and_reads_its_coordinate_units(tmp_path, engine, unit, name):
    path = tmp_path / "x.med"
    _med_write(engine, path, _tets(unit))
    for back in (_core.med_read(str(path)), _med.read(str(path))):
        assert back.unit_coords.split() == [name] * 3
        assert _unit(back) == pytest.approx(unit, rel=1e-15)


@needs_hdf5
@pytest.mark.parametrize("engine", ["core", "python"])
def test_med_explicit_unit_coords_win(tmp_path, engine):
    mesh = _tets(1e-3)
    mesh.unit_coords = "cm".ljust(16) * 3
    path = tmp_path / "x.med"
    _med_write(engine, path, mesh)
    back = mp.read(str(path))
    assert _unit(back) == pytest.approx(0.01)
    # Other tools spell micrometres with the Latin-1 micro sign.
    mesh.unit_coords = "\u00b5m".ljust(16) * 3
    _med_write(engine, path, mesh)
    back = _core.med_read(str(path)) if engine == "core" else _med.read(str(path))
    assert _unit(back) == pytest.approx(1e-6)
    # Mixed or unknown axis units set no length unit.
    mesh = _tets()
    mesh.unit_coords = "m".ljust(16) + "mm".ljust(16) + "m".ljust(16)
    mp.write(str(path), mesh)
    assert _unit(mp.read(str(path))) is None


# --------------------------------------------------------------------------- #
# STEP: the writer declares the unit                                          #
# --------------------------------------------------------------------------- #
def _tet_surface(unit=None):
    tris = np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]])
    pts = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    mesh = mp.Mesh(pts, [("triangle", tris)])
    if unit is not None:
        _units.set_length_unit(mesh, unit)
    return mesh


@pytest.mark.parametrize(
    "unit,entity",
    [
        (None, "SI_UNIT(.MILLI.,.METRE.)"),
        (1.0, "SI_UNIT($,.METRE.)"),
        (0.01, "SI_UNIT(.CENTI.,.METRE.)"),
        (1e-6, "SI_UNIT(.MICRO.,.METRE.)"),
        (0.0254, "CONVERSION_BASED_UNIT('INCH',"),
        (0.3048, "CONVERSION_BASED_UNIT('FOOT',"),
        (0.5, "CONVERSION_BASED_UNIT('LENGTH UNIT',"),
    ],
)
def test_step_declares_the_length_unit_in_both_engines(tmp_path, unit, entity):
    mesh = _tet_surface(unit)
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    core, py = tmp_path / "a" / "x.step", tmp_path / "b" / "x.step"
    _step.write(str(py), mesh)
    text = py.read_text()
    assert entity in text
    if unit == 0.0254:
        assert float(text.split("LENGTH_MEASURE(")[1].split(")")[0]) == 25.4
    if hasattr(_core, "step_write"):
        _core.step_write(str(core), mesh, "facet", 0.0)
        assert core.read_bytes() == py.read_bytes()
        # The B-rep writers declare it too.
        _core.step_write(str(core), mesh, "planar", 0.0)
        assert entity in core.read_text()


# --------------------------------------------------------------------------- #
# STEP and IGES reads: OpenCASCADE returns millimetres, and says so           #
# --------------------------------------------------------------------------- #
def _has_occt():
    from meshioplusplus._fit_cad import has_occt

    return has_occt()


needs_occt = pytest.mark.skipif(not _has_occt(), reason="needs OpenCASCADE")
BOX_IGS = pathlib.Path(__file__).parents[1] / "wasm" / "fixtures" / "box.igs"


@needs_occt
@pytest.mark.parametrize("unit", [None, 1.0, 0.0254, 0.3048, 1e-3])
def test_step_read_states_millimetres(tmp_path, unit):
    path = tmp_path / "x.step"
    mp.write(str(path), _tet_surface(unit), file_format="step")
    back = mp.read(str(path))
    assert _unit(back) == 1e-3
    # The physical size survives: a unit edge is (unit / 1 mm) millimetres.
    scale = (unit if unit is not None else 1e-3) / 1e-3
    assert back.points.max() == pytest.approx(scale, rel=1e-9)
    # transform's scale_units takes it to metres.
    metres = mp.transform(back, scale_units=1e-3)
    assert _unit(metres) == pytest.approx(1.0)
    assert metres.points.max() == pytest.approx(scale * 1e-3, rel=1e-9)


@needs_occt
def test_iges_read_states_millimetres():
    back = mp.read(str(BOX_IGS))
    assert _unit(back) == 1e-3
