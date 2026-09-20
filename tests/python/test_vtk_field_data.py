"""``<FieldData>`` on VTU and VTP (v15.0.0).

Mesh-level ``field_data`` round-trips through both formats in both engines (the
C++ core under strict-core, and the pure-Python reference), and the layouts VTK
itself writes -- on the grid or inside a piece, ``TimeValue`` in particular -- are
read.
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._fallback import set_strict_core
from meshioplusplus.vtp import _vtp
from meshioplusplus.vtu import _vtu

XY = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]


def _mesh():
    mesh = meshioplusplus.Mesh(
        np.array(XY, dtype=np.float64), [("triangle", [[0, 1, 2]])]
    )
    mesh.point_data["u"] = np.arange(3.0)
    mesh.field_data["TimeValue"] = np.array([0.25])
    mesh.field_data["params"] = np.arange(6.0).reshape(2, 3)
    mesh.field_data["ids"] = np.arange(4, dtype=np.int32)
    return mesh


class _Engine:
    def __init__(self, vtu, vtp):
        self.vtu, self.vtp = vtu, vtp


@pytest.fixture(params=["core", "python"])
def engine(request):
    if request.param == "core":
        if not hasattr(_core, "vtu_write_codec"):
            pytest.skip("no C++ core")
        set_strict_core(True)
        yield _Engine(meshioplusplus.vtu, meshioplusplus.vtp)
        set_strict_core(None)
    else:
        yield _Engine(_vtu, _vtp)


def _same(back, mesh):
    assert sorted(back.field_data) == sorted(mesh.field_data)
    for key, value in mesh.field_data.items():
        got = back.field_data[key]
        assert got.dtype == np.asarray(value).dtype, key
        np.testing.assert_array_equal(got, value, err_msg=key)


@pytest.mark.parametrize(
    "binary, compression", [(False, None), (True, None), (True, "zlib")]
)
def test_vtu_roundtrips_field_data(engine, binary, compression, tmp_path):
    mesh = _mesh()
    engine.vtu.write(tmp_path / "a.vtu", mesh, binary=binary, compression=compression)
    _same(engine.vtu.read(tmp_path / "a.vtu"), mesh)


@pytest.mark.parametrize(
    "binary, compression", [(False, None), (True, None), (True, "zlib")]
)
def test_vtp_roundtrips_field_data(engine, binary, compression, tmp_path):
    mesh = _mesh()
    engine.vtp.write(tmp_path / "a.vtp", mesh, binary=binary, compression=compression)
    _same(engine.vtp.read(tmp_path / "a.vtp"), mesh)


def test_the_engines_read_each_others_field_data(tmp_path):
    if not hasattr(_core, "vtu_write_codec"):
        pytest.skip("no C++ core")
    mesh = _mesh()
    for fmt, pkg, twin in (
        ("vtu", meshioplusplus.vtu, _vtu),
        ("vtp", meshioplusplus.vtp, _vtp),
    ):
        set_strict_core(True)
        try:
            pkg.write(tmp_path / f"core.{fmt}", mesh)
        finally:
            set_strict_core(None)
        twin.write(tmp_path / f"py.{fmt}", mesh)
        for path in (tmp_path / f"core.{fmt}", tmp_path / f"py.{fmt}"):
            _same(twin.read(path), mesh)
            set_strict_core(True)
            try:
                _same(pkg.read(path), mesh)
            finally:
                set_strict_core(None)


def test_the_layout_is_what_vtk_reads(engine, tmp_path):
    engine.vtu.write(tmp_path / "a.vtu", _mesh(), binary=False, compression=None)
    text = (tmp_path / "a.vtu").read_text()
    # on the grid, before the piece, with the tuple count VTK requires
    assert text.index("<FieldData>") < text.index("<Piece ")
    assert 'Name="TimeValue" NumberOfTuples="1"' in text
    assert 'Name="params" NumberOfTuples="2" NumberOfComponents="3"' in text


def test_a_mesh_without_field_data_writes_no_field_data_element(engine, tmp_path):
    mesh = _mesh()
    mesh.field_data = {}
    engine.vtu.write(tmp_path / "a.vtu", mesh)
    engine.vtp.write(tmp_path / "a.vtp", mesh)
    assert "FieldData" not in (tmp_path / "a.vtu").read_text()
    assert "FieldData" not in (tmp_path / "a.vtp").read_text()


def test_writing_does_not_touch_the_callers_field_data(engine, tmp_path):
    mesh = _mesh()
    before = {k: (v.dtype, v.tolist()) for k, v in mesh.field_data.items()}
    engine.vtu.write(tmp_path / "a.vtu", mesh)
    engine.vtp.write(tmp_path / "a.vtp", mesh)
    assert before == {k: (v.dtype, v.tolist()) for k, v in mesh.field_data.items()}


def _vtu_text(grid_fd="", piece_fd=""):
    return (
        '<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" '
        'byte_order="LittleEndian">\n<UnstructuredGrid>\n'
        f'{grid_fd}<Piece NumberOfPoints="3" NumberOfCells="1">\n{piece_fd}'
        '<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">'
        "0 0 0 1 0 0 0 1 0</DataArray></Points>\n<Cells>"
        '<DataArray type="Int64" Name="connectivity" format="ascii">0 1 2</DataArray>'
        '<DataArray type="Int64" Name="offsets" format="ascii">3</DataArray>'
        '<DataArray type="UInt8" Name="types" format="ascii">5</DataArray>'
        "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n"
    )


def _fd(name, value, type_="Float64"):
    return (
        f'<FieldData><DataArray type="{type_}" Name="{name}" NumberOfTuples="1" '
        f'format="ascii">{value}</DataArray></FieldData>\n'
    )


def test_reads_vtks_time_value_on_the_grid_and_in_a_piece_the_piece_winning(
    engine, tmp_path
):
    (tmp_path / "g.vtu").write_text(_vtu_text(grid_fd=_fd("TimeValue", "0.25")))
    assert float(engine.vtu.read(tmp_path / "g.vtu").field_data["TimeValue"][0]) == 0.25
    (tmp_path / "b.vtu").write_text(
        _vtu_text(_fd("TimeValue", "0.25"), _fd("TimeValue", "0.75"))
    )
    assert float(engine.vtu.read(tmp_path / "b.vtu").field_data["TimeValue"][0]) == 0.75


def test_a_non_numeric_array_is_skipped_not_fatal(engine, tmp_path):
    """A `vtkStringArray` used to be ignored along with the whole section."""
    fd = (
        '<FieldData><DataArray type="String" Name="label" NumberOfTuples="1" '
        'format="ascii">bracket</DataArray>'
        '<DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" '
        'format="ascii">2</DataArray></FieldData>\n'
    )
    (tmp_path / "s.vtu").write_text(_vtu_text(grid_fd=fd))
    assert sorted(engine.vtu.read(tmp_path / "s.vtu").field_data) == ["TimeValue"]


def test_a_value_with_no_vtk_type_is_skipped_on_write_with_a_warning(
    engine, tmp_path, capsys
):
    mesh = _mesh()
    mesh.field_data["note"] = {"not": "an array"}
    if engine.vtu is _vtu:
        engine.vtu.write(tmp_path / "a.vtu", mesh)
        assert sorted(engine.vtu.read(tmp_path / "a.vtu").field_data) == sorted(
            k for k in mesh.field_data if k != "note"
        )
        assert "field_data 'note' is not a numeric array" in capsys.readouterr().err
    else:
        # the core refuses a non-array value and the shim falls back to the twin
        set_strict_core(None)
        meshioplusplus.vtu.write(tmp_path / "a.vtu", mesh)
        assert "note" not in meshioplusplus.vtu.read(tmp_path / "a.vtu").field_data


def test_metadata_names_field_data_like_a_real_read(tmp_path):
    mesh = _mesh()
    for name, pkg in (("a.vtu", meshioplusplus.vtu), ("a.vtp", meshioplusplus.vtp)):
        pkg.write(tmp_path / name, mesh)
        meta = meshioplusplus.read_metadata(tmp_path / name)
        assert meta["field_data_names"] == sorted(mesh.field_data)


def test_selective_reads_narrow_field_data(tmp_path):
    mesh = _mesh()
    meshioplusplus.vtu.write(tmp_path / "a.vtu", mesh)
    only = meshioplusplus.read(tmp_path / "a.vtu", arrays=["TimeValue"])
    assert sorted(only.field_data) == ["TimeValue"]
    assert meshioplusplus.read(tmp_path / "a.vtu", points_only=True).field_data == {}
