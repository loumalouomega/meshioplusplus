"""Direct native MED multi-mesh/profile coverage, without fallback masking."""

import ctypes
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

import meshioplusplus as mp
from meshioplusplus.med._med import _ensure_med_families
from meshioplusplus.med._medmulti import read_med_multi as reference_read
from meshioplusplus.med._medmulti import write_med_multi as reference_write

h5py = pytest.importorskip("h5py")


def _core():
    core = mp._core
    if not getattr(core, "__has_hdf5__", False) or not hasattr(core, "med_read_named"):
        pytest.skip("requires rebuilt native MED named-mesh API")
    return core


def _mesh(value=1.0):
    return mp.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        [("triangle", [[0, 1, 2], [0, 2, 3]])],
        point_data={"pressure": np.full(4, value)},
        cell_data={"stress": [np.full((2, 2), value)]},
        point_sets={"rim": [0, 2]},
        cell_sets={"domain": [[0, 1]]},
    )


@pytest.mark.parametrize("native_writer", [False, True])
def test_native_multi_cross_reference_and_selection(
    tmp_path, monkeypatch, native_writer
):
    core = _core()
    path = tmp_path / "multi.med"
    inputs = [_mesh(2), _mesh(7)]
    names = ["z_surface", "a_volume"]
    if native_writer:
        core.med_write_multi(str(path), inputs, names, "4.1.0")
    else:
        # The historical multi writer consumes family arrays, not bare sets.
        reference_write(path, [_ensure_med_families(mesh) for mesh in inputs], names)
    assert core.med_mesh_names(str(path)) == names
    expected, ref_names = reference_read(path)
    assert ref_names == names
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    actual, actual_names = mp.med.read_med_multi(path)
    assert actual_names == names
    for index, name in enumerate(names):
        selected = mp.med.read(path, mesh_name=name)
        for mesh in (selected, actual[index]):
            assert_array_equal(mesh.points, expected[index].points)
            assert_array_equal(mesh.cells[0].data, expected[index].cells[0].data)
            assert_allclose(
                mesh.point_data["pressure"].reshape(-1),
                expected[index].point_data["pressure"].reshape(-1),
            )
            assert_allclose(
                mesh.cell_data["stress"][0], expected[index].cell_data["stress"][0]
            )
            assert_array_equal(mesh.point_sets["rim"], [0, 2])
            assert_array_equal(mesh.cell_sets["domain"][0], [0, 1])
    with pytest.raises(Exception, match="no mesh named"):
        core.med_read_named(str(path), "missing")


def _profile_file(path, location):
    reference_write(path, [_mesh()], ["mesh"])
    field, support_name, indices = (
        ("pressure", "NOE", [4, 2])
        if location == "point"
        else ("stress", "MAI.TR3", [2])
    )
    with h5py.File(path, "r+") as file:
        field_group = file[f"CHA/{field}"]
        step = field_group[next(iter(field_group))]
        support = step[support_name]
        support.move("MED_NO_PROFILE_INTERNAL", "subset")
        support.attrs.modify("PFL", np.bytes_("subset"))
        group = support["subset"]
        del group["CO"]
        values = np.array([9.0, 5.0])
        group.create_dataset("CO", data=values)
        definition = file.require_group("PROFILS").create_group("subset")
        definition.attrs["NBR"] = len(indices)
        definition.create_dataset("PFL", data=np.array(indices, dtype=np.int64))


@pytest.mark.parametrize("location", ["point", "cell"])
def test_native_profile_expansion_matches_reference(tmp_path, monkeypatch, location):
    core = _core()
    path = tmp_path / "profile.med"
    _profile_file(path, location)
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    mesh = core.med_read_named(str(path), "mesh")
    expected = reference_read(path)[0][0]
    if location == "point":
        actual = mesh.point_data["pressure"].reshape(-1)
        assert_allclose(actual, [np.nan, 5.0, np.nan, 9.0], equal_nan=True)
        assert_allclose(
            actual, expected.point_data["pressure"].reshape(-1), equal_nan=True
        )
    else:
        actual = mesh.cell_data["stress"][0]
        assert_allclose(actual, [[np.nan, np.nan], [9.0, 5.0]], equal_nan=True)
        assert_allclose(actual, expected.cell_data["stress"][0], equal_nan=True)


@pytest.mark.parametrize("bad_index", [0, -1, 5])
def test_native_profile_rejects_invalid_indices(tmp_path, bad_index):
    core = _core()
    path = tmp_path / "bad-profile.med"
    _profile_file(path, "point")
    with h5py.File(path, "r+") as file:
        file["PROFILS/subset/PFL"][0] = bad_index
    with pytest.raises(Exception, match="profile index is out of range"):
        core.med_read_named(str(path), "mesh")


def test_named_mesh_reference_fallback_return_contract(tmp_path, monkeypatch):
    path = tmp_path / "multi.med"
    reference_write(path, [_mesh(2), _mesh(7)], ["first", "second"])
    monkeypatch.setattr(mp.med, "_HAS_HDF5", False)
    meshes, names = mp.med.read_med_multi(path)
    assert names == ["first", "second"]
    selected = mp.med.read(path, mesh_name="second")
    assert_allclose(selected.point_data["pressure"], meshes[1].point_data["pressure"])
    with pytest.raises(mp.ReadError, match="no mesh named"):
        mp.med.read(path, mesh_name="missing")


def test_med_named_python_cli(tmp_path, monkeypatch):
    from meshioplusplus._cli import main

    source, target = tmp_path / "multi.med", tmp_path / "selected.vtu"
    reference_write(source, [_mesh(2), _mesh(7)], ["fluid", "solid"])
    monkeypatch.setattr(mp.med, "_HAS_HDF5", False)
    main(["convert", str(source), str(target), "--mesh-name=solid"])
    assert_allclose(mp.read(target).point_data["pressure"].reshape(-1), [7, 7, 7, 7])


@pytest.mark.parametrize("native", [False, True])
def test_empty_component_name_metadata_does_not_emit_invalid_vtu(tmp_path, native):
    from meshioplusplus.vtu._vtu import write as reference_vtu_write

    mesh = _mesh()
    mesh.field_data["med:nom"] = [[], []]
    path = tmp_path / "empty-nom.vtu"
    if native:
        core = _core()
        core.vtu_write(str(path), mesh, False, False)
    else:
        reference_vtu_write(path, mesh)
    assert_allclose(mp.read(path).point_data["pressure"].reshape(-1), [1, 1, 1, 1])
    assert "med:nom" not in mp.read(path).field_data


def test_med_named_native_cli(tmp_path):
    _core()
    configured = os.environ.get("MESHIOPLUSPLUS_NATIVE_CLI")
    binary = (
        Path(configured)
        if configured
        else Path(__file__).resolve().parents[2] / "build/cpp-release/meshioplusplus"
    )
    if not binary.is_file():
        pytest.skip("native CLI is not built")
    source, target = tmp_path / "multi.med", tmp_path / "selected.vtu"
    reference_write(source, [_mesh(2), _mesh(7)], ["fluid", "solid"])
    result = subprocess.run(
        [str(binary), "convert", str(source), str(target), "--mesh-name=solid"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert_allclose(mp.read(target).point_data["pressure"].reshape(-1), [7, 7, 7, 7])
    result = subprocess.run(
        [str(binary), "convert", str(source), str(target), "--mesh-name=missing"],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "no mesh named" in result.stderr


def _med_c_library():
    configured = os.environ.get("MESHIOPLUSPLUS_LIB")
    library = (
        Path(configured)
        if configured
        else Path(__file__).resolve().parents[2]
        / "build/cpp-release/libmeshioplusplus.so"
    )
    if not library.is_file():
        pytest.skip("C API shared library is not built")
    lib = ctypes.CDLL(str(library))
    if not hasattr(lib, "mio_med_read_named"):
        pytest.skip("C API MED named-mesh functions are not built")
    lib.mio_format_readable.argtypes = [ctypes.c_char_p]
    lib.mio_format_readable.restype = ctypes.c_int
    if not lib.mio_format_readable(b"med"):
        pytest.skip("C API library was built without HDF5")
    lib.mio_med_mesh_count.argtypes = [ctypes.c_char_p]
    lib.mio_med_mesh_count.restype = ctypes.c_int64
    lib.mio_med_mesh_name.argtypes = [
        ctypes.c_char_p,
        ctypes.c_int64,
        ctypes.c_void_p,
        ctypes.c_int64,
    ]
    lib.mio_med_mesh_name.restype = ctypes.c_int64
    lib.mio_med_read_named.argtypes = [
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_void_p,
    ]
    lib.mio_med_read_named.restype = ctypes.c_void_p
    lib.mio_med_write_multi.argtypes = [
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.POINTER(ctypes.c_char_p),
        ctypes.c_int64,
        ctypes.c_char_p,
    ]
    lib.mio_med_write_multi.restype = ctypes.c_int
    lib.mio_mesh_free.argtypes = [ctypes.c_void_p]
    lib.mio_mesh_num_points.argtypes = [ctypes.c_void_p]
    lib.mio_mesh_num_points.restype = ctypes.c_int64
    lib.mio_mesh_get_point_data.argtypes = [
        ctypes.c_void_p,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int64),
    ]
    lib.mio_mesh_get_point_data.restype = ctypes.c_int
    return lib


def _med_c_point_values(lib, handle, name, count):
    data = ctypes.c_void_p()
    dtype = ctypes.c_int()
    assert (
        lib.mio_mesh_get_point_data(
            handle, name, ctypes.byref(data), ctypes.byref(dtype), None, None
        )
        == 0
    )
    types = [
        ctypes.c_float,
        ctypes.c_double,
        ctypes.c_int8,
        ctypes.c_int16,
        ctypes.c_int32,
        ctypes.c_int64,
        ctypes.c_uint8,
        ctypes.c_uint16,
        ctypes.c_uint32,
        ctypes.c_uint64,
    ]
    values = ctypes.cast(data, ctypes.POINTER(types[dtype.value]))
    return np.ctypeslib.as_array(values, shape=(count,)).copy()


def test_med_c_api_names_buffers_selection_and_multi_write(tmp_path):
    lib = _med_c_library()
    source, target = tmp_path / "multi.med", tmp_path / "rewritten.med"
    reference_write(source, [_mesh(2), _mesh(7)], ["z_mesh", "a_mesh"])
    path = os.fsencode(source)
    assert lib.mio_med_mesh_count(path) == 2
    assert lib.mio_med_mesh_count(None) == -1
    assert lib.mio_med_mesh_name(path, 0, None, 0) == len("z_mesh")
    buffer = ctypes.create_string_buffer(3)
    assert lib.mio_med_mesh_name(path, 0, buffer, len(buffer)) == len("z_mesh")
    assert buffer.value == b"z_"
    assert lib.mio_med_mesh_name(path, 2, None, 0) == -1
    assert lib.mio_med_read_named(path, b"missing", None) is None
    handles = []
    try:
        for name, expected in [(b"z_mesh", 2), (b"a_mesh", 7)]:
            handle = lib.mio_med_read_named(path, name, None)
            assert handle is not None
            handles.append(handle)
            assert lib.mio_mesh_num_points(handle) == 4
            assert_allclose(
                _med_c_point_values(lib, handle, b"pressure", 4), np.full(4, expected)
            )
        inputs = (ctypes.c_void_p * 2)(*handles)
        names = (ctypes.c_char_p * 2)(b"first", b"second")
        assert lib.mio_med_write_multi(os.fsencode(target), inputs, names, 2, None) == 0
        assert lib.mio_med_mesh_count(os.fsencode(target)) == 2
        assert lib.mio_med_write_multi(None, inputs, names, 2, None) != 0
    finally:
        for handle in handles:
            lib.mio_mesh_free(handle)


def test_med_c_api_profile_expansion_and_bounds(tmp_path):
    lib = _med_c_library()
    path = tmp_path / "profile.med"
    _profile_file(path, "point")
    handle = lib.mio_med_read_named(os.fsencode(path), b"mesh", None)
    assert handle is not None
    try:
        assert_allclose(
            _med_c_point_values(lib, handle, b"pressure", 4),
            [np.nan, 5.0, np.nan, 9.0],
            equal_nan=True,
        )
    finally:
        lib.mio_mesh_free(handle)
    with h5py.File(path, "r+") as file:
        file["PROFILS/subset/PFL"][0] = 0
    assert lib.mio_med_read_named(os.fsencode(path), b"mesh", None) is None
