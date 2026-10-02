"""Sets/data semantics, Python mutation contract, and direct native parity."""

from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_array_equal

import meshioplusplus as mp
from meshioplusplus._pipeline import _apply_sets_data


def _mesh():
    return mp.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        [("triangle", [[0, 1, 2], [0, 2, 3]]), ("line", [[0, 1]])],
        point_sets={"z": [0, 1], "a": [1, 2]},
        cell_sets={"z": [[0], [0]], "a": [[1], [0]]},
        field_data={"text": "not a native numeric field"},
        info={"keep": "this"},
    )


def _core():
    core = getattr(mp, "_core", None)
    if core is None or not hasattr(core, "sets_to_data"):
        pytest.skip("requires rebuilt native sets/data kernels")
    return core


@pytest.mark.parametrize("location", ["point", "cell"])
def test_direct_native_matches_python_order_and_preserves_input(location, monkeypatch):
    core = _core()
    mesh = _mesh()
    before = mesh.copy()
    sets = getattr(mesh, f"{location}_sets")
    native = core.sets_to_data(mesh, location, order=list(sets))
    # Force the reference twin independently of native dispatch.
    monkeypatch.setattr(mp, "_core", None)
    getattr(mesh, f"{location}_sets_to_data")()
    key = "z-a"
    expected = getattr(mesh, f"{location}_data")[key]
    actual = getattr(native, f"{location}_data")[key]
    if location == "point":
        assert_array_equal(actual, expected)
    else:
        for a, b in zip(actual, expected):
            assert_array_equal(a, b)
    assert not getattr(native, f"{location}_sets")
    assert_array_equal(native.points, before.points)
    assert getattr(before, f"{location}_sets")
    assert_array_equal(before.cells[1].data, native.cells[1].data)


def test_mesh_native_dispatch_preserves_python_only_metadata(monkeypatch):
    core = _core()
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    mesh = _mesh()
    points, cells = mesh.points, mesh.cells
    info, fields = mesh.info, mesh.field_data
    mesh.cell_sets_to_data("labels")
    mesh.point_sets_to_data(":")
    assert mesh.points is points and mesh.cells is cells
    assert mesh.info is info and mesh.field_data is fields
    assert_array_equal(mesh.cell_data["labels"][0], [0, 1])
    assert_array_equal(mesh.cell_data["labels"][1], [1])
    assert_array_equal(mesh.point_data["z:a"], [0, 1, 1, -1])
    assert core is not None


@pytest.mark.parametrize("dtype", [np.int8, np.uint16, np.int64, np.uint64])
def test_direct_integer_data_to_sets(dtype, monkeypatch):
    core = _core()
    mesh = _mesh()
    mesh.point_sets = {}
    values = [0, 1, 0, 1]
    if dtype == np.uint64:
        values = [0, 2**64 - 1, 0, 2**53 + 1]
    elif dtype == np.int64:
        values = [-(2**63), 2**63 - 1, -(2**63), 0]
    mesh.point_data["ids"] = np.array(values, dtype=dtype)
    native = core.data_to_sets(mesh, "point", "ids")
    monkeypatch.setattr(mp, "_core", None)
    reference = _apply_sets_data(
        mesh, {"Op": "DataToSets", "Location": "point", "Key": "ids"}
    )
    assert list(native.point_sets) == sorted(reference.point_sets)
    for name in native.point_sets:
        assert_array_equal(native.point_sets[name], reference.point_sets[name])
    assert "ids" not in native.point_data
    assert "ids" in mesh.point_data


def test_pipeline_sets_data_order_names_and_validation():
    mesh = _mesh()
    step = {
        "Op": "SetsToData",
        "Location": "cell",
        "Order": ["a", "z"],
        "Name": "labels",
    }
    result = _apply_sets_data(mesh, step)
    assert_array_equal(result.cell_data["labels"][0], [1, 0])
    assert_array_equal(result.cell_data["labels"][1], [1])
    assert mesh.cell_sets
    assert result.info == mesh.info
    assert result.field_data == mesh.field_data
    with pytest.raises(ValueError, match="Order"):
        _apply_sets_data(mesh, dict(step, Order=["a", "a"]))
    with pytest.raises(ValueError, match="Location"):
        _apply_sets_data(mesh, dict(step, Location="field"))
    with pytest.raises(ValueError, match="requires 'Key'"):
        _apply_sets_data(mesh, {"Op": "DataToSets"})


def test_method_dispatch_passes_order_and_keeps_python_metadata(monkeypatch, capsys):
    calls = []

    def sets_to_data(mesh, location, **kwargs):
        calls.append((location, kwargs))
        result = mesh.copy()
        result.point_data["z:a"] = np.array([0, 1, 1, -1])
        result.point_sets = {}
        result.field_data = {}  # Native conversion cannot represent the string.
        return result

    monkeypatch.setattr(mp, "_core", SimpleNamespace(sets_to_data=sets_to_data))
    mesh = _mesh()
    points, cells, fields, info = mesh.points, mesh.cells, mesh.field_data, mesh.info
    mesh.point_sets_to_data(":")
    assert calls == [("point", {"join": ":", "order": ["z", "a"]})]
    assert mesh.points is points and mesh.cells is cells
    assert mesh.field_data is fields and mesh.info is info
    assert not mesh.point_sets
    assert "Not all points" in capsys.readouterr().err


def test_native_argument_errors_do_not_fall_back(monkeypatch):
    def data_to_sets(*args, **kwargs):
        raise ValueError("bad scalar shape")

    monkeypatch.setattr(mp, "_core", SimpleNamespace(data_to_sets=data_to_sets))
    mesh = _mesh()
    mesh.point_data["ids"] = np.array([0, 1, 0, 1])
    with pytest.raises(ValueError, match="bad scalar shape"):
        mesh.point_data_to_sets("ids")
    assert "ids" in mesh.point_data


def test_pipeline_sets_data_roundtrip_and_reports(tmp_path):
    mesh = _mesh()
    # The string field is intentionally Python-only; VTK writes numeric fields.
    mesh.field_data = {}
    source, output = tmp_path / "source.vtu", tmp_path / "output.vtu"
    mp.write(source, mesh)
    # VTU pads two-dimensional coordinates to three dimensions on disk.
    source_points = mp.read(source).points.copy()
    report = mp.run_pipeline(
        {
            "Input": {"Path": str(source)},
            "Operations": [
                {"Op": "SetsToData", "Location": "cell", "Order": ["z", "a"]},
                {"Op": "DataToSets", "Location": "cell", "Key": "z-a"},
            ],
            "Output": {"Path": str(output)},
        }
    )
    assert [step["op"] for step in report["steps"]] == ["SetsToData", "DataToSets"]
    result = mp.read(output)
    assert "z-a" not in result.cell_data
    assert_array_equal(result.cell_sets["z"][0], [0])
    assert_array_equal(result.cell_sets["a"][0], [1])
    assert_array_equal(result.cell_sets["a"][1], [0])
    assert_array_equal(result.points, source_points)
