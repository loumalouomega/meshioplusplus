"""Tests for the data array-management operation (rename / drop / keep)."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._data_manage import _manage_py

from .helpers_data import assert_same_geometry, data_mesh


def test_rename_changes_the_key_and_preserves_values():
    m = data_mesh()
    out = mp.data_rename(m, "point", "T", "temperature")
    assert "T" not in out.point_data
    assert "temperature" in out.point_data
    assert np.array_equal(out.point_data["temperature"], m.point_data["T"])
    # Untouched arrays survive.
    assert "v" in out.point_data
    assert "mat" in out.cell_data
    assert "meta" in out.field_data
    assert_same_geometry(m, out)


def test_rename_multi_block_cell_data():
    m = data_mesh()
    out = mp.data_rename(m, "cell", "mat", "material")
    assert "mat" not in out.cell_data
    # One array per cell block must survive.
    assert len(out.cell_data["material"]) == len(m.cells)
    assert out.cell_data["material"][0][0] == 1.0
    assert out.cell_data["material"][1][0] == 3.0


def test_drop_removes_exactly_the_named_arrays():
    m = data_mesh()
    out = mp.data_drop(m, "point", ["T"])
    assert "T" not in out.point_data
    assert "v" in out.point_data
    # Other locations untouched.
    assert set(out.cell_data) == {"mat", "tag"}
    assert "meta" in out.field_data
    assert_same_geometry(m, out)


def test_drop_multiple_names_and_locations():
    m = data_mesh()
    out = mp.data_manage(m, drop=[("point", "T"), ("point", "v"), ("cell", "tag")])[
        "mesh"
    ]
    assert out.point_data == {}
    assert set(out.cell_data) == {"mat"}


def test_keep_retains_only_the_named_subset():
    m = data_mesh()
    out = mp.data_keep(m, "cell", ["tag"])
    assert set(out.cell_data) == {"tag"}
    # A location the whitelist does not mention is untouched.
    assert set(out.point_data) == {"T", "v"}
    assert "meta" in out.field_data


def test_keep_nothing_drops_everything_at_that_location():
    m = data_mesh()
    out = mp.data_keep(m, "point", [])
    assert out.point_data == {}
    assert "mat" in out.cell_data


def test_unknown_key_raises_a_clear_error():
    m = data_mesh()
    with pytest.raises(ValueError) as exc:
        mp.data_drop(m, "point", ["nope"])
    msg = str(exc.value)
    assert "nope" in msg
    assert "point_data" in msg
    # The message lists every available key.
    assert "T" in msg
    assert "v" in msg


def test_ignore_missing_silences_unknown_keys():
    m = data_mesh()
    out = mp.data_drop(m, "point", ["nope"], ignore_missing=True)
    assert "T" in out.point_data


def test_rename_onto_an_existing_name_raises():
    m = data_mesh()
    with pytest.raises(ValueError):
        mp.data_rename(m, "point", "T", "v")


def test_two_renames_to_the_same_target_raise():
    m = data_mesh()
    with pytest.raises(ValueError):
        mp.data_manage(m, rename=[("point", "T", "x"), ("point", "v", "x")])


def test_swapping_two_names_is_allowed():
    m = data_mesh()
    out = mp.data_manage(m, rename=[("point", "T", "v"), ("point", "v", "T")])["mesh"]
    # The old "v" (a 3-vector) is now called "T" and vice versa.
    assert out.point_data["T"].shape == m.point_data["v"].shape
    assert out.point_data["v"].shape == m.point_data["T"].shape


def test_report_lists_what_was_dropped_and_renamed():
    m = data_mesh()
    r = mp.data_manage(m, drop=[("point", "T")], rename=[("field", "meta", "metadata")])
    assert r["dropped"] == ["point_data:T"]
    assert r["renamed"] == [("field_data:meta", "field_data:metadata")]


def test_unknown_location_raises():
    m = data_mesh()
    with pytest.raises(ValueError):
        mp.data_drop(m, "vertex", ["T"])


def test_sets_survive_unchanged():
    m = data_mesh()
    m.point_sets = {"left": np.array([0, 3])}
    m.cell_sets = {"tris": [np.array([0, 1]), np.array([], dtype=int)]}
    out = mp.data_drop(m, "point", ["T"])
    assert np.array_equal(out.point_sets["left"], np.array([0, 3]))
    assert np.array_equal(out.cell_sets["tris"][0], np.array([0, 1]))


def test_cpp_matches_python():
    core = pytest.importorskip("meshioplusplus._core")
    m = data_mesh()
    cases = [
        ([], [("point", "T")], []),
        ([("cell", "mat")], [], []),
        ([], [], [("point", "T", "temp")]),
    ]
    for keep, drop, rename in cases:
        cpp = core.data_manage(m, keep, drop, rename, False)
        py = _manage_py(m, keep, drop, rename, False)
        assert cpp["dropped"] == py["dropped"]
        assert [tuple(t) for t in cpp["renamed"]] == py["renamed"]
        assert sorted(cpp["mesh"].point_data) == sorted(py["mesh"].point_data)
        assert sorted(cpp["mesh"].cell_data) == sorted(py["mesh"].cell_data)


def test_roundtrip_write_read(tmp_path):
    m = data_mesh()
    out = mp.data_rename(m, "point", "T", "temperature")
    path = tmp_path / "out.vtu"
    mp.write(path, out)
    back = mp.read(path)
    assert "temperature" in back.point_data
    assert np.allclose(back.point_data["temperature"], m.point_data["T"])
