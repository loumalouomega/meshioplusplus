"""Tests for the read-only per-array data summary."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._data_info import _info_py

from .helpers_data import data_mesh


def by_name(arrays, location, name):
    for a in arrays:
        if a["location"] == location and a["name"] == name:
            return a
    raise AssertionError(f"{location}:{name} not in the report")


def test_counts_every_array():
    arrays = mp.data_info(data_mesh())
    assert len(arrays) == 5
    locations = [a["location"] for a in arrays]
    assert locations.count("point_data") == 2
    assert locations.count("cell_data") == 2
    assert locations.count("field_data") == 1


def test_grouped_by_location_then_sorted_by_name():
    arrays = mp.data_info(data_mesh())
    assert [(a["location"], a["name"]) for a in arrays] == [
        ("point_data", "T"),
        ("point_data", "v"),
        ("cell_data", "mat"),
        ("cell_data", "tag"),
        ("field_data", "meta"),
    ]


def test_scalar_statistics_match_hand_computed_values():
    a = by_name(mp.data_info(data_mesh()), "point_data", "T")
    # T = {0, 1, 11, 10, 2, 12}
    assert a["num_entries"] == 6
    assert a["num_components"] == 1
    assert a["num_values"] == 6
    assert a["min"] == pytest.approx(0.0)
    assert a["max"] == pytest.approx(12.0)
    assert a["mean"] == pytest.approx(36.0 / 6.0)
    assert a["num_nan"] == 0
    assert a["num_inf"] == 0
    assert a["num_finite"] == 6


def test_vector_reports_shape_and_per_component_statistics():
    a = by_name(mp.data_info(data_mesh()), "point_data", "v")
    assert a["shape"] == (6, 3)
    assert a["num_components"] == 3
    assert a["num_entries"] == 6
    assert a["num_values"] == 18
    # Component 0 of v is {1,0,0,1,2,0}
    assert a["min_per_component"][0] == pytest.approx(0.0)
    assert a["max_per_component"][0] == pytest.approx(2.0)
    assert a["mean_per_component"][0] == pytest.approx(4.0 / 6.0)
    assert a["max_per_component"][2] == pytest.approx(1.0)


def test_cell_data_spans_every_block():
    a = by_name(mp.data_info(data_mesh()), "cell_data", "mat")
    assert a["num_blocks"] == 2
    assert a["num_entries"] == 3  # 2 triangles + 1 quad
    assert a["min"] == pytest.approx(1.0)
    assert a["max"] == pytest.approx(3.0)
    assert a["mean"] == pytest.approx(2.0)
    assert a["inconsistent_blocks"] is False


def test_dtype_reported_as_stored():
    arrays = mp.data_info(data_mesh())
    assert by_name(arrays, "point_data", "T")["dtype"] == "f8"
    assert by_name(arrays, "cell_data", "tag")["dtype"] == "i4"


def test_nan_and_inf_counts_are_correct():
    m = data_mesh()
    m.point_data["n"] = np.array([1.0, np.nan, np.inf, -np.inf, 2.0, np.nan])
    a = by_name(mp.data_info(m), "point_data", "n")
    assert a["num_nan"] == 2
    assert a["num_inf"] == 2
    assert a["num_finite"] == 2
    assert a["num_values"] == 6
    # Reductions ignore the non-finite values entirely.
    assert a["min"] == pytest.approx(1.0)
    assert a["max"] == pytest.approx(2.0)
    assert a["mean"] == pytest.approx(1.5)


def test_all_non_finite_array_yields_nan_statistics():
    m = data_mesh()
    m.point_data["n"] = np.full(6, np.nan)
    a = by_name(mp.data_info(m), "point_data", "n")
    assert a["num_finite"] == 0
    assert np.isnan(a["min"])
    assert np.isnan(a["max"])
    assert np.isnan(a["mean"])


def test_mesh_with_no_data():
    m = mp.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]]),
        [("triangle", np.array([[0, 1, 2]]))],
    )
    assert mp.data_info(m) == []


def test_does_not_modify_the_mesh():
    m = data_mesh()
    before = sorted(m.point_data)
    mp.data_info(m)
    assert sorted(m.point_data) == before


def test_cpp_matches_python():
    pytest.importorskip("meshioplusplus._core")
    m = data_mesh()
    m.point_data["n"] = np.array([1.0, np.nan, np.inf, -np.inf, 2.0, np.nan])
    cpp = mp.data_info(m)
    py = _info_py(m)
    assert len(cpp) == len(py)
    for a, b in zip(cpp, py):
        assert a["location"] == b["location"]
        assert a["name"] == b["name"]
        for key in (
            "num_blocks",
            "num_entries",
            "num_components",
            "num_values",
            "num_nan",
            "num_inf",
            "num_finite",
            "inconsistent_blocks",
        ):
            assert a[key] == b[key], (a["name"], key)
        for key in ("min", "max", "mean"):
            assert np.isclose(a[key], b[key], equal_nan=True), (a["name"], key)
