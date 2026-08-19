"""Tests for cell-measure-weighted field integration (the ``data_integrate``
operation)."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import data_integrate
from meshioplusplus._data_integrate import _data_integrate_py
from meshioplusplus._regions import Region


def _quad_row(n, shift=0.0):
    """A 1 x n row of unit-area quads in the z = 0 plane."""
    pts = []
    for j in range(2):
        for i in range(n + 1):
            pts.append([i + shift, j + shift, 0.0])
    pts = np.array(pts, float)

    def pid(i, j):
        return j * (n + 1) + i

    cells = np.array(
        [[pid(i, 0), pid(i + 1, 0), pid(i + 1, 1), pid(i, 1)] for i in range(n)]
    )
    return mp.Mesh(pts, [("quad", cells)])


def _domain(report, name="f"):
    for arr in report:
        if arr["name"] == name:
            return arr["domain"]
    raise KeyError(name)


def _region(report, region_name, name="f"):
    for arr in report:
        if arr["name"] == name:
            for r in arr["regions"]:
                if r["name"] == region_name:
                    return r
    raise KeyError(region_name)


def test_uniform_field_total_equals_value_times_domain_measure():
    mesh = _quad_row(4)
    mesh.cell_data["f"] = [np.array([3.0, 3.0, 3.0, 3.0])]

    report = data_integrate(mesh)
    domain = _domain(report)
    assert domain["domain_measure_per_component"] == pytest.approx([4.0])
    assert domain["total_per_component"] == pytest.approx([12.0])


def test_mean_recovers_the_constant():
    mesh = _quad_row(3)
    mesh.cell_data["f"] = [np.array([5.0, 5.0, 5.0])]

    domain = _domain(data_integrate(mesh))
    assert domain["mean_per_component"] == pytest.approx([5.0])


def test_translating_the_mesh_does_not_change_the_total():
    vals = [1.0, 2.0, 3.0, 4.0, 5.0]
    a = _quad_row(5)
    b = _quad_row(5, shift=1.0e8)
    a.cell_data["f"] = [np.array(vals)]
    b.cell_data["f"] = [np.array(vals)]

    total_a = _domain(data_integrate(a))["total_per_component"][0]
    total_b = _domain(data_integrate(b))["total_per_component"][0]
    assert total_a == pytest.approx(total_b, abs=1e-6)


def test_ragged_or_unsupported_cells_are_excluded_not_zeroed():
    mesh = _quad_row(3)
    mesh.cells.append(mp.CellBlock("vertex", np.array([[0], [1]])))
    mesh.cell_data["f"] = [np.array([1.0, 1.0, 1.0]), np.array([99.0, 99.0])]

    domain = _domain(data_integrate(mesh))
    assert domain["num_cells"] == 3
    assert domain["num_skipped"] == 2
    assert domain["total_per_component"] == pytest.approx([3.0])


def test_non_finite_values_are_excluded_from_numerator_and_denominator():
    mesh = _quad_row(3)
    mesh.cell_data["f"] = [np.array([1.0, np.nan, 3.0])]

    domain = _domain(data_integrate(mesh))
    assert domain["num_nan_per_component"] == [1]
    assert domain["total_per_component"] == pytest.approx([4.0])
    assert domain["domain_measure_per_component"] == pytest.approx([2.0])
    assert domain["mean_per_component"] == pytest.approx([2.0])


def test_per_region_totals_are_independent_not_a_partition():
    mesh = _quad_row(4)
    mesh.cell_data["f"] = [np.array([1.0, 2.0, 3.0, 4.0])]
    mesh.regions.append(
        Region(name="a", kind="cell", entries=np.array([0, 1], dtype=np.int64))
    )
    mesh.regions.append(
        Region(name="b", kind="cell", entries=np.array([1, 2], dtype=np.int64))
    )

    report = data_integrate(mesh)
    assert _region(report, "a")["total_per_component"] == pytest.approx([3.0])
    assert _region(report, "b")["total_per_component"] == pytest.approx([5.0])
    assert _domain(report)["total_per_component"] == pytest.approx([10.0])


def test_point_and_side_regions_are_skipped():
    mesh = _quad_row(2)
    mesh.cell_data["f"] = [np.array([1.0, 1.0])]
    mesh.regions.append(
        Region(name="pts", kind="point", entries=np.array([0, 1], dtype=np.int64))
    )

    report = data_integrate(mesh)
    assert report[0]["regions"] == []


def test_empty_arrays_means_every_cell_data_array():
    mesh = _quad_row(2)
    mesh.cell_data["f"] = [np.array([1.0, 1.0])]
    mesh.cell_data["g"] = [np.array([2.0, 2.0])]

    report = data_integrate(mesh)
    assert {a["name"] for a in report} == {"f", "g"}


def test_point_data_name_raises_and_names_the_fix():
    mesh = _quad_row(2)
    mesh.point_data["f"] = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])

    with pytest.raises(ValueError, match="point_data_to_cell_data"):
        data_integrate(mesh, arrays=["f"])


def test_unknown_array_raises_and_lists_what_exists():
    mesh = _quad_row(2)
    mesh.cell_data["f"] = [np.array([1.0, 1.0])]

    with pytest.raises(ValueError, match="f"):
        data_integrate(mesh, arrays=["nope"])


def test_data_integrate_is_exported():
    assert "data_integrate" in mp.__all__
    assert mp.data_integrate is data_integrate


# --- C++ / numpy parity -------------------------------------------------------


@pytest.mark.parametrize(
    "with_regions",
    [False, True],
)
def test_cpp_matches_python(with_regions):
    """Pins parity only -- both sides can be identically wrong -- so it never
    substitutes for the analytic tests above."""
    pytest.importorskip("meshioplusplus._core")
    mesh = _quad_row(6)
    mesh.cell_data["f"] = [np.array([1.0 + 0.37 * i for i in range(6)])]
    mesh.cell_data["v"] = [np.array([[1.0 + i, 2.0 - i] for i in range(6)])]
    if with_regions:
        mesh.regions.append(
            Region(name="a", kind="cell", entries=np.array([0, 1, 4], dtype=np.int64))
        )
        mesh.regions.append(
            Region(name="b", kind="cell", entries=np.array([1, 2, 3], dtype=np.int64))
        )

    cpp_report = mp.data_integrate(mesh)
    py_report = _data_integrate_py(mesh, [])

    assert [a["name"] for a in cpp_report] == [a["name"] for a in py_report]
    for cpp_arr, py_arr in zip(cpp_report, py_report):
        assert cpp_arr["num_components"] == py_arr["num_components"]
        for key in (
            "domain_measure_per_component",
            "total_per_component",
            "mean_per_component",
        ):
            assert cpp_arr["domain"][key] == pytest.approx(py_arr["domain"][key])
        assert (
            cpp_arr["domain"]["num_nan_per_component"]
            == py_arr["domain"]["num_nan_per_component"]
        )
        assert [r["name"] for r in cpp_arr["regions"]] == [
            r["name"] for r in py_arr["regions"]
        ]
        for cpp_r, py_r in zip(cpp_arr["regions"], py_arr["regions"]):
            for key in (
                "domain_measure_per_component",
                "total_per_component",
                "mean_per_component",
            ):
                assert cpp_r[key] == pytest.approx(py_r[key])


def test_cpp_matches_python_on_ragged_and_nan_cells():
    """The parametrized parity test above never exercises a skipped cell or a
    non-finite value."""
    pytest.importorskip("meshioplusplus._core")
    mesh = _quad_row(3)
    mesh.cells.append(mp.CellBlock("vertex", np.array([[0], [1]])))
    mesh.cell_data["f"] = [np.array([1.0, np.nan, 3.0]), np.array([99.0, 99.0])]

    cpp_report = mp.data_integrate(mesh)
    py_report = _data_integrate_py(mesh, [])
    assert cpp_report == py_report
