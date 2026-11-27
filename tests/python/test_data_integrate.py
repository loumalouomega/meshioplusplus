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


# --- mode prefixes: point:, flux:, boundary: --------------------------------


def _unit_cube():
    pts = np.array([[x, y, z] for z in (0, 1) for y in (0, 1) for x in (0, 1)], float)
    mesh = mp.Mesh(pts, [("hexahedron", np.array([[0, 1, 3, 2, 4, 5, 7, 6]]))])
    mesh.point_data["u"] = pts.copy()  # div u = 3
    mesh.point_data["T"] = pts[:, 0] + 2 * pts[:, 1] + 3 * pts[:, 2]
    return mesh


def _flat_two_triangles():
    pts = np.array([[0, 0, 0], [2, 0, 0], [2, 1, 0], [0, 1, 0.0]])
    mesh = mp.Mesh(pts, [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))])
    mesh.point_data["T"] = pts[:, 0] + 10 * pts[:, 1]
    mesh.point_data["w"] = np.tile([0.0, 0.0, 1.0], (4, 1))
    return mesh


@pytest.mark.parametrize("engine", ["core", "twin"])
def test_point_data_is_integrated_through_its_cell_means(engine):
    mesh = _flat_two_triangles()
    run = data_integrate if engine == "core" else _data_integrate_py
    r = run(mesh, ["point:T"])
    assert r[0]["name"] == "point:T"
    # exact for a linear field: the integral of x + 10 y over [0, 2] x [0, 1] is 12
    assert r[0]["domain"]["total_per_component"][0] == pytest.approx(12.0)
    assert r[0]["domain"]["domain_measure_per_component"][0] == pytest.approx(2.0)


@pytest.mark.parametrize("engine", ["core", "twin"])
def test_flux_through_a_surface_and_the_divergence_theorem(engine):
    run = data_integrate if engine == "core" else _data_integrate_py
    flat = run(_flat_two_triangles(), ["flux:w"])
    assert flat[0]["domain"]["total_per_component"][0] == pytest.approx(
        2.0
    )  # area * (w . +z)
    # the flux of u = (x, y, z) through the boundary of the unit cube equals the
    # integral of div u = 3 over its volume 1
    cube = run(_unit_cube(), ["boundary:flux:u"])
    assert cube[0]["domain"]["total_per_component"][0] == pytest.approx(3.0)
    # a constant field has zero flux through a closed surface
    c = _unit_cube()
    c.point_data["k"] = np.tile([1.0, -2.0, 3.0], (8, 1))
    zero = run(c, ["boundary:flux:k"])
    assert zero[0]["domain"]["total_per_component"][0] == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("engine", ["core", "twin"])
def test_boundary_point_integral_of_a_linear_field(engine):
    run = data_integrate if engine == "core" else _data_integrate_py
    r = run(_unit_cube(), ["boundary:point:T"])
    assert r[0]["domain"]["total_per_component"][0] == pytest.approx(18.0)  # six faces


@pytest.mark.parametrize("engine", ["core", "twin"])
def test_modes_keep_the_requested_order_and_mix_with_plain_names(engine):
    run = data_integrate if engine == "core" else _data_integrate_py
    mesh = _flat_two_triangles()
    mesh.cell_data["c"] = [np.array([1.0, 3.0])]
    r = run(mesh, ["point:T", "c", "flux:w"])
    assert [a["name"] for a in r] == ["point:T", "c", "flux:w"]
    assert r[1]["domain"]["total_per_component"][0] == pytest.approx(1.0 * 1 + 1.0 * 3)


@pytest.mark.parametrize("engine", ["core", "twin"])
def test_an_array_with_the_prefixed_name_wins(engine):
    run = data_integrate if engine == "core" else _data_integrate_py
    mesh = _flat_two_triangles()
    mesh.cell_data["point:T"] = [np.array([1.0, 1.0])]
    r = run(mesh, ["point:T"])
    assert r[0]["domain"]["total_per_component"][0] == pytest.approx(
        2.0
    )  # the cell array, not T


@pytest.mark.parametrize("engine", ["core", "twin"])
@pytest.mark.parametrize(
    "name,needle",
    [
        ("boundary:T", "needs 'point:'.*or 'flux:'"),
        ("point:nothing", "names no point_data array"),
        ("flux:T", "3-component"),
    ],
)
def test_mode_errors(engine, name, needle):
    run = data_integrate if engine == "core" else _data_integrate_py
    with pytest.raises(ValueError, match=needle):
        run(_flat_two_triangles(), [name])


def test_flux_needs_surface_cells():
    with pytest.raises(ValueError, match="surface"):
        data_integrate(_unit_cube(), ["flux:u"])
