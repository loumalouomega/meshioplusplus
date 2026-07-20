"""Tests for the elementwise data-expression evaluator."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._data_calc import _calc_py

from .helpers_data import assert_same_geometry, data_mesh


def test_norm_equals_the_euclidean_magnitude():
    m = data_mesh()
    out = mp.data_calc(m, "norm(v)", location="point", output="speed")
    expected = np.linalg.norm(m.point_data["v"], axis=1)
    assert out.point_data["speed"] == pytest.approx(expected)
    # norm collapses to a scalar, so the result is 1-D.
    assert out.point_data["speed"].ndim == 1
    assert_same_geometry(m, out)


def test_linear_combination_matches_a_numpy_reference():
    m = data_mesh()
    out = mp.data_calc(m, "2 * T + 1", location="point", output="y")
    assert out.point_data["y"] == pytest.approx(2 * m.point_data["T"] + 1)


def test_difference_of_two_fields():
    m = data_mesh()
    m.point_data["p_old"] = np.arange(6.0)
    m.point_data["p_new"] = np.arange(6.0) * 3
    out = mp.data_calc(m, "p_new - p_old", location="point", output="dp")
    assert out.point_data["dp"] == pytest.approx(np.arange(6.0) * 2)


@pytest.mark.parametrize(
    "expr,expected",
    [
        ("1 + 2 * 3", 7.0),
        ("(1 + 2) * 3", 9.0),
        ("2 * -3", -6.0),
        ("-(2 + 3)", -5.0),
        ("10 / 4", 2.5),
        ("min(3, 5)", 3.0),
        ("max(3, 5)", 5.0),
        ("abs(0 - 4)", 4.0),
        ("sqrt(9)", 3.0),
    ],
)
def test_arithmetic(expr, expected):
    m = data_mesh()
    out = mp.data_calc(m, expr, location="point", output="o")
    assert out.point_data["o"][0] == pytest.approx(expected)


def test_unary_minus_on_an_array():
    m = data_mesh()
    out = mp.data_calc(m, "-T", location="point", output="o")
    assert out.point_data["o"] == pytest.approx(-m.point_data["T"])
    out2 = mp.data_calc(m, "--T", location="point", output="o")
    assert out2.point_data["o"] == pytest.approx(m.point_data["T"])


def test_scalar_broadcasts_against_a_vector():
    m = data_mesh()
    for expr in ("v * 2", "2 * v"):
        out = mp.data_calc(m, expr, location="point", output="o")
        assert out.point_data["o"].shape == (6, 3)
        assert out.point_data["o"] == pytest.approx(m.point_data["v"] * 2)


def test_result_shape_at_the_correct_location():
    m = data_mesh()
    out = mp.data_calc(m, "v - v", location="point", output="z")
    assert out.point_data["z"].shape == (6, 3)
    assert out.point_data["z"] == pytest.approx(np.zeros((6, 3)))


def test_colon_and_backtick_identifiers():
    m = data_mesh()
    m.point_data["gmsh:physical"] = np.arange(6.0)
    m.point_data["with space"] = np.ones(6)
    a = mp.data_calc(m, "gmsh:physical + 1", location="point", output="a")
    assert a.point_data["a"] == pytest.approx(np.arange(6.0) + 1)
    b = mp.data_calc(m, "`with space` * 3", location="point", output="b")
    assert b.point_data["b"] == pytest.approx(np.full(6, 3.0))


def test_cell_location_evaluates_per_block():
    m = data_mesh()
    out = mp.data_calc(m, "mat * 10", location="cell", output="mat2")
    assert len(out.cell_data["mat2"]) == len(m.cells)
    assert out.cell_data["mat2"][0] == pytest.approx([10.0, 20.0])
    assert out.cell_data["mat2"][1] == pytest.approx([30.0])


def test_field_location():
    m = data_mesh()
    out = mp.data_calc(m, "meta + 1", location="field", output="meta2")
    assert out.field_data["meta2"] == pytest.approx([2.0, 3.0, 4.0])


def test_division_by_zero_is_not_an_error():
    m = data_mesh()
    out = mp.data_calc(m, "1 / 0", location="point", output="o")
    assert np.isinf(out.point_data["o"][0])


def test_unknown_array_name_raises():
    m = data_mesh()
    with pytest.raises(ValueError) as exc:
        mp.data_calc(m, "temp + 1", location="point", output="o")
    msg = str(exc.value)
    assert "temp" in msg
    assert "available" in msg


def test_unknown_function_raises_and_lists_known_ones():
    m = data_mesh()
    with pytest.raises(ValueError) as exc:
        mp.data_calc(m, "log(T)", location="point", output="o")
    msg = str(exc.value)
    assert "log" in msg
    assert "sqrt" in msg


@pytest.mark.parametrize(
    "expr,needle",
    [
        ("norm(v, T)", "takes exactly 1 argument"),
        ("min(T)", "takes exactly 2 arguments"),
        ("T +", "unexpected end of expression"),
        ("(T", "expected ')'"),
        ("T T", "trailing input"),
        ("T # 1", "unexpected character '#'"),
        ("", "empty"),
    ],
)
def test_error_messages(expr, needle):
    m = data_mesh()
    with pytest.raises(ValueError) as exc:
        mp.data_calc(m, expr, location="point", output="o")
    assert needle in str(exc.value)


def test_component_width_mismatch_is_diagnosed():
    m = data_mesh()
    m.point_data["t9"] = np.ones((6, 9))
    with pytest.raises(ValueError) as exc:
        mp.data_calc(m, "v * t9", location="point", output="o")
    assert "cannot combine" in str(exc.value)


def test_depth_guard():
    m = data_mesh()
    expr = "(" * 200 + "T" + ")" * 200
    with pytest.raises(ValueError) as exc:
        mp.data_calc(m, expr, location="point", output="o")
    assert "nests deeper" in str(exc.value)


def test_existing_output_name_needs_overwrite():
    m = data_mesh()
    with pytest.raises(ValueError):
        mp.data_calc(m, "T * 2", location="point", output="T")
    out = mp.data_calc(m, "T * 2", location="point", output="T", overwrite=True)
    assert out.point_data["T"] == pytest.approx(m.point_data["T"] * 2)


def test_cpp_matches_python():
    pytest.importorskip("meshioplusplus._core")
    m = data_mesh()
    exprs = [
        "norm(v)",
        "2 * T + 1",
        "1 + 2 * 3",
        "-T",
        "--T",
        "min(T, 5)",
        "max(T, 5)",
        "abs(0 - T)",
        "sqrt(T)",
        "v * 2",
        "2 * v",
        "v - v",
        "(1 + 2) * 3",
        "T / 3",
        "norm(v * 2) + T",
    ]
    for expr in exprs:
        cpp = mp.data_calc(m, expr, location="point", output="o")
        py = _calc_py(m, expr, "point_data", "o", False)
        assert np.allclose(
            cpp.point_data["o"], py.point_data["o"], equal_nan=True
        ), expr


def test_roundtrip_write_read(tmp_path):
    m = data_mesh()
    out = mp.data_calc(m, "norm(v)", location="point", output="speed")
    path = tmp_path / "out.vtu"
    mp.write(path, out)
    back = mp.read(path)
    assert np.allclose(back.point_data["speed"], out.point_data["speed"])
