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
        mp.data_calc(m, "foo(T)", location="point", output="o")
    msg = str(exc.value)
    assert "foo" in msg
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


_MATH = {
    "exp": np.exp,
    "log": np.log,
    "log10": np.log10,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "atan": np.arctan,
    "sinh": np.sinh,
    "cosh": np.cosh,
    "tanh": np.tanh,
    "floor": np.floor,
    "ceil": np.ceil,
    "sign": np.sign,
}


@pytest.mark.parametrize("name", sorted(_MATH))
def test_unary_math_functions_match_numpy(name):
    m = data_mesh()
    # T - 5.5 spans negative, zero-adjacent and positive values; exp stays finite.
    out = mp.data_calc(m, f"{name}((T - 5.5) / 4)", location="point", output="o")
    t = (m.point_data["T"] - 5.5) / 4
    with np.errstate(all="ignore"):
        expect = _MATH[name](t)
    np.testing.assert_allclose(out.point_data["o"], expect, rtol=1e-12, equal_nan=True)


def test_inverse_trigonometric_functions_and_their_domain():
    m = data_mesh()
    out = mp.data_calc(m, "asin(T / 12)", location="point", output="o")
    np.testing.assert_allclose(out.point_data["o"], np.arcsin(m.point_data["T"] / 12))
    out = mp.data_calc(m, "acos(T / 6)", location="point", output="o")
    got = out.point_data["o"]
    # |T / 6| > 1 for T = 11, 10, 12: a domain error is NaN, not an exception.
    assert np.isnan(got[[2, 3, 5]]).all() and np.isfinite(got[[0, 1, 4]]).all()


def test_log_of_zero_and_of_a_negative():
    m = data_mesh()
    out = mp.data_calc(m, "log(T)", location="point", output="o")
    got = out.point_data["o"]
    assert got[0] == -np.inf
    out = mp.data_calc(m, "log(-T)", location="point", output="o")
    assert np.isnan(out.point_data["o"][1])


def test_binary_math_functions_broadcast_like_min_and_max():
    m = data_mesh()
    out = mp.data_calc(m, "pow(T, 2)", location="point", output="o")
    np.testing.assert_allclose(out.point_data["o"], m.point_data["T"] ** 2)
    out = mp.data_calc(m, "pow(2, T / 4)", location="point", output="o")
    np.testing.assert_allclose(out.point_data["o"], 2.0 ** (m.point_data["T"] / 4))
    out = mp.data_calc(m, "atan2(v, 2)", location="point", output="o")
    assert out.point_data["o"].shape == (6, 3)
    np.testing.assert_allclose(out.point_data["o"], np.arctan2(m.point_data["v"], 2.0))
    # pow of a negative base to a fractional exponent is a domain error: NaN.
    out = mp.data_calc(m, "pow(0 - T, 0.5)", location="point", output="o")
    assert np.isnan(out.point_data["o"][1])


def test_math_functions_keep_a_vector_a_vector():
    m = data_mesh()
    out = mp.data_calc(m, "exp(v)", location="point", output="o")
    assert out.point_data["o"].shape == (6, 3)
    np.testing.assert_allclose(out.point_data["o"], np.exp(m.point_data["v"]))


def test_math_function_arity():
    m = data_mesh()
    with pytest.raises(ValueError, match="takes exactly 2 arguments"):
        mp.data_calc(m, "pow(T)", location="point", output="o")
    with pytest.raises(ValueError, match="takes exactly 1 argument"):
        mp.data_calc(m, "sin(T, T)", location="point", output="o")


def test_math_functions_match_between_the_core_and_the_twin():
    pytest.importorskip("meshioplusplus._core")
    m = data_mesh()
    exprs = [f"{n}((T - 5.5) / 4)" for n in sorted(_MATH)] + [
        "asin(T / 12)",
        "acos(T / 6)",
        "atan2(T - 5, v)",
        "pow(T, 0.5)",
        "pow(T - 5, 2)",
        "log10(T + 1) * exp(-T / 8)",
    ]
    for expr in exprs:
        cpp = mp.data_calc(m, expr, location="point", output="o")
        py = _calc_py(m, expr, "point_data", "o", False)
        np.testing.assert_allclose(
            cpp.point_data["o"],
            py.point_data["o"],
            rtol=1e-12,
            equal_nan=True,
            err_msg=expr,
        )


def tensor_mesh():
    m = data_mesh()
    t9 = np.array(
        [
            [1, 2, 3, 0, 1, 4, 5, 6, 0],
            [2, 0, 0, 0, 3, 0, 0, 0, 5],
            [1, 0, 0, 0, 1, 0, 0, 0, 1],
            [0, 1, 0, -1, 0, 0, 0, 0, 2],
            [4, 1, 2, 1, 5, 3, 2, 3, 6],
            [1, 1, 1, 1, 1, 1, 1, 1, 1],
        ],
        dtype=float,
    )
    m.point_data["t9"] = t9
    # The same symmetric tensors as 6 components: xx yy zz xy yz zx.
    m.point_data["t6"] = (
        np.array([[2, 3, 5, 0.5, 0.25, 0.75]] * 6, dtype=float)
        * np.arange(1, 7)[:, None]
    )
    return m, t9


def test_dot_and_cross_match_numpy():
    m, _ = tensor_mesh()
    v = m.point_data["v"]
    out = mp.data_calc(m, "dot(v, v)", location="point", output="o")
    np.testing.assert_allclose(out.point_data["o"], (v * v).sum(axis=1))
    out = mp.data_calc(m, "cross(v, v * 0 + 1)", location="point", output="o")
    np.testing.assert_allclose(out.point_data["o"], np.cross(v, np.ones((6, 3))))
    out = mp.data_calc(m, "dot(cross(v, 1 + v), v)", location="point", output="o")
    np.testing.assert_allclose(out.point_data["o"], 0.0, atol=1e-12)  # a x b . a = 0


def test_det_trace_and_transpose_of_a_nine_component_tensor():
    m, t9 = tensor_mesh()
    a = t9.reshape(-1, 3, 3)
    out = mp.data_calc(m, "det(t9)", location="point", output="o")
    np.testing.assert_allclose(out.point_data["o"], np.linalg.det(a), atol=1e-10)
    assert out.point_data["o"][0] == pytest.approx(1.0)  # the integer-inverse matrix
    out = mp.data_calc(m, "trace(t9)", location="point", output="o")
    np.testing.assert_allclose(out.point_data["o"], np.trace(a, axis1=1, axis2=2))
    out = mp.data_calc(m, "transpose(t9)", location="point", output="o")
    assert out.point_data["o"].shape == (6, 9)
    np.testing.assert_allclose(out.point_data["o"], a.transpose(0, 2, 1).reshape(-1, 9))
    out = mp.data_calc(m, "transpose(transpose(t9))", location="point", output="o")
    np.testing.assert_allclose(out.point_data["o"], t9)


def test_det_and_trace_of_a_symmetric_six_component_tensor():
    m, _ = tensor_mesh()
    s = m.point_data["t6"]
    full = np.zeros((6, 3, 3))
    full[:, 0, 0], full[:, 1, 1], full[:, 2, 2] = s[:, 0], s[:, 1], s[:, 2]
    full[:, 0, 1] = full[:, 1, 0] = s[:, 3]
    full[:, 1, 2] = full[:, 2, 1] = s[:, 4]
    full[:, 0, 2] = full[:, 2, 0] = s[:, 5]
    out = mp.data_calc(m, "det(t6)", location="point", output="o")
    np.testing.assert_allclose(out.point_data["o"], np.linalg.det(full), rtol=1e-10)
    out = mp.data_calc(m, "trace(t6)", location="point", output="o")
    np.testing.assert_allclose(out.point_data["o"], np.trace(full, axis1=1, axis2=2))


def test_mag_is_norm():
    m, _ = tensor_mesh()
    a = mp.data_calc(m, "mag(v)", location="point", output="o")
    b = mp.data_calc(m, "norm(v)", location="point", output="o")
    np.testing.assert_array_equal(a.point_data["o"], b.point_data["o"])


def test_tuple_function_width_errors():
    m, _ = tensor_mesh()
    for expr, needle in [
        ("dot(v, t9)", "same width"),
        ("cross(v, t9)", "3-component"),
        ("det(v)", "6-component"),
        ("trace(T)", "6-component"),
        ("transpose(t6)", "9-component"),
    ]:
        with pytest.raises(ValueError, match=needle):
            mp.data_calc(m, expr, location="point", output="o")


def test_tuple_functions_match_between_the_core_and_the_twin():
    pytest.importorskip("meshioplusplus._core")
    m, _ = tensor_mesh()
    for expr in [
        "dot(v, v * 2)",
        "cross(v, 1 + v)",
        "det(t9)",
        "trace(t9) * 2",
        "transpose(t9)",
        "det(t6)",
        "trace(t6)",
        "mag(cross(v, 2 + v))",
        "dot(transpose(t9), t9)",
    ]:
        cpp = mp.data_calc(m, expr, location="point", output="o")
        py = _calc_py(m, expr, "point_data", "o", False)
        np.testing.assert_allclose(
            cpp.point_data["o"],
            py.point_data["o"],
            rtol=1e-12,
            equal_nan=True,
            err_msg=expr,
        )


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
