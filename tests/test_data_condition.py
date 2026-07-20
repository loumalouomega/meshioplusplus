"""Tests for data value conditioning (clamp / normalize / standardize)."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._data_condition import _condition_py

from .helpers_data import assert_same_geometry, data_mesh


def test_clamp_bounds_values_correctly():
    m = data_mesh()
    out = mp.data_condition(m, "point", ["T"], mode="clamp", lo=1.0, hi=10.0)
    # T = {0, 1, 11, 10, 2, 12}
    assert out.point_data["T"] == pytest.approx([1.0, 1.0, 10.0, 10.0, 2.0, 10.0])
    assert_same_geometry(m, out)


def test_normalize_maps_min_to_zero_and_max_to_one():
    m = data_mesh()
    out = mp.data_condition(m, "point", ["T"], mode="normalize")
    got = out.point_data["T"]
    assert got.min() == pytest.approx(0.0)
    assert got.max() == pytest.approx(1.0)
    assert got == pytest.approx(m.point_data["T"] / 12.0)


def test_normalize_to_a_custom_range():
    m = data_mesh()
    out = mp.data_condition(m, "point", ["T"], mode="normalize", lo=-1.0, hi=1.0)
    got = out.point_data["T"]
    assert got.min() == pytest.approx(-1.0)
    assert got.max() == pytest.approx(1.0)


def test_standardize_gives_zero_mean_unit_std():
    m = data_mesh()
    out = mp.data_condition(m, "point", ["T"], mode="standardize")
    got = out.point_data["T"]
    assert got.mean() == pytest.approx(0.0, abs=1e-12)
    assert got.std() == pytest.approx(1.0)


def test_constant_array_normalizes_to_the_target_lower_bound():
    m = data_mesh()
    m.point_data["k"] = np.full(6, 5.0)
    out = mp.data_condition(m, "point", ["k"], mode="normalize")
    assert out.point_data["k"] == pytest.approx(np.zeros(6))


def test_component_scope_conditions_each_component_independently():
    m = data_mesh()
    out = mp.data_condition(m, "point", ["v"], mode="normalize")
    got = out.point_data["v"]
    # Component 0 of v is {1,0,0,1,2,0}: min 0, max 2.
    assert got[:, 0] == pytest.approx([0.5, 0.0, 0.0, 0.5, 1.0, 0.0])


def test_magnitude_scope_preserves_direction():
    m = data_mesh()
    out = mp.data_condition(m, "point", ["v"], mode="normalize", scope="magnitude")
    got = out.point_data["v"]
    # Row 3 was (1,1,0): it must still point along (1,1,0).
    assert got[3, 0] == pytest.approx(got[3, 1])
    assert got[3, 2] == pytest.approx(0.0)


def test_clamp_preserves_integer_dtype():
    # MESHIO-backend-pinned: NATIVE/KRATOS canonicalize widths, but the Python
    # path always uses the MESHIO backend.
    m = data_mesh()
    out = mp.data_condition(m, "cell", ["tag"], mode="clamp", lo=0.0, hi=25.0)
    assert np.issubdtype(np.asarray(out.cell_data["tag"][0]).dtype, np.integer)
    assert out.cell_data["tag"][1][0] == 25


def test_normalize_always_produces_float():
    m = data_mesh()
    out = mp.data_condition(m, "cell", ["tag"], mode="normalize")
    assert np.issubdtype(np.asarray(out.cell_data["tag"][0]).dtype, np.floating)


def test_cell_data_statistics_are_joint_across_blocks():
    m = data_mesh()
    out = mp.data_condition(m, "cell", ["mat"], mode="normalize")
    # mat spans both blocks: {1,2} and {3}. Joint min 1, max 3.
    assert out.cell_data["mat"][0] == pytest.approx([0.0, 0.5])
    assert out.cell_data["mat"][1] == pytest.approx([1.0])
    assert len(out.cell_data["mat"]) == len(m.cells)


def test_nan_is_excluded_from_reductions_and_passed_through():
    m = data_mesh()
    m.point_data["n"] = np.array([0.0, np.nan, 4.0, 2.0, 1.0, 3.0])
    out = mp.data_condition(m, "point", ["n"], mode="normalize")
    got = out.point_data["n"]
    # min 0, max 4 with the NaN ignored.
    assert got[0] == pytest.approx(0.0)
    assert got[2] == pytest.approx(1.0)
    assert np.isnan(got[1])


def test_nan_replace_policy():
    m = data_mesh()
    m.point_data["n"] = np.array([0.0, np.nan, 4.0, 2.0, 1.0, 3.0])
    out = mp.data_condition(
        m, "point", ["n"], mode="normalize", nan_policy="replace", nan_replacement=-1.0
    )
    assert out.point_data["n"][1] == pytest.approx(-1.0)


def test_nan_fail_policy_raises():
    m = data_mesh()
    m.point_data["n"] = np.array([0.0, np.nan, 4.0, 2.0, 1.0, 3.0])
    with pytest.raises(ValueError):
        mp.data_condition(m, "point", ["n"], mode="normalize", nan_policy="fail")


def test_inf_is_excluded_from_reductions():
    m = data_mesh()
    m.point_data["n"] = np.array([0.0, np.inf, 4.0, 2.0, 1.0, 3.0])
    out = mp.data_condition(m, "point", ["n"], mode="normalize")
    # max is 4, not inf.
    assert out.point_data["n"][2] == pytest.approx(1.0)


def test_suffix_leaves_the_original_alone():
    m = data_mesh()
    out = mp.data_condition(m, "point", ["T"], mode="normalize", suffix="_n")
    assert out.point_data["T_n"].max() == pytest.approx(1.0)
    assert out.point_data["T"] == pytest.approx(m.point_data["T"])


def test_inverted_bounds_raise():
    m = data_mesh()
    with pytest.raises(ValueError):
        mp.data_condition(m, "point", ["T"], mode="clamp", lo=10.0, hi=1.0)


def test_unknown_name_raises():
    m = data_mesh()
    with pytest.raises(ValueError):
        mp.data_condition(m, "point", ["nope"], mode="clamp")


def test_cpp_matches_python():
    pytest.importorskip("meshioplusplus._core")
    m = data_mesh()
    for mode in ("clamp", "normalize", "standardize"):
        for scope in ("component", "magnitude"):
            cpp = mp.data_condition(m, "point", ["v"], mode=mode, scope=scope)
            py = _condition_py(
                m, "point_data", ["v"], mode, scope, 0.0, 1.0, "ignore", 0.0, "", True
            )
            assert np.allclose(
                cpp.point_data["v"], py.point_data["v"], equal_nan=True
            ), (mode, scope)


def test_roundtrip_write_read(tmp_path):
    m = data_mesh()
    out = mp.data_condition(m, "point", ["T"], mode="normalize")
    path = tmp_path / "out.vtu"
    mp.write(path, out)
    back = mp.read(path)
    assert np.allclose(back.point_data["T"], out.point_data["T"])
    assert back.point_data["T"].min() == pytest.approx(0.0)
    assert back.point_data["T"].max() == pytest.approx(1.0)
