"""Tests for the tensor_invariants data operation."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._tensor_invariants import _tensor_invariants_py

from .helpers_data import assert_same_geometry


def tensor_mesh():
    """A single-triangle mesh with a 6-component and a 9-component array."""
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    m = mp.Mesh(points, [("triangle", np.array([[0, 1, 2]]))])
    m.point_data["s"] = np.array(
        [
            [1.0, 2.0, 3.0, 0.5, 0.6, 0.7],
            [2.0, 2.0, 2.0, 0.0, 0.0, 0.0],
            [np.nan, 1.0, 1.0, 0.0, 0.0, 0.0],
        ]
    )
    m.point_data["g9"] = np.array(
        [
            [1.0, 0.0, 0.0, 0.0, 2.0, 0.5, 0.0, 0.5, 3.0],
            [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
            [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        ]
    )
    return m


def test_mises_matches_closed_form():
    m = tensor_mesh()
    out = mp.tensor_invariants(m, "point", ["s"], outputs=["mises"])
    xx, yy, zz, xy, yz, zx = 1.0, 2.0, 3.0, 0.5, 0.6, 0.7
    expect = np.sqrt(
        0.5
        * (
            (xx - yy) ** 2
            + (yy - zz) ** 2
            + (zz - xx) ** 2
            + 6.0 * (xy**2 + yz**2 + zx**2)
        )
    )
    assert out.point_data["s_mises"][0] == pytest.approx(expect)
    assert out.point_data["s_mises"][1] == pytest.approx(0.0, abs=1e-12)  # isotropic
    assert np.isnan(out.point_data["s_mises"][2])
    assert_same_geometry(m, out)


def test_principal_ascending_and_sums_to_trace():
    m = tensor_mesh()
    out = mp.tensor_invariants(m, "point", ["s"], outputs=["principal"])
    p = out.point_data["s_principal"]
    assert (np.diff(p[0]) >= 0).all()
    assert p[0].sum() == pytest.approx(1.0 + 2.0 + 3.0)
    assert p[1] == pytest.approx([2.0, 2.0, 2.0])  # isotropic
    assert np.isnan(p[2]).all()


def test_hydrostatic_is_mean_diagonal():
    m = tensor_mesh()
    out = mp.tensor_invariants(m, "point", ["s"], outputs=["hydrostatic"])
    assert out.point_data["s_hydrostatic"][0] == pytest.approx((1.0 + 2.0 + 3.0) / 3.0)
    assert out.point_data["s_hydrostatic"][1] == pytest.approx(2.0)


def test_deviatoric_subtracts_hydrostatic_from_diagonal_only():
    m = tensor_mesh()
    out = mp.tensor_invariants(m, "point", ["s"], outputs=["deviatoric"])
    d = out.point_data["s_deviatoric"]
    h = (1.0 + 2.0 + 3.0) / 3.0
    assert d[0] == pytest.approx([1.0 - h, 2.0 - h, 3.0 - h, 0.5, 0.6, 0.7])
    assert d[1] == pytest.approx([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # isotropic


def test_nine_component_uses_symmetric_part():
    m = tensor_mesh()
    out = mp.tensor_invariants(m, "point", ["g9"], outputs=["mises", "deviatoric"])
    expect = np.sqrt(
        0.5 * ((1 - 2) ** 2 + (2 - 3) ** 2 + (3 - 1) ** 2 + 6.0 * (0.5**2))
    )
    assert out.point_data["g9_mises"][0] == pytest.approx(expect)
    dev = out.point_data["g9_deviatoric"]
    assert dev.shape[1] == 9
    assert dev[0, 1] == pytest.approx(0.0)  # xy unchanged
    assert dev[0, 5] == pytest.approx(0.5)  # yz unchanged
    assert dev[0, 7] == pytest.approx(0.5)  # zy unchanged (not symmetrized)


def test_field_location_raises():
    m = tensor_mesh()
    with pytest.raises(ValueError):
        mp.tensor_invariants(m, "field")


def test_unknown_name_raises():
    m = tensor_mesh()
    with pytest.raises(ValueError):
        mp.tensor_invariants(m, "point", ["nope"])


def test_wrong_component_count_raises_when_named_explicitly():
    m = tensor_mesh()
    m.point_data["T"] = np.array([0.0, 1.0, 2.0])
    with pytest.raises(ValueError):
        mp.tensor_invariants(m, "point", ["T"])


def test_auto_select_skips_non_tensor_shaped_arrays():
    m = tensor_mesh()
    m.point_data["T"] = np.array([0.0, 1.0, 2.0])
    out = mp.tensor_invariants(m, "point", outputs=["mises"])
    assert "s_mises" in out.point_data
    assert "g9_mises" in out.point_data
    assert "T_mises" not in out.point_data


def test_overwrite_false_refuses_collision():
    m = tensor_mesh()
    m.point_data["s_mises"] = np.zeros(3)
    with pytest.raises(ValueError):
        mp.tensor_invariants(m, "point", ["s"], outputs=["mises"], overwrite=False)


def test_cell_data_one_output_per_block():
    points = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [2, 0, 0]])
    m = mp.Mesh(
        points,
        [("triangle", np.array([[0, 1, 2]])), ("triangle", np.array([[0, 2, 3]]))],
    )
    m.cell_data["g"] = [
        np.array([[1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0]]),
        np.array([[2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0]]),
    ]
    out = mp.tensor_invariants(m, "cell", ["g"], outputs=["hydrostatic"])
    assert len(out.cell_data["g_hydrostatic"]) == 2
    assert out.cell_data["g_hydrostatic"][0][0] == pytest.approx(2.0)
    assert out.cell_data["g_hydrostatic"][1][0] == pytest.approx(2.0)
    assert_same_geometry(m, out)


def test_prefix_and_suffix():
    m = tensor_mesh()
    out = mp.tensor_invariants(
        m, "point", ["s"], outputs=["mises"], prefix="p_", suffix="_q"
    )
    assert "p_s_mises_q" in out.point_data


def test_cpp_matches_python():
    pytest.importorskip("meshioplusplus._core")
    m = tensor_mesh()
    for loc, names in (("point_data", ["s"]), ("point_data", ["g9"])):
        cpp = mp.tensor_invariants(m, loc.replace("_data", ""), names)
        py = _tensor_invariants_py(m, loc, names, [], "", "", True)
        for out_name in ("mises", "principal", "hydrostatic", "deviatoric"):
            key = f"{names[0]}_{out_name}"
            assert np.allclose(
                cpp.point_data[key], py.point_data[key], equal_nan=True
            ), (loc, names, out_name)
