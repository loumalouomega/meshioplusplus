"""Tests for roadmap §8.1 complex-field convention and helpers."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._data_calc import _calc_py
from meshioplusplus._data_info import _info_py
from meshioplusplus._interpolate import _interpolate_py
from meshioplusplus.ensight import _ensight as py_ensight
from meshioplusplus.unv import _unv as py_unv


def _make_mesh():
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        dtype=np.float64,
    )
    cells = [("triangle", np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64))]
    m = mp.Mesh(points, cells)
    mp.from_complex(m, "p", np.array([3 + 4j, 6 + 8j, 9 + 12j, 12 + 16j]))
    mp.from_complex(
        m,
        "U",
        np.array(
            [[1 + 0j, 0 + 2j, 0 + 0j]] * 4,
            dtype=np.complex128,
        ),
    )
    mp.from_complex(
        m,
        "flux",
        [np.array([1.5 + 0.5j, -2.5 + 3.5j], dtype=np.complex128)],
        location="cell",
    )
    return m


def test_complex_helpers_roundtrip():
    m = _make_mesh()
    pairs = mp.complex_pairs(m)
    assert pairs["point_data"] == ["U", "p"]
    assert pairs["cell_data"] == ["flux"]

    p_c = mp.to_complex(m, "p", location="point")
    assert p_c.dtype == np.complex128
    np.testing.assert_allclose(p_c, [3 + 4j, 6 + 8j, 9 + 12j, 12 + 16j])

    flux_c = mp.to_complex(m, "flux", location="cell")
    assert len(flux_c) == 1
    np.testing.assert_allclose(flux_c[0], [1.5 + 0.5j, -2.5 + 3.5j])

    all_c = mp.to_complex(m)
    assert set(all_c["point_data"]) == {"U", "p"}
    assert set(all_c["cell_data"]) == {"flux"}

    # float32 pair preserves complex64
    m32 = mp.Mesh(m.points, m.cells)
    mp.from_complex(
        m32, "z", np.array([1 + 2j, 3 + 4j, 5 + 6j, 7 + 8j], dtype=np.complex64)
    )
    assert mp.to_complex(m32, "z").dtype == np.complex64


@pytest.mark.parametrize("engine", ["core", "py"])
def test_unv_complex_roundtrip(tmp_path, engine):
    m = _make_mesh()
    out = tmp_path / f"complex_{engine}.unv"
    if engine == "core":
        mp.write(out, m)
        back = mp.read(out)
    else:
        py_unv.write(str(out), m)
        back = py_unv.read(str(out))
    np.testing.assert_allclose(mp.to_complex(back, "p"), mp.to_complex(m, "p"))
    np.testing.assert_allclose(mp.to_complex(back, "U"), mp.to_complex(m, "U"))
    np.testing.assert_allclose(
        mp.to_complex(back, "flux", location="cell")[0],
        mp.to_complex(m, "flux", location="cell")[0],
    )


@pytest.mark.parametrize("engine", ["core", "py"])
@pytest.mark.parametrize("binary", [False, True])
def test_ensight_complex_roundtrip(tmp_path, engine, binary):
    m = _make_mesh()
    out = tmp_path / f"complex_{engine}_{binary}.case"
    if engine == "core":
        mp.write(out, m, binary=binary)
        back = mp.read(out)
    else:
        py_ensight.write(str(out), m, binary=binary)
        back = py_ensight.read(str(out))
    case_text = out.read_text()
    assert "complex scalar per node: p p_real.scl p_imag.scl" in case_text
    assert "complex vector per node: U U_real.vec U_imag.vec" in case_text
    assert "complex scalar per element: flux flux_real.escl flux_imag.escl" in case_text
    np.testing.assert_allclose(
        mp.to_complex(back, "p"), mp.to_complex(m, "p"), rtol=1e-5
    )
    np.testing.assert_allclose(
        mp.to_complex(back, "U"), mp.to_complex(m, "U"), rtol=1e-5
    )
    np.testing.assert_allclose(
        mp.to_complex(back, "flux", location="cell")[0],
        mp.to_complex(m, "flux", location="cell")[0],
        rtol=1e-5,
    )


def test_cross_format_op2_unv_complex_roundtrip(tmp_path):
    """Complex OP2 <-> UNV round-trip yields identical names, shapes, and values."""
    m = mp.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2]], dtype=np.int64))],
    )
    mp.from_complex(
        m,
        "DISPLACEMENT",
        np.array(
            [
                [1 + 2j, 3 + 4j, 5 + 6j],
                [2 - 1j, 0 + 3j, -1 + 2j],
                [4 + 0j, 1 - 1j, 2 + 2j],
            ],
            dtype=np.complex128,
        ),
    )
    op2_path = tmp_path / "c.op2"
    unv_path = tmp_path / "c.unv"
    mp.write(op2_path, m)
    m_op2 = mp.read(op2_path)
    assert "DISPLACEMENT" in mp.complex_pairs(m_op2)["point_data"]
    mp.write(unv_path, m_op2)
    m_unv = mp.read(unv_path)
    np.testing.assert_allclose(
        mp.to_complex(m_unv, "DISPLACEMENT"),
        mp.to_complex(m, "DISPLACEMENT"),
        rtol=1e-5,
    )


def test_data_calc_abs_phase_both_engines():
    m = _make_mesh()
    for fn in (
        lambda expr, out: mp.data_calc(m, expr, output=out),
        lambda expr, out: _calc_py(m, expr, "point_data", out, True),
    ):
        res_mag = fn("abs(p)", "p_mag")
        np.testing.assert_allclose(res_mag.point_data["p_mag"], [5.0, 10.0, 15.0, 20.0])
        res_ph = fn("phase(p)", "p_phase")
        np.testing.assert_allclose(
            res_ph.point_data["p_phase"], np.full(4, np.arctan2(4.0, 3.0))
        )
        res_hy = fn("hypot(p_real, p_imag)", "p_hy")
        np.testing.assert_allclose(res_hy.point_data["p_hy"], [5.0, 10.0, 15.0, 20.0])


def test_interpolate_complex_base_name_both_engines():
    src = _make_mesh()
    dst = mp.Mesh(src.points.copy(), src.cells)
    out_core = mp.interpolate(src, dst, arrays=["p"])
    out_py = _interpolate_py(
        src,
        dst,
        method="barycentric",
        arrays=["p"],
        extrapolate=False,
        default_value=0.0,
        on_conflict="overwrite",
    )
    for out in (out_core, out_py):
        np.testing.assert_allclose(mp.to_complex(out, "p"), mp.to_complex(src, "p"))


def test_data_info_reports_complex_pairs():
    m = _make_mesh()
    for info in (mp.data_info(m), _info_py(m)):
        by_key = {(a["location"], a["name"]): a for a in info}
        assert by_key[("point_data", "p_real")]["complex_pair"] == "p"
        assert by_key[("point_data", "p_real")]["complex_part"] == "real"
        assert by_key[("point_data", "p_imag")]["complex_pair"] == "p"
        assert by_key[("point_data", "p_imag")]["complex_part"] == "imag"
        assert by_key[("cell_data", "flux_real")]["complex_pair"] == "flux"
