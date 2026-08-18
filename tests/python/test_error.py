"""Tests for the ZZ recovery-based error estimator ``meshioplusplus.estimate_error``.

The central oracle: ``gradient``'s Green-Gauss is exact for a linear field, so
the raw and recovered gradients agree everywhere and ``error:zz`` must be (up to
rounding) zero mesh-wide. A quadratic field breaks that exactness on purpose,
giving a nonzero, non-uniform indicator to mark against.
"""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._error import _estimate_error_py


def _hex_grid(n=3):
    """An ``n x n x n`` unit-cube hexahedron grid over ``[0, n]^3``."""
    m = n + 1
    pts = np.array(
        [[i, j, k] for k in range(m) for j in range(m) for i in range(m)], dtype=float
    )

    def idx(i, j, k):
        return (k * m + j) * m + i

    conn = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                conn.append(
                    [
                        idx(i, j, k),
                        idx(i + 1, j, k),
                        idx(i + 1, j + 1, k),
                        idx(i, j + 1, k),
                        idx(i, j, k + 1),
                        idx(i + 1, j, k + 1),
                        idx(i + 1, j + 1, k + 1),
                        idx(i, j + 1, k + 1),
                    ]
                )
    return mp.Mesh(pts, [("hexahedron", np.array(conn, dtype=np.int64))])


def _linear_mesh():
    m = _hex_grid()
    m.point_data["f"] = (
        2.0 * m.points[:, 0] - 3.0 * m.points[:, 1] + 0.5 * m.points[:, 2] + 7.0
    )
    return m


def _quadratic_mesh():
    m = _hex_grid()
    m.point_data["f"] = (m.points**2).sum(axis=1)
    return m


def _flat_zz(mesh, name="error:zz"):
    return np.concatenate(
        [np.asarray(b, dtype=float).reshape(-1) for b in mesh.cell_data[name]]
    )


def _flat_marked(mesh, name="error:marked"):
    return np.concatenate(
        [np.asarray(b, dtype=np.int64).reshape(-1) for b in mesh.cell_data[name]]
    )


# --------------------------------------------------------------------------- #
# validation                                                                   #
# --------------------------------------------------------------------------- #


def test_empty_array_name_raises():
    m = _linear_mesh()
    with pytest.raises(ValueError):
        mp.estimate_error(m, "")


def test_unknown_array_raises():
    m = _linear_mesh()
    with pytest.raises(ValueError):
        mp.estimate_error(m, "nope")


def test_cell_data_array_names_the_fix():
    m = _linear_mesh()
    m.cell_data["f_cell"] = [np.zeros(len(b)) for b in m.cells]
    with pytest.raises(ValueError, match="cell_data_to_point_data"):
        mp.estimate_error(m, "f_cell")


def test_unknown_method_raises():
    m = _linear_mesh()
    with pytest.raises(ValueError):
        mp.estimate_error(m, "f", method="kelly")


def test_unknown_marking_raises():
    m = _linear_mesh()
    with pytest.raises(ValueError):
        mp.estimate_error(m, "f", marking="median")


@pytest.mark.parametrize("marking", ["fraction", "dorfler"])
@pytest.mark.parametrize("value", [0.0, -0.1, 1.5])
def test_marking_value_out_of_range_raises(marking, value):
    m = _linear_mesh()
    with pytest.raises(ValueError):
        mp.estimate_error(m, "f", marking=marking, marking_value=value)


def test_output_collision_without_overwrite_raises():
    m = _linear_mesh()
    out = mp.estimate_error(m, "f")
    with pytest.raises(ValueError):
        mp.estimate_error(out, "f")
    # overwrite=True is fine.
    mp.estimate_error(out, "f", overwrite=True)


# --------------------------------------------------------------------------- #
# the exactness oracle                                                         #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("engine", ["auto", "python"])
def test_linear_field_gives_zero_indicator_everywhere(engine):
    m = _linear_mesh()
    if engine == "python":
        out, ge, skipped, marked = _estimate_error_py(
            m, "f", "none", 0.0, None, None, False
        )
    else:
        out, report = mp.estimate_error(m, "f", return_report=True)
        ge, skipped, marked = (
            report["global_error"],
            report["num_skipped"],
            report["num_marked"],
        )
    assert skipped == 0
    assert marked == 0
    assert ge == pytest.approx(0.0, abs=1e-9)
    assert _flat_zz(out) == pytest.approx(np.zeros(27), abs=1e-9)


def test_linear_field_marks_nothing_under_a_real_absolute_threshold():
    # Dorfler/Fraction always mark *something* on a near-zero-but-not-exactly-
    # zero indicator (their target is proportional to the total, however
    # tiny) -- not a meaningful invariant to pin. A real threshold well above
    # the rounding-noise floor is.
    m = _linear_mesh()
    out, report = mp.estimate_error(
        m, "f", marking="absolute", marking_value=1e-6, return_report=True
    )
    assert report["num_marked"] == 0
    assert _flat_marked(out).sum() == 0


# --------------------------------------------------------------------------- #
# a genuine, non-uniform indicator                                             #
# --------------------------------------------------------------------------- #


def test_quadratic_field_gives_nonzero_nonuniform_indicator():
    m = _quadratic_mesh()
    out, report = mp.estimate_error(m, "f", return_report=True)
    assert report["num_skipped"] == 0
    assert report["global_error"] > 0.0
    zz = _flat_zz(out)
    assert zz.shape == (27,)
    assert (zz >= 0.0).all()
    assert zz.max() - zz.min() > 1e-9
    assert report["global_error"] == pytest.approx(np.sqrt(np.sum(zz**2)), rel=1e-6)


def test_marking_none_attaches_no_marked_array():
    m = _quadratic_mesh()
    out = mp.estimate_error(m, "f")
    assert "error:marked" not in out.cell_data


def test_absolute_marks_exactly_cells_above_threshold():
    m = _quadratic_mesh()
    plain = mp.estimate_error(m, "f")
    zz = _flat_zz(plain)
    threshold = float(np.median(zz))

    out, report = mp.estimate_error(
        m, "f", marking="absolute", marking_value=threshold, return_report=True
    )
    zz2 = _flat_zz(out)
    marked = _flat_marked(out)
    expected = zz2 > threshold
    assert (marked == expected.astype(np.int64)).all()
    assert report["num_marked"] == int(expected.sum())


def test_fraction_marks_the_largest_indicators_and_the_right_count():
    m = _quadratic_mesh()
    fraction = 0.25
    out, report = mp.estimate_error(
        m, "f", marking="fraction", marking_value=fraction, return_report=True
    )
    zz = _flat_zz(out)
    marked = _flat_marked(out)
    expected_count = int(fraction * len(zz) + 0.5)
    assert report["num_marked"] == expected_count
    if 0 < expected_count < len(zz):
        assert zz[marked == 1].min() >= zz[marked == 0].max()


def test_dorfler_marks_the_minimal_bulk_prefix():
    m = _quadratic_mesh()
    theta = 0.6
    out, report = mp.estimate_error(
        m, "f", marking="dorfler", marking_value=theta, return_report=True
    )
    zz = _flat_zz(out)
    marked = _flat_marked(out)
    total_sq = float(np.sum(zz**2))
    marked_sq = float(np.sum(zz[marked == 1] ** 2))
    assert report["num_marked"] == int(marked.sum())
    assert marked_sq >= theta * total_sq - 1e-9
    if marked.sum() > 0:
        min_marked = float(zz[marked == 1].min())
        assert marked_sq - min_marked**2 < theta * total_sq - 1e-9


def test_custom_output_names():
    m = _quadratic_mesh()
    out = mp.estimate_error(
        m,
        "f",
        marking="absolute",
        marking_value=0.5,
        output="my_ind",
        marked_name="my_marks",
    )
    assert "my_ind" in out.cell_data
    assert "my_marks" in out.cell_data
    assert "error:zz" not in out.cell_data
    assert "error:marked" not in out.cell_data


def test_geometry_and_unrelated_data_pass_through_unchanged():
    m = _quadratic_mesh()
    m.cell_data["tag"] = [np.full(len(b), 42, dtype=np.int64) for b in m.cells]
    out = mp.estimate_error(m, "f")
    assert len(out.points) == len(m.points)
    assert len(out.cells) == len(m.cells)
    assert list(out.cell_data["tag"][0]) == [42] * len(m.cells[0])
    assert "f" in out.point_data


# --------------------------------------------------------------------------- #
# C++ / numpy parity                                                           #
# --------------------------------------------------------------------------- #


def test_cpp_matches_python():
    m = _quadratic_mesh()
    cpp, cpp_report = mp.estimate_error(
        m, "f", marking="dorfler", marking_value=0.6, return_report=True
    )
    out, ge, skipped, marked = _estimate_error_py(
        m, "f", "dorfler", 0.6, None, None, False
    )

    # Same, already-accepted tolerance data_average's own Measure weighting
    # has (test_data_location.py::test_cpp_matches_python): the recovery step
    # composes that exact weighting.
    assert _flat_zz(out) == pytest.approx(_flat_zz(cpp), abs=1e-8)
    assert _flat_marked(out).tolist() == _flat_marked(cpp).tolist()
    assert ge == pytest.approx(cpp_report["global_error"], rel=1e-6)
    assert skipped == cpp_report["num_skipped"]
    assert marked == cpp_report["num_marked"]


def test_python_twin_raises_on_polyhedron_block():
    pts = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ]
    faces = [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ]
    m = mp.Mesh(pts, [("polyhedron6", [faces])])
    m.point_data["f"] = np.asarray(pts)[:, 0]
    with pytest.raises(NotImplementedError):
        _estimate_error_py(m, "f", "none", 0.0, None, None, False)


def test_python_twin_raises_on_quadratic_3d_type():
    # tetra10: the 4 tetra corners plus 6 edge midpoints.
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [0.5, 0, 0],
            [0.5, 0.5, 0],
            [0, 0.5, 0],
            [0, 0, 0.5],
            [0.5, 0, 0.5],
            [0, 0.5, 0.5],
        ],
        dtype=float,
    )
    conn = np.array([list(range(10))], dtype=np.int64)
    m = mp.Mesh(pts, [("tetra10", conn)])
    m.point_data["f"] = m.points[:, 0]
    with pytest.raises(NotImplementedError):
        _estimate_error_py(m, "f", "none", 0.0, None, None, False)
