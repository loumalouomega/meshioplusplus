"""Tests for the Hessian composition ``meshioplusplus.hessian``.

The one MESH-SHAPE-INDEPENDENT guarantee is "a linear field has an exactly
zero Hessian" (see ``operations/hessian.hpp``'s own doc comment): a linear
field's gradient is EXACTLY constant after the first Green-Gauss + Point-
averaging pass, and Green-Gauss of a spatially constant field is trivially
exact, so the second pass gives exactly zero regardless of mesh regularity.

For a genuinely quadratic field, the tests below measure (not assume) the
honest picture: on a regular axis-aligned hex grid, an INTERIOR cell's
neighbourhood is symmetric enough for the intermediate Uniform-weighted
Point-averaging step to be exact too, so the composed Hessian comes back
exact to machine precision there -- mirrors
``tests/cpp/test_hessian.cpp``'s identical finding.

``_hessian_py`` is not a separate numpy kernel: it composes the already
C++/numpy-parity-tested public :func:`meshioplusplus.gradient` twice, so
``test_cpp_matches_python`` below is exact equality, not a tolerance --
verified empirically before being pinned that way.
"""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._hessian import _hessian_py

# --------------------------------------------------------------------------- #
# fixtures                                                                     #
# --------------------------------------------------------------------------- #


def _hex_grid(n):
    """A regular n x n x n hexahedron grid on [0,n]^3."""
    pts = []
    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                pts.append([float(i), float(j), float(k)])
    pts = np.array(pts)

    def vid(i, j, k):
        return (k * (n + 1) + j) * (n + 1) + i

    cells = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                cells.append(
                    [
                        vid(i, j, k),
                        vid(i + 1, j, k),
                        vid(i + 1, j + 1, k),
                        vid(i, j + 1, k),
                        vid(i, j, k + 1),
                        vid(i + 1, j, k + 1),
                        vid(i + 1, j + 1, k + 1),
                        vid(i, j + 1, k + 1),
                    ]
                )
    cells = np.array(cells, dtype=np.int64)
    return mp.Mesh(pts, [("hexahedron", cells)])


def _tapered_hex():
    """A frustum: a 2x2 base at z=0 tapering to a 1x1 top at z=1 -- an
    irregular single cell, used to show the quadratic-field composition is
    approximate (not exact) off a structured grid."""
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.5, 0.5, 1.0],
            [1.5, 0.5, 1.0],
            [1.5, 1.5, 1.0],
            [0.5, 1.5, 1.0],
        ]
    )
    return mp.Mesh(pts, [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))])


def _linear_scalar(points, a, b, c, d):
    return points @ np.array([a, b, c]) + d


def _quadratic_scalar(points, H):
    return 0.5 * np.einsum("ni,ij,nj->n", points, H, points)


def _cell_idx(n, i, j, k):
    return (k * n + j) * n + i


# --------------------------------------------------------------------------- #
# linear field: exactly zero Hessian, any mesh                                #
# --------------------------------------------------------------------------- #


def test_linear_field_on_a_tapered_hex_is_exactly_zero():
    m = _tapered_hex()
    m.point_data["f"] = _linear_scalar(m.points, 3.0, -2.0, 5.0, 7.0)
    out, report = mp.hessian(m, "f", return_report=True)
    h = np.asarray(out.cell_data["f:hessian"][0]).reshape(-1)
    assert h == pytest.approx(np.zeros(9), abs=1e-9)
    assert report["num_skipped"] == 0


def test_linear_field_on_a_hex_grid_is_exactly_zero():
    m = _hex_grid(3)
    m.point_data["f"] = _linear_scalar(m.points, 0.5, -1.25, 2.0, 4.0)
    out = mp.hessian(m, "f")
    h = np.asarray(out.cell_data["f:hessian"][0])
    assert h == pytest.approx(np.zeros_like(h), abs=1e-9)


# --------------------------------------------------------------------------- #
# quadratic field: exact away from a structured mesh's own boundary           #
# --------------------------------------------------------------------------- #


def test_quadratic_field_on_a_structured_grid_is_exact_away_from_the_boundary():
    H = np.array([[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 4.0]])
    n = 8
    m = _hex_grid(n)
    m.point_data["f"] = _quadratic_scalar(m.points, H)
    out = mp.hessian(m, "f")
    h = np.asarray(out.cell_data["f:hessian"][0]).reshape(-1, 3, 3)
    # At least two rings in from every face.
    for k in range(2, n - 2):
        for j in range(2, n - 2):
            for i in range(2, n - 2):
                c = _cell_idx(n, i, j, k)
                assert h[c] == pytest.approx(H, abs=1e-9)


def test_quadratic_field_on_a_structured_grid_has_bounded_boundary_error():
    H = np.array([[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 4.0]])
    n = 8
    m = _hex_grid(n)
    m.point_data["f"] = _quadratic_scalar(m.points, H)
    out = mp.hessian(m, "f")
    h = np.asarray(out.cell_data["f:hessian"][0]).reshape(-1, 3, 3)
    max_err = np.max(np.abs(h - H[None, :, :]))
    assert max_err > 1e-6  # proves the boundary error is real
    assert max_err < 3.0  # ... but bounded; a regression trip-wire


def test_quadratic_field_on_an_irregular_mesh_is_approximate_but_bounded():
    H = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
    m = _tapered_hex()
    m.point_data["f"] = _quadratic_scalar(m.points, H)
    out = mp.hessian(m, "f")
    h = np.asarray(out.cell_data["f:hessian"][0]).reshape(3, 3)
    max_err = np.max(np.abs(h - H))
    assert max_err > 1e-6
    assert max_err < 3.0


# --------------------------------------------------------------------------- #
# shape / contract                                                             #
# --------------------------------------------------------------------------- #


def test_scalar_input_gives_nine_components():
    m = _hex_grid(2)
    m.point_data["f"] = _linear_scalar(m.points, 1.0, 0.0, 0.0, 0.0)
    out = mp.hessian(m, "f")
    assert np.asarray(out.cell_data["f:hessian"][0]).shape[-1] == 9


def test_rejects_a_multi_component_input_by_name():
    m = _hex_grid(2)
    m.point_data["u"] = np.zeros((m.points.shape[0], 3))
    with pytest.raises(ValueError, match="scalar fields only"):
        mp.hessian(m, "u")


def test_rejects_a_cell_data_input_by_name():
    m = _hex_grid(2)
    m.cell_data["c"] = [np.zeros(len(m.cells[0]))]
    with pytest.raises(ValueError, match="cell_data_to_point_data"):
        mp.hessian(m, "c")


def test_rejects_an_unknown_array_name():
    m = _hex_grid(2)
    m.point_data["f"] = _linear_scalar(m.points, 1.0, 0.0, 0.0, 0.0)
    with pytest.raises(ValueError):
        mp.hessian(m, "nope")


def test_point_location_attaches_point_data():
    m = _hex_grid(3)
    m.point_data["f"] = _linear_scalar(m.points, 1.0, -1.0, 2.0, 0.0)
    out = mp.hessian(m, "f", location="point")
    assert "f:hessian" in out.point_data
    assert "f:hessian" not in out.cell_data
    assert np.asarray(out.point_data["f:hessian"]).shape == (m.points.shape[0], 9)


def test_output_name_and_overwrite():
    m = _hex_grid(2)
    m.point_data["f"] = _linear_scalar(m.points, 1.0, 0.0, 0.0, 0.0)
    out1 = mp.hessian(m, "f", output="curv")
    assert "curv" in out1.cell_data

    with pytest.raises(ValueError):
        mp.hessian(out1, "f", output="curv")
    mp.hessian(out1, "f", output="curv", overwrite=True)  # does not raise


def test_geometry_and_existing_data_pass_through_unchanged():
    m = _hex_grid(2)
    m.point_data["f"] = _linear_scalar(m.points, 1.0, 2.0, 3.0, 0.0)
    m.cell_data["mat"] = [np.arange(8.0)]
    out = mp.hessian(m, "f")
    assert out.points.shape == m.points.shape
    assert np.asarray(out.cell_data["mat"][0]).tolist() == list(range(8))


def test_least_squares_method_runs():
    m = _hex_grid(3)
    m.point_data["f"] = _linear_scalar(m.points, 1.0, -2.0, 0.5, 3.0)
    out = mp.hessian(m, "f", method="least-squares")
    h = np.asarray(out.cell_data["f:hessian"][0])
    assert h == pytest.approx(np.zeros_like(h), abs=1e-6)


def test_hessian_is_exported():
    assert "hessian" in mp.__all__
    assert mp.hessian is mp.hessian


# --------------------------------------------------------------------------- #
# C++ / composition parity                                                     #
# --------------------------------------------------------------------------- #


def test_cpp_matches_python():
    """Pins parity only. Not a tolerance: ``_hessian_py`` composes the same
    public ``gradient()`` C++/numpy-parity-tested function `_core.hessian`
    itself calls internally, so the two are expected -- and verified -- to be
    exactly equal, not merely close."""
    pytest.importorskip("meshioplusplus._core")
    m = _hex_grid(6)
    H = np.array([[2.0, 1.0, 0.3], [1.0, 3.0, 1.0], [0.3, 1.0, 4.0]])
    m.point_data["f"] = _quadratic_scalar(m.points, H)

    cpp_out, cpp_report = mp.hessian(m, "f", return_report=True)
    py_out, skipped, fallback = _hessian_py(m, "f", "green-gauss", "cell", None, False)

    cpp_h = np.asarray(cpp_out.cell_data["f:hessian"][0])
    py_h = np.asarray(py_out.cell_data["f:hessian"][0])
    assert cpp_h.tolist() == py_h.tolist()
    assert skipped == cpp_report["num_skipped"]
    assert fallback == cpp_report["num_fallback"]
