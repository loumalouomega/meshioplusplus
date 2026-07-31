"""Tests for the field differential operators ``meshioplusplus.gradient``.

Several "obvious" tests are inert for this operation and are deliberately absent
or strengthened -- the fixtures here exist to defeat exactly those traps:

* "the gradient of a constant field is zero" proves nothing: the numerator is
  ``f_bar * sum(A_j)`` and ``sum(A_j) == 0`` over **any** closed surface with
  **any** quadrature weights and **any** volume.
* a cube or an axis-aligned hex grid cannot distinguish the face fan from a naive
  corner-average face value, because a parallelogram's corner average **is** its
  area centroid. Hence :func:`_tapered_hex` (planar trapezoidal faces) and
  :func:`_warped_hex` (a genuinely non-planar face).
* an all-tetra suite cannot see any quad-face bug, since a triangle's fan is
  redundant.
* a solenoidal field cannot see a swapped divergence/curl, and a curl of the form
  ``(0, 0, c)`` cannot see an index permutation. The vector fixtures use three
  distinct nonzero answers.
* ``test_cpp_matches_python`` pins **parity only** -- both sides can be
  identically wrong -- so it never substitutes for the analytic assertions.
"""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._gradient import _gradient_py

from .helpers_data import assert_same_geometry

# --------------------------------------------------------------------------- #
# geometry builders with hand-computable answers                              #
# --------------------------------------------------------------------------- #


def _tapered_hex_points():
    """A frustum: a 2x2 base at z=0 tapering to a 1x1 top at z=1.

    Every side face is a **planar trapezoid**, whose corner average is not its
    area centroid -- the fixture a naive (no-fan) face value fails.
    """
    return np.array(
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


def _warped_hex_points():
    """A cube with one top corner lifted, so face (4,5,6,7) is non-planar."""
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.4],
            [0.0, 1.0, 1.0],
        ]
    )


def _one_hex(points):
    return mp.Mesh(points, [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))])


def _tapered_hex():
    return _one_hex(_tapered_hex_points())


def _warped_hex():
    return _one_hex(_warped_hex_points())


def _hex_grid(n=3, off=0.0):
    """A structured n x n x n hexahedron grid, optionally translated by ``off``."""
    pts = []
    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                pts.append([i + off, j + off, k + off])

    def vid(i, j, k):
        return (k * (n + 1) + j) * (n + 1) + i

    cells = [
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
        for k in range(n)
        for j in range(n)
        for i in range(n)
    ]
    return mp.Mesh(np.array(pts, dtype=float), [("hexahedron", np.array(cells))])


def _quad_sheet():
    """A planar trapezoid plus a triangle-ish quad, in the xy-plane."""
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [5.0, 0.0, 0.0],
            [4.0, 1.0, 0.0],
        ]
    )
    return mp.Mesh(pts, [("quad", np.array([[0, 1, 2, 3], [1, 4, 5, 2]]))])


def _linear(mesh, a, b, c, d=0.0):
    p = np.asarray(mesh.points, dtype=float)
    return a * p[:, 0] + b * p[:, 1] + c * p[:, 2] + d


def _linear_vector(mesh, matrix):
    p = np.asarray(mesh.points, dtype=float)
    return p @ np.asarray(matrix, dtype=float).T


# --------------------------------------------------------------------------- #
# the headline invariant                                                       #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("method", ["green-gauss", "least-squares"])
@pytest.mark.parametrize("location", ["cell", "point"])
def test_linear_field_gradient_is_exact(method, location):
    """The reason this operation exists: grad of ax+by+cz+d is exactly (a,b,c).

    A 3x3x3 grid gives least-squares a full-rank neighbourhood everywhere; the
    tapered/warped single-cell fixtures below cover the Green-Gauss quadrature.
    """
    mesh = _hex_grid(3)
    mesh.point_data["f"] = _linear(mesh, 3.0, -2.0, 5.0, 7.0)
    out, report = mp.gradient(
        mesh, "f", method=method, location=location, return_report=True
    )
    assert report["num_skipped"] == 0
    if location == "point":
        got = np.asarray(out.point_data["f:gradient"])
    else:
        got = np.concatenate([np.asarray(a) for a in out.cell_data["f:gradient"]])
    assert got.shape[1] == 3
    # Not `==`: for the point location, summing n copies of g and dividing by n
    # is not bit-exact in IEEE even though every contributing cell value is.
    np.testing.assert_allclose(
        got, np.tile([3.0, -2.0, 5.0], (got.shape[0], 1)), atol=1e-11
    )


@pytest.mark.parametrize("builder", [_tapered_hex, _warped_hex])
def test_green_gauss_is_exact_on_trapezoidal_and_warped_faces(builder):
    """Exactness does not depend on face planarity.

    Two faces sharing an edge contribute oppositely-wound fan triangles there, so
    the fan surface is closed and the divergence theorem applies verbatim. A cube
    fixture would pass even with a broken quadrature -- these two do not.
    """
    mesh = builder()
    mesh.point_data["f"] = _linear(mesh, 3.0, -2.0, 5.0, 7.0)
    out = mp.gradient(mesh, "f")
    np.testing.assert_allclose(
        np.asarray(out.cell_data["f:gradient"][0])[0], [3.0, -2.0, 5.0], atol=1e-12
    )


def test_linear_field_on_a_planar_trapezoid_is_exact_in_plane():
    mesh = _quad_sheet()
    mesh.point_data["f"] = _linear(mesh, 1.25, -0.75, 0.0, 9.0)
    out = mp.gradient(mesh, "f")
    got = np.asarray(out.cell_data["f:gradient"][0])
    np.testing.assert_allclose(got[0], [1.25, -0.75, 0.0], atol=1e-12)
    # The z component is projected onto the cell's own plane, so it is exactly 0.
    assert got[0, 2] == 0.0


@pytest.mark.parametrize("method", ["green-gauss", "least-squares"])
def test_translating_the_mesh_does_not_move_the_gradient(method):
    """The gate for recentring.

    ``V = (1/3) sum x_j . A_j`` only telescopes because ``sum A_j == 0``; without
    recentring the terms are ~off*area and cancel to ~area, costing log10(off)
    digits of a quantity we then divide by. Verified inert at off=1e6/atol=1e-8
    and firing at 1e8/1e-9.
    """
    off = 1.0e8
    mesh = _hex_grid(2, off=off)
    mesh.point_data["f"] = _linear(mesh, 3.0, -2.0, 5.0)
    out = mp.gradient(mesh, "f", method=method)
    got = np.asarray(out.cell_data["f:gradient"][0])
    np.testing.assert_allclose(
        got, np.tile([3.0, -2.0, 5.0], (got.shape[0], 1)), atol=1e-9
    )


# --------------------------------------------------------------------------- #
# vector operators                                                             #
# --------------------------------------------------------------------------- #


def test_divergence_of_an_anisotropic_field():
    """u = (2x, 3y, 5z) -> 10, from three DISTINCT nonzero diagonal terms.

    A solenoidal field (div = 0) would be inert against a swapped diagonal.
    """
    mesh = _tapered_hex()
    mesh.point_data["u"] = _linear_vector(mesh, [[2, 0, 0], [0, 3, 0], [0, 0, 5]])
    out = mp.gradient(mesh, "u", operator="divergence")
    got = np.asarray(out.cell_data["u:divergence"][0])
    assert got.shape == (1,)
    np.testing.assert_allclose(got[0], 10.0, atol=1e-12)


def test_curl_has_three_distinct_nonzero_components():
    """u = (7z, 11x, 13y) -> curl = (13, 7, 11).

    Distinct and nonzero in every slot, so any index permutation or sign flip in
    the curl formula fails.
    """
    mesh = _warped_hex()
    mesh.point_data["u"] = _linear_vector(mesh, [[0, 0, 7], [11, 0, 0], [0, 13, 0]])
    out = mp.gradient(mesh, "u", operator="curl")
    np.testing.assert_allclose(
        np.asarray(out.cell_data["u:curl"][0])[0], [13.0, 7.0, 11.0], atol=1e-12
    )


def test_tensor_layout_is_component_major_then_derivative():
    """All nine entries distinct, so the i*3+j <-> j*3+i transpose fails."""
    c = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
    mesh = _tapered_hex()
    mesh.point_data["u"] = _linear_vector(mesh, c)
    out = mp.gradient(mesh, "u")
    got = np.asarray(out.cell_data["u:gradient"][0])
    assert got.shape == (1, 9)
    np.testing.assert_allclose(got[0], c.reshape(9), atol=1e-12)


def test_component_selection_yields_three_components():
    c = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
    mesh = _tapered_hex()
    mesh.point_data["u"] = _linear_vector(mesh, c)
    out = mp.gradient(mesh, "u", component=1)
    got = np.asarray(out.cell_data["u:gradient"][0])
    assert got.shape == (1, 3)
    np.testing.assert_allclose(got[0], [4.0, 5.0, 6.0], atol=1e-12)


def test_two_component_field_is_read_as_u_v_zero():
    """A 2-component field pads to (u, v, 0), like 2-D point coordinates."""
    mesh = _quad_sheet()
    p = np.asarray(mesh.points, dtype=float)
    mesh.point_data["w"] = np.stack([2.0 * p[:, 0], 3.0 * p[:, 1]], axis=1)
    div = mp.gradient(mesh, "w", operator="divergence")
    np.testing.assert_allclose(
        np.asarray(div.cell_data["w:divergence"][0])[0], 5.0, atol=1e-12
    )
    # Only the z component of the curl is nonzero for a 2-component field.
    curl = mp.gradient(mesh, "w", operator="curl")
    got = np.asarray(curl.cell_data["w:curl"][0])[0]
    assert got[0] == 0.0 and got[1] == 0.0


def test_identities_curl_of_gradient_and_divergence_of_curl_vanish():
    mesh = _hex_grid(3)
    mesh.point_data["f"] = _linear(mesh, 2.0, -3.0, 4.0)
    mesh.point_data["u"] = _linear_vector(mesh, [[0, 0, 7], [11, 0, 0], [0, 13, 0]])

    with_grad = mp.gradient(mesh, "f", location="point")
    cg = mp.gradient(with_grad, "f:gradient", operator="curl")
    np.testing.assert_allclose(
        np.asarray(cg.cell_data["f:gradient:curl"][0]), 0.0, atol=1e-9
    )

    with_curl = mp.gradient(mesh, "u", operator="curl", location="point")
    dc = mp.gradient(with_curl, "u:curl", operator="divergence")
    np.testing.assert_allclose(
        np.asarray(dc.cell_data["u:curl:divergence"][0]), 0.0, atol=1e-9
    )


# --------------------------------------------------------------------------- #
# structural invariants                                                        #
# --------------------------------------------------------------------------- #


def test_an_inverted_cell_matches_its_positively_wound_twin():
    """Numerator and denominator flip together, so the answers must be equal --
    not merely finite, which would be inert."""
    pts = _warped_hex_points()
    f = 3.0 * pts[:, 0] - 2.0 * pts[:, 1] + 5.0 * pts[:, 2]
    good = mp.Mesh(pts, [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))])
    good.point_data["f"] = f
    flipped = mp.Mesh(pts, [("hexahedron", np.array([[4, 5, 6, 7, 0, 1, 2, 3]]))])
    flipped.point_data["f"] = f
    np.testing.assert_allclose(
        np.asarray(mp.gradient(good, "f").cell_data["f:gradient"][0]),
        np.asarray(mp.gradient(flipped, "f").cell_data["f:gradient"][0]),
        atol=1e-12,
    )


def test_reversing_a_quad_ring_leaves_the_gradient_unchanged():
    """A triangle's reversal is nearly symmetric; a quad exercises ring order."""
    pts = np.array([[0.0, 0, 0], [4.0, 0, 0], [2.0, 1, 0], [0.0, 1, 0]])
    f = 1.25 * pts[:, 0] - 0.75 * pts[:, 1] + 9.0
    fwd = mp.Mesh(pts, [("quad", np.array([[0, 1, 2, 3]]))])
    fwd.point_data["f"] = f
    rev = mp.Mesh(pts, [("quad", np.array([[0, 3, 2, 1]]))])
    rev.point_data["f"] = f
    np.testing.assert_allclose(
        np.asarray(mp.gradient(fwd, "f").cell_data["f:gradient"][0]),
        np.asarray(mp.gradient(rev, "f").cell_data["f:gradient"][0]),
        atol=1e-12,
    )


def test_geometry_regions_and_other_data_are_untouched():
    mesh = _hex_grid(2)
    mesh.point_data["f"] = _linear(mesh, 1.0, 2.0, 3.0)
    mesh.point_data["other"] = np.arange(len(mesh.points), dtype=np.float64)
    mesh.cell_data["mat"] = [np.arange(len(mesh.cells[0].data), dtype=np.int32)]
    mesh.field_data["meta"] = np.array([1.0, 2.0])
    mesh.point_sets["inlet"] = np.array([0, 1, 2])

    out = mp.gradient(mesh, "f")
    assert_same_geometry(mesh, out)
    np.testing.assert_array_equal(
        np.asarray(out.point_data["other"]), np.asarray(mesh.point_data["other"])
    )
    np.testing.assert_array_equal(
        np.asarray(out.cell_data["mat"][0]), np.asarray(mesh.cell_data["mat"][0])
    )
    assert np.asarray(out.cell_data["mat"][0]).dtype == np.int32
    np.testing.assert_array_equal(
        np.asarray(out.field_data["meta"]), np.asarray(mesh.field_data["meta"])
    )
    np.testing.assert_array_equal(np.asarray(out.point_sets["inlet"]), [0, 1, 2])


def test_output_is_always_float64():
    mesh = _hex_grid(2)
    mesh.point_data["f"] = _linear(mesh, 1.0, 2.0, 3.0).astype(np.float32)
    out = mp.gradient(mesh, "f")
    assert np.asarray(out.cell_data["f:gradient"][0]).dtype == np.float64


def test_output_name_can_be_overridden():
    mesh = _hex_grid(2)
    mesh.point_data["f"] = _linear(mesh, 1.0, 2.0, 3.0)
    out = mp.gradient(mesh, "f", output="dTdx")
    assert "dTdx" in out.cell_data
    assert "f:gradient" not in out.cell_data


# --------------------------------------------------------------------------- #
# the gates must actually fire                                                 #
# --------------------------------------------------------------------------- #


def test_least_squares_falls_back_on_a_degenerate_neighbourhood():
    """Asserting ``num_fallback == 0`` on a nice mesh would be inert -- it passes
    if the conditioning test never fires. A lone cell has no neighbours at all,
    so ``M == 0``; the fallback must produce the Green-Gauss answer, not NaN."""
    mesh = _tapered_hex()
    mesh.point_data["f"] = _linear(mesh, 3.0, -2.0, 5.0, 7.0)
    lsq, report = mp.gradient(mesh, "f", method="least-squares", return_report=True)
    assert report["num_fallback"] == 1
    assert report["num_skipped"] == 0
    gg = mp.gradient(mesh, "f", method="green-gauss")
    np.testing.assert_array_equal(
        np.asarray(lsq.cell_data["f:gradient"][0]),
        np.asarray(gg.cell_data["f:gradient"][0]),
    )


def test_least_squares_falls_back_on_a_collinear_strip():
    """A row of cells has a rank-1 normal matrix in 3-D."""
    mesh = _hex_grid(1)  # a single cell would be the lone-cell case
    pts = []
    for k in range(2):
        for j in range(2):
            for i in range(5):
                pts.append([i, j, k])
    pts = np.array(pts, dtype=float)

    def vid(i, j, k):
        return (k * 2 + j) * 5 + i

    cells = [
        [
            vid(i, 0, 0),
            vid(i + 1, 0, 0),
            vid(i + 1, 1, 0),
            vid(i, 1, 0),
            vid(i, 0, 1),
            vid(i + 1, 0, 1),
            vid(i + 1, 1, 1),
            vid(i, 1, 1),
        ]
        for i in range(4)
    ]
    mesh = mp.Mesh(pts, [("hexahedron", np.array(cells))])
    mesh.point_data["f"] = pts[:, 0]
    _out, report = mp.gradient(mesh, "f", method="least-squares", return_report=True)
    assert report["num_fallback"] > 0


def test_unsupported_cells_are_nan_and_counted():
    """Boundary triangles on a tet mesh are below the mesh's dimension."""
    mesh = mp.Mesh(
        np.array([[0.0, 0, 0], [1.0, 0, 0], [0.0, 1, 0], [0.0, 0, 1]]),
        [
            ("tetra", np.array([[0, 1, 2, 3]])),
            ("triangle", np.array([[0, 1, 2], [0, 1, 3]])),
        ],
    )
    mesh.point_data["f"] = np.array([0.0, 1.0, 2.0, 3.0])
    out, report = mp.gradient(mesh, "f", return_report=True)
    assert report["num_skipped"] == 2
    assert np.all(np.isnan(np.asarray(out.cell_data["f:gradient"][1])))
    assert not np.any(np.isnan(np.asarray(out.cell_data["f:gradient"][0])))


def test_ragged_blocks_are_skipped_not_guessed():
    mesh = mp.Mesh(
        np.array([[0.0, 0, 0], [1.0, 0, 0], [1.0, 1, 0], [0.0, 1, 0], [2.0, 0.5, 0]]),
        [("polygon", [[0, 1, 2, 3], [1, 4, 2]])],
    )
    mesh.point_data["f"] = np.arange(5, dtype=float)
    out, report = mp.gradient(mesh, "f", return_report=True)
    assert report["num_skipped"] == 2
    assert np.all(np.isnan(np.asarray(out.cell_data["f:gradient"][0])))


# --------------------------------------------------------------------------- #
# errors                                                                       #
# --------------------------------------------------------------------------- #


def test_cell_data_field_raises_by_name():
    mesh = _hex_grid(2)
    mesh.cell_data["mat"] = [np.arange(8, dtype=float)]
    with pytest.raises(ValueError) as exc:
        mp.gradient(mesh, "mat")
    msg = str(exc.value)
    assert "to-point" in msg or "cell_data_to_point_data" in msg


def test_unknown_array_lists_what_exists():
    mesh = _hex_grid(2)
    mesh.point_data["temperature"] = _linear(mesh, 1.0, 0.0, 0.0)
    with pytest.raises(ValueError) as exc:
        mp.gradient(mesh, "nope")
    assert "temperature" in str(exc.value)


def test_bad_arguments_raise():
    mesh = _hex_grid(2)
    mesh.point_data["f"] = _linear(mesh, 1.0, 0.0, 0.0)
    mesh.point_data["u"] = _linear_vector(mesh, np.eye(3))

    with pytest.raises(ValueError):
        mp.gradient(mesh, "f", operator="laplacian")
    with pytest.raises(ValueError):
        mp.gradient(mesh, "f", method="magic")
    with pytest.raises(ValueError):
        mp.gradient(mesh, "f", location="field")
    with pytest.raises(ValueError):
        mp.gradient(mesh, "u", component=7)
    with pytest.raises(ValueError):
        # A component selection is meaningless for divergence.
        mp.gradient(mesh, "u", operator="divergence", component=0)
    with pytest.raises(ValueError):
        # A scalar has no divergence.
        mp.gradient(mesh, "f", operator="divergence")


def test_output_name_collision_is_rejected_unless_overwriting():
    mesh = _hex_grid(2)
    mesh.point_data["f"] = _linear(mesh, 1.0, 0.0, 0.0)
    mesh.cell_data["taken"] = [np.zeros(8)]
    with pytest.raises(ValueError):
        mp.gradient(mesh, "f", output="taken")
    out = mp.gradient(mesh, "f", output="taken", overwrite=True)
    assert "taken" in out.cell_data


@pytest.mark.parametrize(
    "kwargs",
    [
        {"array": "mat"},  # cell_data
        {"array": "nope"},  # unknown
        {"array": "u", "component": 7},  # out of range
        {"array": "f", "operator": "divergence"},  # scalar has no divergence
        {"array": "f", "output": "mat"},  # collision without overwrite
    ],
)
def test_the_error_surface_does_not_depend_on_the_core(kwargs):
    """Validation happens before either path runs.

    The shim re-raises ValueError/TypeError from the core rather than falling
    through to numpy, so a user error must be rejected identically whether or not
    the extension is built. This drives the shared validator through the public
    entry point and asserts the numpy reference rejects the same inputs.
    """
    mesh = _hex_grid(2)
    mesh.point_data["f"] = _linear(mesh, 1.0, 0.0, 0.0)
    mesh.point_data["u"] = _linear_vector(mesh, np.eye(3))
    mesh.cell_data["mat"] = [np.arange(8, dtype=float)]
    with pytest.raises(ValueError):
        mp.gradient(mesh, **kwargs)


# --------------------------------------------------------------------------- #
# C++ / numpy byte-parity and determinism                                      #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("method", ["green-gauss", "least-squares"])
@pytest.mark.parametrize("location", ["cell", "point"])
@pytest.mark.parametrize(
    "array,operator,component",
    [
        ("f", "gradient", None),  # scalar
        ("r2", "gradient", None),  # NON-linear: the two paths could easily diverge
        ("u", "gradient", None),  # (n,9) tensor
        ("u", "gradient", 1),  # one component
        ("u", "divergence", None),
        ("u", "curl", None),
    ],
)
def test_cpp_matches_python(method, location, array, operator, component):
    """Byte-identity across the C++/numpy boundary.

    This pins **parity only** -- both sides can be identically wrong -- so it
    never substitutes for the analytic tests above.
    """
    core = pytest.importorskip("meshioplusplus._core")
    mesh = _hex_grid(3)
    p = np.asarray(mesh.points, dtype=float)
    mesh.point_data["f"] = 0.5 * p[:, 0] - 1.25 * p[:, 1] + 2.0 * p[:, 2] + 4.0
    mesh.point_data["r2"] = (p**2).sum(axis=1)
    mesh.point_data["u"] = (
        p @ np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float).T
    )

    res = core.gradient(
        mesh,
        array,
        operator,
        method,
        location,
        "",
        -1 if component is None else component,
        False,
    )
    got = res["mesh"]
    ref, skipped, fell = _gradient_py(
        mesh,
        array,
        operator=operator,
        method=method,
        location=location,
        component=component,
    )

    assert (res["num_skipped"], res["num_fallback"]) == (skipped, fell)
    assert np.array_equal(np.asarray(got.points), np.asarray(ref.points))
    assert [cb.type for cb in got.cells] == [cb.type for cb in ref.cells]
    for ca, cb in zip(got.cells, ref.cells):
        assert np.array_equal(np.asarray(ca.data), np.asarray(cb.data))
    assert set(got.point_data) == set(ref.point_data)
    for key in got.point_data:
        a, b = np.asarray(got.point_data[key]), np.asarray(ref.point_data[key])
        assert a.dtype == b.dtype and a.shape == b.shape
        assert a.tobytes() == b.tobytes(), key
    assert set(got.cell_data) == set(ref.cell_data)
    for key in got.cell_data:
        for a, b in zip(got.cell_data[key], ref.cell_data[key]):
            a, b = np.asarray(a), np.asarray(b)
            assert a.dtype == b.dtype and a.shape == b.shape
            assert a.tobytes() == b.tobytes(), key


def test_cpp_matches_python_on_skipped_and_fallback_paths():
    """The parity cases above never exercise a NaN row or a fallback."""
    core = pytest.importorskip("meshioplusplus._core")
    mesh = mp.Mesh(
        np.array([[0.0, 0, 0], [1.0, 0, 0], [0.0, 1, 0], [0.0, 0, 1]]),
        [
            ("tetra", np.array([[0, 1, 2, 3]])),
            ("triangle", np.array([[0, 1, 2], [0, 1, 3]])),
        ],
    )
    mesh.point_data["f"] = np.array([0.0, 1.0, 2.0, 3.0])
    for method in ("green-gauss", "least-squares"):
        res = core.gradient(mesh, "f", "gradient", method, "cell", "", -1, False)
        ref, skipped, fell = _gradient_py(mesh, "f", method=method)
        assert (res["num_skipped"], res["num_fallback"]) == (skipped, fell)
        for a, b in zip(
            res["mesh"].cell_data["f:gradient"], ref.cell_data["f:gradient"]
        ):
            a, b = np.asarray(a), np.asarray(b)
            assert a.dtype == b.dtype and a.shape == b.shape
            # NaN != NaN, so compare the raw bytes rather than the values.
            assert a.tobytes() == b.tobytes()


@pytest.mark.parametrize("method", ["green-gauss", "least-squares"])
def test_determinism_two_runs(method):
    mesh = _hex_grid(3)
    mesh.point_data["f"] = _linear(mesh, 0.5, -1.25, 2.0, 4.0)
    a = mp.gradient(mesh, "f", method=method)
    b = mp.gradient(mesh, "f", method=method)
    assert (
        np.asarray(a.cell_data["f:gradient"][0]).tobytes()
        == np.asarray(b.cell_data["f:gradient"][0]).tobytes()
    )


# --------------------------------------------------------------------------- #
# the compositions this operation exists for                                   #
# --------------------------------------------------------------------------- #


def test_composes_with_isosurface_to_contour_a_derived_quantity():
    """gradient -> data to-point -> isosurface on |grad f|."""
    mesh = _hex_grid(4)
    p = np.asarray(mesh.points, dtype=float)
    mesh.point_data["T"] = (p**2).sum(axis=1)
    g = mp.gradient(mesh, "T", location="point")
    grad = np.asarray(g.point_data["T:gradient"])
    g.point_data["gmag"] = np.sqrt((grad**2).sum(axis=1))
    contour = mp.isosurface(g, "gmag", [4.0])
    assert len(contour.points) > 0
    assert sum(len(cb.data) for cb in contour.cells) > 0


def test_composes_with_refine_as_an_error_indicator():
    """gradient -> refine(where=...) : gradient-magnitude driven adaptivity."""
    mesh = _hex_grid(3)
    p = np.asarray(mesh.points, dtype=float)
    mesh.point_data["T"] = p[:, 0] ** 2
    g = mp.gradient(mesh, "T")
    grad = [np.asarray(a) for a in g.cell_data["T:gradient"]]
    g.cell_data["gmag"] = [np.sqrt((a**2).sum(axis=1)) for a in grad]
    before = sum(len(cb.data) for cb in g.cells)
    refined = mp.refine(g, where="gmag > 3.0")
    after = sum(len(cb.data) for cb in refined.cells)
    assert after > before, "the indicator selected no cells"
