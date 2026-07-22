"""Cross-mesh field transfer (``meshioplusplus.interpolate``)."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._interpolate import _interpolate_py

from . import helpers_data


def _hex_grid(n, shift=0.0):
    """An n x n x n block of unit-cube hexahedra over [shift, n + shift]^3."""
    m = n + 1
    # Point id (k * m + j) * m + i carries coordinates (i, j, k).
    points = np.array(
        [[i, j, k] for k in range(m) for j in range(m) for i in range(m)], dtype=float
    )
    points += shift

    def pid(i, j, k):
        return (k * m + j) * m + i

    cells = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                cells.append(
                    [
                        pid(i, j, k),
                        pid(i + 1, j, k),
                        pid(i + 1, j + 1, k),
                        pid(i, j + 1, k),
                        pid(i, j, k + 1),
                        pid(i + 1, j, k + 1),
                        pid(i + 1, j + 1, k + 1),
                        pid(i, j + 1, k + 1),
                    ]
                )
    return mp.Mesh(points, [("hexahedron", np.array(cells, dtype=np.int64))])


def _quad_grid_2d(n):
    """An n x n grid of unit quads in the z = 0 plane over [0, n]^2."""
    m = n + 1
    jj, ii = np.meshgrid(np.arange(m), np.arange(m), indexing="ij")
    points = np.column_stack([ii.ravel(), jj.ravel()]).astype(float)

    def pid(i, j):
        return j * m + i

    cells = []
    for j in range(n):
        for i in range(n):
            cells.append([pid(i, j), pid(i + 1, j), pid(i + 1, j + 1), pid(i, j + 1)])
    return mp.Mesh(points, [("quad", np.array(cells, dtype=np.int64))])


def _point_cloud(pts):
    return mp.Mesh(np.asarray(pts, dtype=float), [])


def _linear(points):
    p = np.asarray(points, dtype=float)
    z = p[:, 2] if p.shape[1] > 2 else 0.0
    return 2.0 * p[:, 0] + 3.0 * p[:, 1] + 5.0 * z + 7.0


def test_nearest_piecewise_constant():
    src = _point_cloud([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    src.point_data["f"] = np.array([10.0, 20.0, 30.0, 40.0])
    tgt = _point_cloud([[0.1, 0.1, 0], [0.9, 0.2, 0], [0.2, 0.8, 0], [1.4, 1.3, 0]])
    out = mp.interpolate(src, tgt)
    assert np.array_equal(out.point_data["f"], [10.0, 20.0, 30.0, 40.0])


def test_nearest_tie_breaks_to_lowest_index():
    src = _point_cloud([[0, 0, 0], [1, 0, 0]])
    src.point_data["f"] = np.array([10.0, 20.0])
    tgt = _point_cloud([[0.5, 0, 0]])
    out = mp.interpolate(src, tgt)
    assert out.point_data["f"][0] == 10.0


def test_nearest_preserves_dtype():
    src = _point_cloud([[0, 0, 0], [1, 0, 0]])
    src.point_data["tag"] = np.array([3, 9], dtype=np.int32)
    tgt = _point_cloud([[0.1, 0, 0], [0.9, 0, 0]])
    out = mp.interpolate(src, tgt)
    assert out.point_data["tag"].dtype.kind == "i"
    assert np.array_equal(out.point_data["tag"], [3, 9])


def test_barycentric_linear_exact():
    src = _hex_grid(2)
    src.point_data["f"] = _linear(src.points)
    probes = np.array(
        [
            [0.25, 0.5, 0.75],
            [1.1, 0.2, 1.9],
            [0.0, 0.0, 0.0],
            [2.0, 2.0, 2.0],
            [1.0, 1.0, 1.0],
            [1.5, 0.5, 0.5],
            [0.7, 1.3, 0.2],
        ]
    )
    tgt = _point_cloud(probes)
    out = mp.interpolate(src, tgt, method="barycentric")
    assert out.point_data["f"].dtype == np.float64
    np.testing.assert_allclose(out.point_data["f"], _linear(probes), rtol=0, atol=1e-12)


def test_barycentric_triangle_source_2d():
    src = _quad_grid_2d(3)
    src.point_data["f"] = _linear(src.points)
    probes = np.array([[0.5, 0.5], [2.7, 1.2], [0.0, 3.0], [1.0, 1.0]])
    tgt = _point_cloud(probes)
    out = mp.interpolate(src, tgt, method="barycentric")
    np.testing.assert_allclose(out.point_data["f"], _linear(probes), rtol=0, atol=1e-12)


def test_barycentric_vector_components():
    src = _hex_grid(2)
    p = np.asarray(src.points, dtype=float)
    src.point_data["v"] = np.column_stack(
        [p[:, 0] + 2 * p[:, 1], 3 * p[:, 0] - p[:, 1]]
    )
    tgt = _point_cloud([[0.3, 1.4, 0.6]])
    out = mp.interpolate(src, tgt, method="barycentric")
    np.testing.assert_allclose(
        out.point_data["v"], [[0.3 + 2 * 1.4, 3 * 0.3 - 1.4]], rtol=0, atol=1e-12
    )


def test_outside_default_value():
    src = _hex_grid(1)
    src.point_data["f"] = _linear(src.points)
    tgt = _point_cloud([[0.5, 0.5, 0.5], [5.0, 5.0, 5.0]])
    out = mp.interpolate(src, tgt, method="barycentric", default_value=-123.0)
    np.testing.assert_allclose(
        out.point_data["f"][0], _linear([[0.5, 0.5, 0.5]])[0], rtol=0, atol=1e-12
    )
    assert out.point_data["f"][1] == -123.0


def test_outside_extrapolate():
    src = _hex_grid(1)
    src.point_data["f"] = _linear(src.points)
    tgt = _point_cloud([[5.0, 5.0, 5.0]])
    out = mp.interpolate(src, tgt, method="barycentric", extrapolate=True)
    # Nearest source point to (5,5,5) is the cube corner (1,1,1).
    assert out.point_data["f"][0] == _linear([[1.0, 1.0, 1.0]])[0]


def test_cell_data_nearest_cell():
    src = mp.Mesh(
        np.array(
            [[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0]],
            dtype=float,
        ),
        [("quad", np.array([[0, 1, 4, 3], [1, 2, 5, 4]]))],
    )
    src.cell_data["mat"] = [np.array([1.0, 2.0])]
    tgt = mp.Mesh(
        np.array(
            [[1.2, 0.2, 0], [1.8, 0.2, 0], [1.8, 0.8, 0], [1.2, 0.8, 0]], dtype=float
        ),
        [("quad", np.array([[0, 1, 2, 3]]))],
    )
    for method in ("nearest", "barycentric"):
        out = mp.interpolate(src, tgt, method=method, arrays=["mat"])
        assert np.array_equal(out.cell_data["mat"][0], [2.0])


def test_default_arrays_is_point_data_only():
    src = _quad_grid_2d(1)
    src.point_data["f"] = _linear(src.points)
    src.cell_data["mat"] = [np.array([4.0])]
    out = mp.interpolate(src, _quad_grid_2d(1))
    assert "f" in out.point_data
    assert "mat" not in out.cell_data


def test_explicit_array_in_both_locations_transfers_both():
    src = _quad_grid_2d(1)
    src.point_data["m"] = _linear(src.points)
    src.cell_data["m"] = [np.array([4.0])]
    out = mp.interpolate(src, _quad_grid_2d(1), arrays=["m"])
    assert "m" in out.point_data
    assert "m" in out.cell_data


def test_on_conflict_behaviours():
    src = _point_cloud([[0, 0, 0]])
    src.point_data["f"] = np.array([99.0])
    tgt = _point_cloud([[0.25, 0, 0]])
    tgt.point_data["f"] = np.array([1.0])

    with pytest.raises(ValueError):
        mp.interpolate(src, tgt)

    out = mp.interpolate(src, tgt, on_conflict="overwrite")
    assert out.point_data["f"][0] == 99.0

    out = mp.interpolate(src, tgt, on_conflict="suffix")
    assert out.point_data["f"][0] == 1.0
    assert out.point_data["f_interp"][0] == 99.0

    tgt.point_data["f_interp"] = np.array([2.0])
    with pytest.raises(ValueError):
        mp.interpolate(src, tgt, on_conflict="suffix")


def test_target_geometry_and_own_data_unchanged():
    src = _quad_grid_2d(2)
    src.point_data["srcfield"] = _linear(src.points)
    tgt = helpers_data.data_mesh()
    tgt.point_sets["left"] = np.array([0, 3])
    tgt.cell_sets["tris"] = [np.array([0, 1]), np.array([], dtype=int)]

    out = mp.interpolate(src, tgt)
    helpers_data.assert_same_geometry(tgt, out)
    assert "srcfield" in out.point_data
    for name, arr in tgt.point_data.items():
        assert np.array_equal(out.point_data[name], arr)
    for name, blocks in tgt.cell_data.items():
        for got, want in zip(out.cell_data[name], blocks):
            assert np.array_equal(got, want)
            assert got.dtype == want.dtype
    assert "meta" in out.field_data
    assert np.array_equal(out.point_sets["left"], [0, 3])
    assert np.array_equal(out.cell_sets["tris"][0], [0, 1])


def test_source_sets_not_transferred():
    src = _point_cloud([[0, 0, 0]])
    src.point_data["f"] = np.array([1.0])
    src.point_sets["marked"] = np.array([0])
    out = mp.interpolate(src, _point_cloud([[0, 0, 0]]))
    assert "marked" not in out.point_sets


def test_unknown_array_raises():
    src = _point_cloud([[0, 0, 0]])
    src.point_data["f"] = np.array([1.0])
    with pytest.raises(ValueError, match="nope"):
        mp.interpolate(src, _point_cloud([[0, 0, 0]]), arrays=["nope"])


def test_bad_method_and_conflict_names():
    src = _point_cloud([[0, 0, 0]])
    src.point_data["f"] = np.array([1.0])
    tgt = _point_cloud([[0, 0, 0]])
    with pytest.raises(ValueError):
        mp.interpolate(src, tgt, method="bogus")
    with pytest.raises(ValueError):
        mp.interpolate(src, tgt, on_conflict="bogus")


def test_empty_source_raises():
    src = mp.Mesh(np.zeros((0, 3)), [])
    with pytest.raises(ValueError):
        mp.interpolate(src, _point_cloud([[0, 0, 0]]))


def test_barycentric_without_cells_raises():
    src = _point_cloud([[0, 0, 0], [1, 0, 0]])
    src.point_data["f"] = np.array([1.0, 2.0])
    with pytest.raises(ValueError):
        mp.interpolate(src, _point_cloud([[0.5, 0, 0]]), method="barycentric")


def test_barycentric_integer_data_becomes_float64():
    src = _hex_grid(1)
    src.point_data["tag"] = np.arange(8, dtype=np.int32)
    out = mp.interpolate(src, _point_cloud([[0.5, 0.5, 0.5]]), method="barycentric")
    assert out.point_data["tag"].dtype == np.float64


@pytest.mark.parametrize("method", ["nearest", "barycentric"])
def test_cpp_matches_python(method):
    """The C++ core and the numpy fallback must agree byte for byte."""
    core = pytest.importorskip("meshioplusplus._core")

    src = _hex_grid(3)
    src.point_data["f"] = _linear(src.points)
    p = np.asarray(src.points, dtype=float)
    src.point_data["v"] = np.column_stack([p[:, 0] - p[:, 2], p[:, 1] * 0.5])
    src.cell_data["mat"] = [np.arange(27, dtype=np.int32)]
    # Off-grid shift plus points outside the source (the extrapolate branch).
    tgt = _hex_grid(2, shift=0.37)
    tgt.points = np.vstack([np.asarray(tgt.points), [[9.0, -3.0, 11.0]]])

    for extrapolate in (False, True):
        got = core.interpolate(
            src, tgt, method, ["f", "v", "mat"], extrapolate, -7.5, "error"
        )
        want = _interpolate_py(
            src,
            tgt,
            method=method,
            arrays=["f", "v", "mat"],
            extrapolate=extrapolate,
            default_value=-7.5,
            on_conflict="error",
        )
        for name in ("f", "v"):
            a = np.asarray(got.point_data[name])
            b = np.asarray(want.point_data[name])
            assert a.dtype == b.dtype
            assert np.array_equal(a, b), (method, extrapolate, name)
        for a, b in zip(got.cell_data["mat"], want.cell_data["mat"]):
            assert np.asarray(a).dtype == np.asarray(b).dtype
            assert np.array_equal(a, b)


def test_determinism_two_runs():
    src = _hex_grid(3)
    src.point_data["f"] = _linear(src.points)
    tgt = _hex_grid(2, shift=0.37)
    for method in ("nearest", "barycentric"):
        a = mp.interpolate(src, tgt, method=method)
        b = mp.interpolate(src, tgt, method=method)
        assert a.point_data["f"].tobytes() == b.point_data["f"].tobytes()
