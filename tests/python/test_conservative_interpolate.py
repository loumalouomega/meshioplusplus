"""Tests for mass-preserving cross-mesh field transfer
(the ``conservative_interpolate`` operation)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import conservative_interpolate

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(
    _core is None, reason="conservative_interpolate has no pure-Python fallback"
)


def _quad_grid(n, scale):
    """An n x n grid of quads over [0, n*scale]^2."""
    m = n + 1
    pts = np.array(
        [[i * scale, j * scale, 0.0] for j in range(m) for i in range(m)], float
    )

    def pid(i, j):
        return j * m + i

    cells = np.array(
        [
            [pid(i, j), pid(i + 1, j), pid(i + 1, j + 1), pid(i, j + 1)]
            for j in range(n)
            for i in range(n)
        ]
    )
    return meshioplusplus.Mesh(pts, [("quad", cells)])


def _hex_grid(n, scale):
    """An n x n x n grid of hexahedra over [0, n*scale]^3."""
    m = n + 1
    pts = []
    for k in range(m):
        for j in range(m):
            for i in range(m):
                pts.append([i * scale, j * scale, k * scale])
    pts = np.array(pts, float)

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
    return meshioplusplus.Mesh(pts, [("hexahedron", np.array(cells))])


def _quad_weighted_sum(mesh, name):
    total = 0.0
    for block, data in zip(mesh.cells, mesh.cell_data[name]):
        conn = np.asarray(block.data)
        pts = np.asarray(mesh.points)
        for c in range(len(conn)):
            xs = pts[conn[c], 0]
            ys = pts[conn[c], 1]
            a2 = 0.0
            for k in range(4):
                k2 = (k + 1) % 4
                a2 += xs[k] * ys[k2] - xs[k2] * ys[k]
            total += abs(a2) / 2.0 * np.asarray(data)[c]
    return total


def _one_triangle(p0, p1, p2):
    return meshioplusplus.Mesh(
        np.array([p0, p1, p2], float), [("triangle", np.array([[0, 1, 2]]))]
    )


def _one_polyhedron_tet(p0, p1, p2, p3):
    """A single tetrahedron stored as a polyhedron block, to exercise the
    "ragged/polyhedron accepted for free" scope decision."""
    pts = np.array([p0, p1, p2, p3], float)
    faces = [[0, 1, 3], [1, 2, 3], [2, 0, 3], [0, 2, 1]]
    return meshioplusplus.Mesh(pts, [("polyhedron4", [faces])])


@needs_core
def test_cell_data_conserves_total_integral_across_different_partitions():
    src = _quad_grid(3, 2.0)
    tgt = _quad_grid(4, 1.5)
    src.cell_data["f"] = [
        np.array([1.0 + 2 * i + 3 * j for j in range(3) for i in range(3)])
    ]

    out = conservative_interpolate(src, tgt)
    assert _quad_weighted_sum(out, "f") == pytest.approx(
        _quad_weighted_sum(src, "f"), abs=1e-8
    )


@needs_core
def test_point_data_is_transferred_via_composition():
    src = _quad_grid(3, 2.0)
    tgt = _quad_grid(4, 1.5)
    src.point_data["T"] = np.array([1.0 + i + j for j in range(4) for i in range(4)])

    out = conservative_interpolate(src, tgt)
    assert "T" in out.point_data
    assert np.all(np.isfinite(out.point_data["T"]))


@needs_core
def test_present_in_both_locations_transfers_both():
    src = _one_triangle([0, 0, 0], [1, 0, 0], [0, 1, 0])
    src.cell_data["f"] = [np.array([2.0])]
    src.point_data["f"] = np.array([1.0, 2.0, 3.0])
    tgt = _one_triangle([0, 0, 0], [1, 0, 0], [0, 1, 0])

    out = conservative_interpolate(src, tgt)
    assert "f" in out.cell_data
    assert "f" in out.point_data


@needs_core
def test_default_value_fills_uncovered_cells():
    src = _one_triangle([0, 0, 0], [1, 0, 0], [0, 1, 0])
    src.cell_data["f"] = [np.array([9.0])]
    tgt = _one_triangle([10, 10, 0], [11, 10, 0], [10, 11, 0])

    out = conservative_interpolate(src, tgt, default_value=-1.0)
    assert out.cell_data["f"][0][0] == pytest.approx(-1.0)


@needs_core
def test_conflict_policies():
    src = _one_triangle([0, 0, 0], [1, 0, 0], [0, 1, 0])
    src.cell_data["f"] = [np.array([2.0])]
    tgt = _one_triangle([0, 0, 0], [1, 0, 0], [0, 1, 0])
    tgt.cell_data["f"] = [np.array([99.0])]

    with pytest.raises(ValueError):
        conservative_interpolate(src, tgt)

    out_ow = conservative_interpolate(src, tgt, on_conflict="overwrite")
    assert out_ow.cell_data["f"][0][0] == pytest.approx(2.0)

    out_sfx = conservative_interpolate(src, tgt, on_conflict="suffix")
    assert "f_interp" in out_sfx.cell_data
    assert out_sfx.cell_data["f"][0][0] == pytest.approx(99.0)


@needs_core
def test_unknown_array_raises():
    src = _one_triangle([0, 0, 0], [1, 0, 0], [0, 1, 0])
    tgt = _one_triangle([0, 0, 0], [1, 0, 0], [0, 1, 0])
    with pytest.raises(ValueError):
        conservative_interpolate(src, tgt, arrays=["nope"])


@needs_core
def test_mismatched_topological_dimension_raises():
    src = _one_triangle([0, 0, 0], [1, 0, 0], [0, 1, 0])
    src.cell_data["f"] = [np.array([1.0])]
    tgt = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    with pytest.raises(ValueError):
        conservative_interpolate(src, tgt)


@needs_core
def test_target_geometry_and_own_data_untouched():
    src = _one_triangle([0, 0, 0], [1, 0, 0], [0, 1, 0])
    src.cell_data["f"] = [np.array([2.0])]
    tgt = _one_triangle([0, 0, 0], [1, 0, 0], [0, 1, 0])
    tgt.cell_data["g"] = [np.array([42.0])]

    out = conservative_interpolate(src, tgt)
    assert len(out.points) == len(tgt.points)
    assert out.cell_data["g"][0][0] == pytest.approx(42.0)


@needs_core
def test_ragged_polyhedron_source_and_target_are_accepted():
    src = _one_polyhedron_tet([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1])
    src.cell_data["f"] = [np.array([3.5])]
    tgt = _one_polyhedron_tet([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1])

    out = conservative_interpolate(src, tgt)
    assert out.cell_data["f"][0][0] == pytest.approx(3.5, abs=1e-9)


def test_conservative_interpolate_is_exported():
    assert "conservative_interpolate" in meshioplusplus.__all__
    assert meshioplusplus.conservative_interpolate is conservative_interpolate


def _tet_volumes(mesh):
    p = np.asarray(mesh.points, dtype=float)
    c = p[np.concatenate([np.asarray(b.data) for b in mesh.cells])]
    cross = np.cross(c[:, 1] - c[:, 0], c[:, 2] - c[:, 0])
    return np.einsum("ij,ij->i", cross, c[:, 3] - c[:, 0]) / 6.0


@pytest.mark.parametrize(
    "height, depth, layers",
    [(1.05928007, 1.0, 2), (1.0, 1.0, 2), (0.3, 1.0, 2), (1.0, 1.5, 3)],
)
def test_aligned_lattices_conserve_mass_exactly(height, depth, layers):
    # Found by tests/python/test_properties.py: two simplexified lattices of
    # one box whose tetrahedra share cutting planes. A source corner lying ON
    # a target face plane was left out of the clip's capping polygon, so one
    # overlap came out ~6% small and the transferred mass 0.7% short. Then a
    # corner on a target plane only up to rounding (a layer at z = 1/3) was
    # clipped as outside but capped as on-plane, losing a sliver.
    src = meshioplusplus.convert_cells(
        meshioplusplus.grid((1, 1, 1), spacing=(1.0, depth, height)), "simplexify"
    )
    if layers == 2:
        tgt_grid = meshioplusplus.grid((1, 2, 1), spacing=(1.0, depth / 2, height))
    else:
        tgt_grid = meshioplusplus.grid(
            (1, 1, layers), spacing=(1.0, depth, height / layers)
        )
    tgt = meshioplusplus.convert_cells(tgt_grid, "simplexify")
    src = meshioplusplus.Mesh(src.points, [(b.type, b.data) for b in src.cells])
    tgt = meshioplusplus.Mesh(tgt.points, [(b.type, b.data) for b in tgt.cells])
    values = np.array(
        [6.4059207, 2.77088847, 0.50563789, 0.26362359, 8.15137537, 9.13628022]
    )
    src.cell_data["m"] = [values]
    out = conservative_interpolate(src, tgt, arrays=["m"])
    mass_in = float(values @ _tet_volumes(src))
    mass_out = float(np.concatenate(out.cell_data["m"]) @ _tet_volumes(out))
    assert mass_out == pytest.approx(mass_in, rel=1e-12)
