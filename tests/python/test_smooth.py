"""Mesh smoothing (``meshioplusplus.smooth``)."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._smooth import _smooth_py

from . import helpers


def _quad_grid(n=9):
    """An n x n regular grid of unit quads in 2D."""
    xs, ys = np.meshgrid(
        np.arange(n, dtype=float), np.arange(n, dtype=float), indexing="ij"
    )
    points = np.column_stack([xs.ravel(), ys.ravel()])

    def pid(i, j):
        return i * n + j

    cells = np.array(
        [
            [pid(i, j), pid(i + 1, j), pid(i + 1, j + 1), pid(i, j + 1)]
            for i in range(n - 1)
            for j in range(n - 1)
        ]
    )
    return mp.Mesh(points, [("quad", cells)])


def _hex_block(n=4):
    """An n x n x n block of unit hexahedra."""
    g = np.meshgrid(*[np.arange(n + 1, dtype=float)] * 3, indexing="ij")
    points = np.column_stack([a.ravel() for a in g])

    def pid(i, j, k):
        return (i * (n + 1) + j) * (n + 1) + k

    cells = np.array(
        [
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
            for i in range(n)
            for j in range(n)
            for k in range(n)
        ]
    )
    return mp.Mesh(points, [("hexahedron", cells)])


def _boundary_ids(n=9):
    return [i for i in range(n * n) if i // n in (0, n - 1) or i % n in (0, n - 1)]


def _interior_ids(n=9):
    return np.array([i for i in range(n * n) if i not in set(_boundary_ids(n))])


# --- the core promise -------------------------------------------------------


def test_displaced_interior_node_returns_to_its_centroid():
    """A single node pushed off-centre in a regular grid is pulled back."""
    mesh = _quad_grid()
    node = 4 * 9 + 4  # the middle of the grid
    target = mesh.points[node].copy()
    mesh.points[node] += [0.45, -0.35]
    before = np.linalg.norm(mesh.points[node] - target)

    out = mp.smooth(mesh, method="laplacian", iterations=8)
    after = np.linalg.norm(out.points[node] - target)

    assert before > 0.5
    assert after < before / 10


def test_taubin_does_not_shrink_where_laplacian_does():
    """The whole reason Taubin is the default."""
    rng = np.random.default_rng(0)
    mesh = _quad_grid()
    interior = _interior_ids()
    mesh.points[interior] += rng.normal(0, 0.18, (len(interior), 2))
    extent = np.ptp(mesh.points, axis=0)

    lap = mp.smooth(
        mesh,
        method="laplacian",
        iterations=40,
        fix_boundary=False,
        preserve_features=False,
        guard_inversion=False,
    )
    tau = mp.smooth(
        mesh,
        method="taubin",
        iterations=40,
        fix_boundary=False,
        preserve_features=False,
        guard_inversion=False,
    )

    lap_keep = (np.ptp(lap.points, axis=0) / extent).mean()
    tau_keep = (np.ptp(tau.points, axis=0) / extent).mean()

    assert lap_keep < 0.6, "laplacian is expected to shrink the mesh substantially"
    assert tau_keep > 0.95, "taubin is expected to preserve the bounding box"


def test_fix_boundary_pins_boundary_nodes_exactly():
    rng = np.random.default_rng(1)
    mesh = _quad_grid()
    interior = _interior_ids()
    mesh.points[interior] += rng.normal(0, 0.15, (len(interior), 2))
    boundary = _boundary_ids()

    out = mp.smooth(mesh, iterations=12)

    # Bit-identical, not merely close: pinned means pinned.
    assert np.array_equal(out.points[boundary], mesh.points[boundary])
    assert not np.array_equal(out.points[interior], mesh.points[interior])


def test_geometry_only_connectivity_and_data_values_unchanged():
    rng = np.random.default_rng(2)
    mesh = _quad_grid()
    mesh.points[_interior_ids()] += rng.normal(0, 0.1, (len(_interior_ids()), 2))
    mesh.point_data["h"] = mesh.points[:, 0].copy()
    mesh.cell_data["tag"] = [np.arange(len(mesh.cells[0].data))]
    mesh.field_data["meta"] = np.array([1, 2])

    out = mp.smooth(mesh, iterations=5)

    assert len(out.points) == len(mesh.points)
    assert len(out.cells) == len(mesh.cells)
    assert np.array_equal(out.cells[0].data, mesh.cells[0].data)
    assert np.array_equal(out.point_data["h"], mesh.point_data["h"])
    assert np.array_equal(out.cell_data["tag"][0], mesh.cell_data["tag"][0])
    assert np.array_equal(out.field_data["meta"], mesh.field_data["meta"])
    assert not np.array_equal(out.points, mesh.points)


def test_points_dtype_is_preserved():
    mesh = _quad_grid()
    mesh.points = mesh.points.astype(np.float32)
    out = mp.smooth(mesh, iterations=3)
    assert out.points.dtype == np.float32


# --- the inversion guard ----------------------------------------------------


def test_guard_prevents_new_inversions():
    """Smoothing a badly tangled block must not leave inverted cells behind."""
    rng = np.random.default_rng(3)
    mesh = _hex_block()
    n = 4
    interior = [
        i for i in range(len(mesh.points)) if all(0 < c < n for c in mesh.points[i])
    ]
    mesh.points[interior] += rng.normal(0, 0.75, (len(interior), 3))

    out, report = mp.smooth(mesh, iterations=15, return_report=True)

    assert mp.compute_stats(out)["num_inverted"] == 0
    assert report["num_skipped_inversion"] > 0


def test_guard_does_not_lock_in_pre_existing_inversions():
    """The guard is 'do no harm', not 'preserve the sign'.

    An earlier version rejected any sign change, which pinned every cell that
    arrived inverted and so prevented smoothing from repairing a tangle -- the
    opposite of the point.
    """
    rng = np.random.default_rng(3)
    mesh = _hex_block()
    n = 4
    interior = [
        i for i in range(len(mesh.points)) if all(0 < c < n for c in mesh.points[i])
    ]
    mesh.points[interior] += rng.normal(0, 0.75, (len(interior), 3))
    assert mp.compute_stats(mesh)["num_inverted"] > 0, "fixture must start tangled"

    out = mp.smooth(mesh, iterations=15)
    assert mp.compute_stats(out)["num_inverted"] == 0


def test_quality_improves():
    rng = np.random.default_rng(3)
    mesh = _hex_block()
    n = 4
    interior = [
        i for i in range(len(mesh.points)) if all(0 < c < n for c in mesh.points[i])
    ]
    mesh.points[interior] += rng.normal(0, 0.35, (len(interior), 3))

    out = mp.smooth(mesh, iterations=20)

    before = mp.compute_quality(mesh)["metrics"]
    after = mp.compute_quality(out)["metrics"]
    assert (
        after["quality:scaled_jacobian"]["min"]
        > before["quality:scaled_jacobian"]["min"]
    )
    assert after["quality:skewness"]["mean"] < before["quality:skewness"]["mean"]


def test_smoothing_never_increases_the_inverted_count():
    """Guards against the two face-fan copies (smooth's and quality's) drifting.

    smooth_facefan_volume is a deliberate flat-buffer duplicate of quality.cpp's
    quality_facefan_volume; if the two ever disagree in sign, this fails.
    """
    rng = np.random.default_rng(7)
    for sigma in (0.1, 0.3, 0.5):
        mesh = _hex_block()
        n = 4
        interior = [
            i for i in range(len(mesh.points)) if all(0 < c < n for c in mesh.points[i])
        ]
        mesh.points[interior] += rng.normal(0, sigma, (len(interior), 3))
        before = mp.compute_stats(mesh)["num_inverted"]
        after = mp.compute_stats(mp.smooth(mesh, iterations=10))["num_inverted"]
        assert after <= before, f"sigma={sigma}: {before} -> {after}"


# --- pinning ----------------------------------------------------------------


def test_frozen_node_ids_are_pinned():
    rng = np.random.default_rng(4)
    mesh = _quad_grid()
    interior = _interior_ids()
    mesh.points[interior] += rng.normal(0, 0.15, (len(interior), 2))
    frozen = [int(interior[0]), int(interior[5])]

    out = mp.smooth(mesh, iterations=10, frozen=frozen)
    assert np.array_equal(out.points[frozen], mesh.points[frozen])


def test_frozen_accepts_a_point_set_name():
    mesh = _quad_grid()
    interior = _interior_ids()
    mesh.point_sets = {"pinned": np.asarray(interior[:3], dtype=np.int64)}
    out = mp.smooth(mesh, iterations=5, frozen="pinned")
    assert np.array_equal(out.points[interior[:3]], mesh.points[interior[:3]])
    # Sets are carried through untouched (smooth renumbers nothing).
    assert "pinned" in out.point_sets


def test_unknown_point_set_name_raises():
    with pytest.raises(ValueError, match="no point_set named"):
        mp.smooth(_quad_grid(), frozen="nope")


def test_higher_order_cells_are_pinned_not_distorted():
    """A block with no known edge topology holds still rather than being guessed."""
    points = np.array(
        [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    )
    mesh = mp.Mesh(points, [("triangle6", np.array([[0, 1, 2, 3, 4, 5]]))])
    out = mp.smooth(mesh, iterations=10, guard_inversion=False)
    assert np.array_equal(out.points, mesh.points)


# --- validation -------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"method": "bogus"}, "unknown method"),
        ({"lambda_": 1.5}, "lambda must lie"),
        ({"method": "taubin", "lambda_": 0.4, "mu": -0.2}, "mu < -lambda"),
        ({"frozen": [10_000]}, "out of range"),
    ],
)
def test_invalid_options_raise(kwargs, match):
    with pytest.raises(ValueError, match=match):
        mp.smooth(_quad_grid(), **kwargs)


def test_zero_iterations_is_a_no_op():
    mesh = _quad_grid()
    out = mp.smooth(mesh, iterations=0)
    assert np.array_equal(out.points, mesh.points)


# --- cross-implementation + determinism -------------------------------------


@pytest.mark.parametrize("method, lam", [("laplacian", 0.5), ("taubin", 0.33)])
def test_cpp_matches_python(method, lam):
    """The C++ core and the numpy fallback must agree.

    The guard is off: it is a discrete branch on the sign of a cell measure, so
    near a degenerate cell the two implementations could legitimately land on
    opposite sides and then diverge macroscopically rather than by round-off.
    The continuous part of the algorithm is what this pins.
    """
    from meshioplusplus import _core

    rng = np.random.default_rng(5)
    mesh = _quad_grid()
    interior = _interior_ids()
    mesh.points[interior] += rng.normal(0, 0.15, (len(interior), 2))

    got = _core.smooth(mesh, method, 6, lam, -0.34, True, True, 30.0, False)["mesh"]
    ref, _ = _smooth_py(mesh, method, 6, lam, -0.34, True, True, 30.0, False)

    np.testing.assert_allclose(
        np.asarray(got.points), np.asarray(ref.points), rtol=1e-9, atol=1e-12
    )


def test_repeated_runs_are_byte_identical():
    rng = np.random.default_rng(6)
    mesh = _hex_block()
    mesh.points += rng.normal(0, 0.05, mesh.points.shape)
    first = mp.smooth(mesh, iterations=7).points
    for _ in range(4):
        assert np.array_equal(mp.smooth(mesh, iterations=7).points, first)


def test_roundtrip_write_read(tmp_path):
    mesh = _hex_block()
    rng = np.random.default_rng(8)
    mesh.points += rng.normal(0, 0.05, mesh.points.shape)
    out = mp.smooth(mesh, iterations=4)

    path = tmp_path / "smoothed.vtu"
    mp.write(path, out)
    back = mp.read(path)

    np.testing.assert_allclose(back.points, out.points, rtol=1e-12, atol=1e-14)
    assert np.array_equal(back.cells[0].data, out.cells[0].data)


def test_helpers_meshes_survive_smoothing():
    """Smoothing must not corrupt the shared fixtures or change their topology."""
    for mesh in (helpers.tri_mesh, helpers.quad_mesh, helpers.tet_mesh):
        import copy

        m = copy.deepcopy(mesh)
        out = mp.smooth(m, iterations=3)
        assert len(out.points) == len(m.points)
        assert [b.type for b in out.cells] == [b.type for b in m.cells]


# --- ODT smoothing -----------------------------------------------------------

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core_odt = pytest.mark.skipif(
    _core is None,
    reason="method='odt' is C++-core only and has no pure-Python fallback",
)


def _tet_block(n=4):
    """A jittered n x n x n block of tetrahedra (a hex block simplexified)."""
    hexes = _hex_block(n)
    return mp.convert_cells(hexes, mode="simplexify")


@needs_core_odt
def test_odt_moves_a_single_tets_vertex_toward_its_circumcenter():
    pts = np.array([[0, 0, 0], [4, 0, 0], [0, 4, 0], [1, 1, 5]], dtype=float)
    conn = np.array([[0, 1, 2, 3]], dtype=np.int64)
    mesh = mp.Mesh(pts, [("tetra", conn)])

    out, report = mp.smooth(
        mesh,
        method="odt",
        iterations=1,
        fix_boundary=False,
        preserve_features=False,
        return_report=True,
    )
    assert out.points.shape == (4, 3)
    assert out.cells[0].data.shape == (1, 4)
    assert report["num_nodes_moved"] > 0
    assert report["max_displacement"] > 0.0


@needs_core_odt
def test_odt_rejects_non_tet_blocks_by_name():
    with pytest.raises(Exception, match="tet"):
        mp.smooth(_hex_block(2), method="odt", iterations=1)


@needs_core_odt
def test_odt_improves_dihedral_angles_more_than_taubin():
    import copy

    from meshioplusplus import compute_quality

    rng = np.random.default_rng(11)
    base = _tet_block(4)
    lo, hi = base.points.min(axis=0), base.points.max(axis=0)
    is_interior = np.all((base.points > lo + 0.1) & (base.points < hi - 0.1), axis=1)

    jittered = copy.deepcopy(base)
    jittered.points[is_interior] += rng.normal(0, 0.15, (int(is_interior.sum()), 3))

    taubin = mp.smooth(jittered, method="taubin", iterations=10)
    odt = mp.smooth(jittered, method="odt", iterations=10)

    q_taubin = compute_quality(taubin)
    q_odt = compute_quality(odt)
    min_dihedral_taubin = q_taubin["metrics"]["quality:min_dihedral"]["min"]
    min_dihedral_odt = q_odt["metrics"]["quality:min_dihedral"]["min"]
    assert min_dihedral_odt >= min_dihedral_taubin - 1e-9


@needs_core_odt
def test_odt_is_deterministic_across_repeated_runs():
    import copy

    rng = np.random.default_rng(12)
    base = _tet_block(3)
    jittered = copy.deepcopy(base)
    jittered.points += rng.normal(0, 0.05, jittered.points.shape)

    first = mp.smooth(jittered, method="odt", iterations=5).points
    for _ in range(3):
        again = mp.smooth(jittered, method="odt", iterations=5).points
        np.testing.assert_array_equal(again, first)
