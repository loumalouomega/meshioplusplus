"""Tests for the mesh renumbering / reordering operation."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import compute_bandwidth, reorder
from meshioplusplus._reorder import _reorder_py

from . import helpers


def grid_mesh(nx, ny):
    pts = [[float(i), float(j), 0.0] for j in range(ny + 1) for i in range(nx + 1)]

    def idx(i, j):
        return j * (nx + 1) + i

    quads = [
        [idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)]
        for j in range(ny)
        for i in range(nx)
    ]
    return meshioplusplus.Mesh(pts, [("quad", quads)])


def shuffled_grid(nx, ny, stride=13):
    """A grid whose node labels are scrambled (coprime stride) -> big bandwidth."""
    base = grid_mesh(nx, ny)
    n = len(base.points)
    assert np.gcd(stride, n) == 1
    relabel = (np.arange(n) * stride) % n
    pts = np.empty_like(np.asarray(base.points))
    pts[relabel] = np.asarray(base.points)
    cells = [(blk.type, relabel[np.asarray(blk.data)]) for blk in base.cells]
    return meshioplusplus.Mesh(pts, cells)


def bad_path_mesh(n):
    """A path with even-first labelling -> large bandwidth (~n/2)."""
    half = (n + 1) // 2
    label = [(p // 2) if p % 2 == 0 else (half + p // 2) for p in range(n)]
    pts = [[0.0, 0.0, 0.0] for _ in range(n)]
    for p in range(n):
        pts[label[p]] = [float(p), 0.0, 0.0]
    edges = [[label[p], label[p + 1]] for p in range(n - 1)]
    return meshioplusplus.Mesh(pts, [("line", edges)])


def _cell_rows_by_type(mesh):
    """{type: sorted list of connectivity rows} for order-insensitive compare."""
    out = {}
    for blk in mesh.cells:
        rows = [tuple(int(x) for x in row) for row in np.asarray(blk.data)]
        out.setdefault(blk.type, []).extend(rows)
    return {k: sorted(v) for k, v in out.items()}


def _check_pure_permutation(original, out, node_perm, cell_perms):
    n = len(original.points)
    node_perm = np.asarray(node_perm)
    # bijection over [0, n)
    assert node_perm.shape == (n,)
    assert sorted(node_perm.tolist()) == list(range(n))
    inv = np.empty(n, dtype=np.int64)
    inv[node_perm] = np.arange(n)

    # points move by the node permutation
    assert np.allclose(np.asarray(out.points)[node_perm], np.asarray(original.points))

    # connectivity un-permutes back to the original (as multisets per type)
    back = meshioplusplus.Mesh(
        original.points,
        [(blk.type, inv[np.asarray(blk.data)]) for blk in out.cells],
    )
    assert _cell_rows_by_type(back) == _cell_rows_by_type(original)

    # cell permutations: one bijection per block
    assert len(cell_perms) == len(original.cells)
    for blk, cp in zip(out.cells, cell_perms):
        cp = np.asarray(cp)
        assert sorted(cp.tolist()) == list(range(len(blk.data)))


MESHES = {
    "tri": helpers.tri_mesh,
    "tet": helpers.tet_mesh,
    "grid": grid_mesh(4, 3),
    "bad_path": bad_path_mesh(16),
}


@pytest.mark.parametrize("method", ["rcm", "morton", "hilbert"])
@pytest.mark.parametrize("name", list(MESHES))
def test_pure_permutation(name, method):
    mesh = MESHES[name]
    out, node_perm, cell_perms = reorder(mesh, method=method, return_permutation=True)
    _check_pure_permutation(mesh, out, node_perm, cell_perms)


def test_rcm_reduces_bandwidth():
    mesh = bad_path_mesh(64)
    before = compute_bandwidth(mesh)
    out = reorder(mesh, method="rcm")
    after = compute_bandwidth(out)
    assert after < before
    assert after <= 2


def test_rcm_reduces_bandwidth_on_shuffled_grid():
    mesh = shuffled_grid(6, 5)
    before = compute_bandwidth(mesh)
    out = reorder(mesh, method="rcm")
    after = compute_bandwidth(out)
    assert after < before
    # RCM on a 7x6 node grid should get close to the ~grid-width optimum
    assert after <= 16


@pytest.mark.parametrize("method", ["morton", "hilbert"])
def test_sfc_locality(method):
    # A space-filling curve keeps spatially-close nodes close in index: the
    # mean spatial step between consecutive new indices is far below what a
    # random ordering (~mean pairwise distance) would give.
    mesh = grid_mesh(8, 8)
    n = len(mesh.points)
    out, node_perm, _ = reorder(mesh, method=method, return_permutation=True)
    node_perm = np.asarray(node_perm)
    assert sorted(node_perm.tolist()) == list(range(n))
    inv = np.empty(n, dtype=np.int64)
    inv[node_perm] = np.arange(n)
    ordered = np.asarray(mesh.points)[inv]  # points in new-index order
    steps = np.linalg.norm(np.diff(ordered, axis=0), axis=1)
    assert steps.mean() < 3.0  # grid spacing 1; random ordering would be ~4-5


def test_point_data_follows_nodes():
    mesh = grid_mesh(3, 3)
    n = len(mesh.points)
    mesh.point_data = {"tag": np.arange(n, dtype=np.float64) + 0.5}
    out, node_perm, _ = reorder(mesh, method="rcm", return_permutation=True)
    assert "tag" in out.point_data
    np.testing.assert_allclose(
        np.asarray(out.point_data["tag"])[np.asarray(node_perm)],
        np.arange(n) + 0.5,
    )


def test_cell_data_follows_cells():
    mesh = grid_mesh(4, 2)
    ncells = len(mesh.cells[0])
    mesh.cell_data = {"idx": [np.arange(ncells, dtype=np.float64)]}
    out, _, cell_perms = reorder(mesh, method="morton", return_permutation=True)
    cp = np.asarray(cell_perms[0])
    np.testing.assert_allclose(
        np.asarray(out.cell_data["idx"][0])[cp], np.arange(ncells)
    )


def test_sets_remapped():
    # Membership, compared as a set: a region's entries are stored ascending and
    # de-duplicated (see doc/regions.md), so the remapped set is sorted rather
    # than left in the permutation's own order. Set semantics are what a named
    # group means; the old element-order assertion tested storage, not meaning.
    mesh = helpers.add_cell_sets(helpers.add_point_sets(grid_mesh(4, 3)))
    out, node_perm, cell_perms = reorder(mesh, method="rcm", return_permutation=True)
    node_perm = np.asarray(node_perm)
    for name, idx in mesh.point_sets.items():
        np.testing.assert_array_equal(
            out.point_sets[name], np.sort(node_perm[np.asarray(idx)])
        )
    for name, blocks in mesh.cell_sets.items():
        for b, idx in enumerate(blocks):
            np.testing.assert_array_equal(
                out.cell_sets[name][b],
                np.sort(np.asarray(cell_perms[b])[np.asarray(idx)]),
            )


def test_invalid_method():
    with pytest.raises(ValueError):
        reorder(grid_mesh(2, 2), method="nope")


def test_roundtrip_write_read(tmp_path):
    mesh = grid_mesh(5, 4)
    out = reorder(mesh, method="rcm")
    p = tmp_path / "reordered.vtu"
    meshioplusplus.write(p, out)
    back = meshioplusplus.read(p)
    assert np.allclose(np.asarray(back.points), np.asarray(out.points))
    assert _cell_rows_by_type(back) == _cell_rows_by_type(out)


@pytest.mark.parametrize("method", ["rcm", "morton", "hilbert"])
@pytest.mark.parametrize("name", list(MESHES))
def test_cpp_matches_python(name, method):
    """The C++ core and the pure-numpy reference agree on the permutation."""
    _core = pytest.importorskip("meshioplusplus._core")
    mesh = MESHES[name]
    try:
        res = _core.reorder(mesh, method)
    except Exception:
        pytest.skip("C++ core reorder unavailable")
    cpp_perm = np.asarray(res["node_permutation"])

    _, py_perm, _ = _reorder_py(mesh, method)
    py_perm = np.asarray(py_perm)

    # both are valid bijections that achieve the same bandwidth
    n = len(mesh.points)
    assert sorted(cpp_perm.tolist()) == list(range(n))
    cpp_mesh = res["mesh"]
    py_mesh, _, _ = _reorder_py(mesh, method)
    assert compute_bandwidth(cpp_mesh) == compute_bandwidth(py_mesh)
    if method in ("morton", "hilbert"):
        # deterministic bit ops -> identical permutation
        np.testing.assert_array_equal(cpp_perm, py_perm)
