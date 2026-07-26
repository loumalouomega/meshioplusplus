import copy
import pathlib

import numpy as np
import pytest

import meshioplusplus


def _quad_grid(nx=6, ny=4):
    pts = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            pts.append([float(i), float(j), 0.0])
    cells = []
    for j in range(ny):
        for i in range(nx):
            v = j * (nx + 1) + i
            cells.append([v, v + 1, v + nx + 2, v + nx + 1])
    return meshioplusplus.Mesh(
        np.array(pts, dtype=float), [("quad", np.array(cells, dtype=np.int64))]
    )


def _mixed_tri_quad():
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ]
    )
    return meshioplusplus.Mesh(
        points,
        [
            ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
            ("quad", np.array([[1, 4, 5, 2]])),
        ],
    )


def _total_cells(mesh):
    return sum(len(cb.data) for cb in mesh.cells)


def _flat(labels):
    return np.concatenate([np.asarray(a).reshape(-1) for a in labels])


def _has_kahip():
    try:
        from meshioplusplus import _core

        if getattr(_core, "__has_kahip__", False):
            return True
    except Exception:
        pass
    try:
        import kahip  # noqa: F401

        return True
    except ImportError:
        return False


def test_labels_block_aligned_and_in_range():
    mesh = _quad_grid()
    labels = meshioplusplus.partition_labels(mesh, 4, method="sfc")
    assert len(labels) == len(mesh.cells)
    flat = _flat(labels)
    assert len(flat) == _total_cells(mesh)
    assert flat.min() >= 0
    assert flat.max() < 4


def test_unweighted_sizes_differ_by_at_most_one():
    mesh = _quad_grid(6, 4)  # 24 cells
    for nparts in (2, 3, 5, 7):
        flat = _flat(meshioplusplus.partition_labels(mesh, nparts, method="sfc"))
        sizes = np.bincount(flat.astype(int), minlength=nparts)
        assert sizes.max() - sizes.min() <= 1


def test_pieces_are_partition_of_unity():
    mesh = _quad_grid()
    pieces = meshioplusplus.partition(mesh, 3, method="sfc")
    assert len(pieces) == 3
    assert sum(_total_cells(p) for p in pieces) == _total_cells(mesh)
    # Blocks are kept 1:1 with the input (unlike split).
    for piece in pieces:
        assert len(piece.cells) == len(mesh.cells)
        assert [cb.type for cb in piece.cells] == [cb.type for cb in mesh.cells]


def test_pieces_agree_with_labels():
    mesh = _mixed_tri_quad()
    labels = meshioplusplus.partition_labels(mesh, 2, method="sfc")
    pieces = meshioplusplus.partition(mesh, 2, method="sfc", record_ids=True)
    for part, piece in enumerate(pieces):
        for b, ids in enumerate(piece.cell_data["partition:original_cell_id"]):
            for orig in np.asarray(ids, dtype=np.int64):
                assert int(np.asarray(labels[b])[orig]) == part


def test_cpp_matches_python():
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._partition import _labels_py

    mesh = _quad_grid(7, 5)
    got = np.concatenate(
        [
            np.asarray(a)
            for a in core.partition_labels(mesh, 4, "sfc", 0.03, "eco", 0, "")
        ]
    )
    want = _labels_py(mesh, 4, "sfc", 0.03, "eco", 0, None)
    assert np.array_equal(got, want)


def test_weighted_cut_balances_weight():
    mesh = _quad_grid(8, 1)  # 8 cells in a row
    w = np.ones(8)
    w[0] = 100.0
    mesh.cell_data["w"] = [w]
    flat = _flat(meshioplusplus.partition_labels(mesh, 2, method="sfc", weights="w"))
    wsum = np.zeros(2)
    for c in range(8):
        wsum[flat[c]] += w[c]
    assert wsum[0] <= w.sum() / 2 + w.max()
    assert wsum[1] <= w.sum() / 2 + w.max()


def test_weights_errors():
    mesh = _quad_grid(2, 2)
    with pytest.raises(ValueError):
        meshioplusplus.partition_labels(mesh, 2, weights="nope")
    mesh.cell_data["v"] = [np.ones((4, 2))]
    with pytest.raises(ValueError):
        meshioplusplus.partition_labels(mesh, 2, weights="v")
    mesh.cell_data["neg"] = [np.array([1.0, -1.0, 1.0, 1.0])]
    with pytest.raises(ValueError):
        meshioplusplus.partition_labels(mesh, 2, weights="neg")


def test_invalid_nparts_raises():
    mesh = _quad_grid(2, 2)
    with pytest.raises(ValueError):
        meshioplusplus.partition(mesh, 0)
    with pytest.raises(ValueError):
        meshioplusplus.partition_labels(mesh, -1)


def test_negative_ghost_layers_raises():
    mesh = _quad_grid(2, 2)
    with pytest.raises(ValueError, match="ghost_layers"):
        meshioplusplus.partition(mesh, 2, ghost_layers=-1)


def test_ghost_layers_with_python_kahip_route(monkeypatch):
    # The CI "kahip PyPI wheel + KaHIP-less _core" leg: with the wheel importable
    # and no KaHIP in the core, method="auto"/"kahip" normally bypass _core for
    # the pure-Python quality path — which builds disjoint pieces only. Ghosting
    # must then route "auto" back to the C++ core (auto is a preference order,
    # and the core is the one ghost-capable path), while an EXPLICIT "kahip"
    # fails by name rather than silently downgrading to SFC.
    import sys
    import types

    from meshioplusplus import _core, _partition

    if getattr(_core, "__has_kahip__", False):
        pytest.skip("core has KaHIP: the python-kahip bypass never engages")

    monkeypatch.setitem(sys.modules, "kahip", types.ModuleType("kahip"))
    mesh = _quad_grid(4, 4)
    assert _partition._use_python_kahip("auto"), "stub should engage the bypass"

    pieces = meshioplusplus.partition(mesh, 2, ghost_layers=1)  # auto
    assert len(pieces) == 2
    assert all("partition:ghost" in p.cell_data for p in pieces)

    with pytest.raises(NotImplementedError, match="kahip"):
        meshioplusplus.partition(mesh, 2, method="kahip", ghost_layers=1)


def test_ghost_layers_grow_pieces_and_tag_the_halo():
    mesh = _quad_grid(4, 4)
    plain = meshioplusplus.partition(mesh, 2)
    ghosted = meshioplusplus.partition(mesh, 2, ghost_layers=1)
    assert len(ghosted) == len(plain)

    for owned_piece, halo_piece in zip(plain, ghosted):
        owned = sum(len(cb.data) for cb in owned_piece.cells)
        with_halo = sum(len(cb.data) for cb in halo_piece.cells)
        assert with_halo > owned, "a halo must add cells on a connected grid"

        tag = halo_piece.cell_data["partition:ghost"]
        assert len(tag) == len(halo_piece.cells)
        flat = np.concatenate([np.asarray(a).ravel() for a in tag])
        assert flat.min() >= 0 and flat.max() <= 1  # one layer requested
        assert int((flat == 0).sum()) == owned  # depth 0 is exactly the owned set
        assert len(flat) == with_halo


def test_ghost_layers_are_monotone():
    mesh = _quad_grid(5, 5)
    prev = 0
    for layers in range(4):
        piece = meshioplusplus.partition(mesh, 2, ghost_layers=layers)[0]
        n = sum(len(cb.data) for cb in piece.cells)
        assert n >= prev
        assert n <= 25
        prev = n


def test_labels_take_no_ghost_layers():
    # A flat per-cell label array is the ownership map; a cell can be a ghost of
    # several parts at once, so there is nowhere to put that. The Python API
    # simply does not expose the parameter (the C++ one validates it instead).
    mesh = _quad_grid(2, 2)
    with pytest.raises(TypeError):
        meshioplusplus.partition_labels(mesh, 2, ghost_layers=1)


def test_unknown_method_and_mode_raise():
    mesh = _quad_grid(2, 2)
    with pytest.raises(ValueError, match="method"):
        meshioplusplus.partition_labels(mesh, 2, method="metis")
    with pytest.raises(ValueError, match="mode"):
        meshioplusplus.partition_labels(mesh, 2, mode="turbo")


def test_kahip_absent_fails_by_name():
    if _has_kahip():
        pytest.skip("a KaHIP backend is available")
    mesh = _quad_grid(2, 2)
    with pytest.raises(ValueError, match="(?i)kahip"):
        meshioplusplus.partition_labels(mesh, 2, method="kahip")


def test_sets_are_remapped_per_piece():
    mesh = _mixed_tri_quad()
    mesh.point_sets = {"left": np.array([0, 3], dtype=np.int64)}
    mesh.cell_sets = {
        "first": [np.array([0], dtype=np.int64), np.array([], dtype=np.int64)]
    }
    pieces = meshioplusplus.partition(mesh, 2, method="sfc")
    seen_points = 0
    seen_cells = 0
    for piece in pieces:
        assert "left" in piece.point_sets
        seen_points += len(piece.point_sets["left"])
        assert "first" in piece.cell_sets
        seen_cells += sum(len(b) for b in piece.cell_sets["first"])
    # A set point may be shared by pieces; the set cell lives in exactly one.
    assert seen_points >= 2
    assert seen_cells == 1


def test_record_ids_round_trip():
    mesh = _quad_grid(3, 3)
    pieces = meshioplusplus.partition(mesh, 2, method="sfc", record_ids=True)
    for piece in pieces:
        assert "partition:original_point_id" in piece.point_data
        ids = np.asarray(
            piece.point_data["partition:original_point_id"], dtype=np.int64
        )
        # Original ids point back at identical coordinates.
        assert np.allclose(np.asarray(mesh.points)[ids], np.asarray(piece.points))


def test_point_and_cell_data_carried():
    mesh = _quad_grid(3, 2)
    mesh.point_data["t"] = np.arange(len(mesh.points), dtype=float)
    mesh.cell_data["c"] = [np.arange(_total_cells(mesh), dtype=float)]
    pieces = meshioplusplus.partition(mesh, 2, method="sfc")
    for piece in pieces:
        assert "t" in piece.point_data
        assert "c" in piece.cell_data
        assert len(piece.cell_data["c"][0]) == _total_cells(piece)


def test_numpy_fallback_matches_core_pieces():
    pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._partition import _build_piece, _labels_py

    mesh = _quad_grid(4, 3)
    flat = _labels_py(mesh, 3, "sfc", 0.03, "eco", 0, None)
    pieces = meshioplusplus.partition(mesh, 3, method="sfc")
    for part, piece in enumerate(pieces):
        out, _pm, _cm = _build_piece(mesh, flat, part, False)
        assert np.allclose(np.asarray(out.points), np.asarray(piece.points))
        for a, b in zip(out.cells, piece.cells):
            assert a.type == b.type
            assert np.array_equal(np.asarray(a.data), np.asarray(b.data))


def test_roundtrip_write_read(tmp_path: pathlib.Path):
    mesh = _quad_grid(4, 2)
    pieces = meshioplusplus.partition(mesh, 2, method="sfc")
    total = 0
    for part, piece in enumerate(pieces):
        path = tmp_path / f"piece_{part}.vtu"
        meshioplusplus.write(path, piece)
        back = meshioplusplus.read(path)
        total += _total_cells(back)
    assert total == _total_cells(mesh)


def test_deterministic_across_runs():
    mesh = _quad_grid(5, 5)
    a = _flat(meshioplusplus.partition_labels(mesh, 4, method="sfc"))
    b = _flat(meshioplusplus.partition_labels(mesh, 4, method="sfc"))
    assert np.array_equal(a, b)


def test_input_mesh_is_unchanged():
    mesh = _quad_grid(3, 3)
    ref = copy.deepcopy(mesh)
    meshioplusplus.partition(mesh, 3, method="sfc", record_ids=True)
    assert np.array_equal(np.asarray(mesh.points), np.asarray(ref.points))
    assert not mesh.point_data and not mesh.cell_data


# --------------------------------------------------------------------------- #
# KaHIP (balance/coverage only -- assignment varies by KaHIP version)          #
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not _has_kahip(), reason="no KaHIP backend available")
def test_kahip_coverage_and_balance():
    mesh = _quad_grid(8, 8)  # 64 cells
    flat = _flat(
        meshioplusplus.partition_labels(mesh, 4, method="kahip", imbalance=0.03)
    )
    assert len(flat) == 64
    sizes = np.bincount(flat.astype(int), minlength=4)
    assert sizes.min() > 0
    assert sizes.max() <= 64 / 4 * 1.03 + 1


@pytest.mark.skipif(not _has_kahip(), reason="no KaHIP backend available")
def test_kahip_pieces_partition_of_unity():
    mesh = _quad_grid(6, 6)
    pieces = meshioplusplus.partition(mesh, 3, method="kahip")
    assert len(pieces) == 3
    assert sum(_total_cells(p) for p in pieces) == _total_cells(mesh)


@pytest.mark.skipif(not _has_kahip(), reason="no KaHIP backend available")
def test_kahip_fixed_seed_reproducible():
    mesh = _quad_grid(6, 6)
    a = _flat(meshioplusplus.partition_labels(mesh, 4, method="kahip", seed=42))
    b = _flat(meshioplusplus.partition_labels(mesh, 4, method="kahip", seed=42))
    assert np.array_equal(a, b)


@pytest.mark.skipif(not _has_kahip(), reason="no KaHIP backend available")
def test_kahip_invalid_imbalance_raises():
    mesh = _quad_grid(2, 2)
    with pytest.raises(ValueError, match="imbalance"):
        meshioplusplus.partition_labels(mesh, 2, method="kahip", imbalance=1.5)
