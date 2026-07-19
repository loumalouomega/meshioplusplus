"""Tests for the split operation (by type / component / region)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import split


def _mixed_tri_quad():
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0]], float
    )
    return meshioplusplus.Mesh(
        pts,
        [
            ("triangle", np.array([[0, 1, 3]])),
            ("quad", np.array([[1, 2, 5, 4]])),
        ],
    )


def test_split_by_type():
    mesh = _mixed_tri_quad()
    pieces = split(mesh, by="type")
    assert set(pieces) == {"triangle", "quad"}
    # each piece has exactly its own block
    assert [cb.type for cb in pieces["triangle"].cells] == ["triangle"]
    assert [cb.type for cb in pieces["quad"].cells] == ["quad"]
    # cells sum to the original
    total = sum(sum(len(cb.data) for cb in p.cells) for p in pieces.values())
    assert total == 2
    # each piece pruned to its own points
    for p in pieces.values():
        used = set()
        for cb in p.cells:
            used |= set(cb.data.reshape(-1).tolist())
        assert used == set(range(len(p.points)))


def test_split_by_component_two_blobs():
    # two disjoint triangles
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [5, 0, 0], [6, 0, 0], [5, 1, 0]], float
    )
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2], [3, 4, 5]]))])
    pieces = split(mesh, by="component")
    assert len(pieces) == 2
    for p in pieces.values():
        assert len(p.points) == 3
        assert sum(len(cb.data) for cb in p.cells) == 1


def test_split_by_region_tag():
    mesh = _mixed_tri_quad()
    mesh.cell_data["region"] = [np.array([7], np.int64), np.array([9], np.int64)]
    pieces = split(mesh, by="region", tag="region")
    assert set(pieces) == {"7", "9"}


def test_split_by_region_cell_sets():
    mesh = _mixed_tri_quad()
    # a cell_set selecting the triangle block's only cell
    mesh.cell_sets = {"grp": [np.array([0]), np.array([], dtype=np.int64)]}
    pieces = split(mesh, by="region")
    assert "grp" in pieces
    assert sum(len(cb.data) for cb in pieces["grp"].cells) == 1


def test_cpp_matches_python():
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._split import _split_py

    mesh = _mixed_tri_quad()
    got = {p["key"]: p["mesh"] for p in core.split(mesh, "type", "")}
    ref = {k: mesh_ for k, mesh_, _, _ in _split_py(mesh, "type", "")}
    assert set(got) == set(ref)
    for k in got:
        assert len(got[k].points) == len(ref[k].points)
        assert sum(len(cb.data) for cb in got[k].cells) == sum(
            len(cb.data) for cb in ref[k].cells
        )


def test_roundtrip_write_read(tmp_path):
    mesh = _mixed_tri_quad()
    pieces = split(mesh, by="type")
    for key, piece in pieces.items():
        p = tmp_path / f"s_{key}.vtu"
        meshioplusplus.write(p, piece)
        back = meshioplusplus.read(p)
        assert len(back.points) == len(piece.points)
