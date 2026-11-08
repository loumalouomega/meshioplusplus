"""Tests for ``_set_remap.py``'s piece remapping, shared by split and partition.

``remap_piece_sets`` carries ``point_sets``/``cell_sets`` onto one piece;
``keep_invalid`` and ``drop_empty_blocks`` are the two places split
(``False``/``True``) and partition (``True``/``False``) differ.
"""

import numpy as np
from numpy.testing import assert_array_equal

from meshioplusplus import Mesh
from meshioplusplus._set_remap import (
    cell_set_is_remappable,
    point_set_is_remappable,
    remap_piece_sets,
)


def _mesh():
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], dtype=float)
    cells = [
        ("triangle", np.array([[0, 1, 2], [1, 3, 2]])),
        ("line", np.array([[0, 1], [1, 3], [3, 2]])),
    ]
    mesh = Mesh(points, cells)
    mesh.point_sets = {"corner": np.array([0, 3]), "tags": np.array([-7, 2])}
    mesh.cell_sets = {
        "sel": [np.array([1]), np.array([0, 2])],
        # gmsh's bounding entities: per-block tags, not cell indices
        "gmsh:bounding_entities": [np.array([-1, 5]), np.array([9])],
    }
    return mesh


def test_point_set_predicate():
    assert point_set_is_remappable(np.array([0, 3]), 4)
    assert point_set_is_remappable(np.array([], dtype=np.int64), 4)
    assert not point_set_is_remappable(np.array([0, 4]), 4)
    assert not point_set_is_remappable(np.array([-1]), 4)
    assert not point_set_is_remappable(np.array([0.5]), 4)


def test_cell_set_predicate_is_strict_about_range():
    cells = _mesh().cells
    assert cell_set_is_remappable(cells, [np.array([1]), np.array([0, 2])])
    assert cell_set_is_remappable(cells, [None, np.array([2])])
    assert cell_set_is_remappable(cells, [np.array([], dtype=np.int64)])
    assert not cell_set_is_remappable(cells, [np.array([2]), None])  # past the end
    assert not cell_set_is_remappable(cells, [np.array([-1]), None])
    assert not cell_set_is_remappable(cells, [np.array([0.0]), None])


def _piece():
    # A piece keeping points 1, 2, 3 (as 0, 1, 2), triangle 1 and lines 1, 2.
    point_map = np.array([-1, 0, 1, 2])
    cell_maps = [np.array([-1, 0]), np.array([-1, 0, 1])]
    return point_map, cell_maps


def _out(blocks=("triangle", "line")):
    # The piece itself: the cell sets are views over its blocks.
    cells = {
        "triangle": ("triangle", np.array([[0, 2, 1]])),
        "line": ("line", np.array([[0, 2], [2, 1]])),
    }
    return Mesh(np.zeros((3, 3)), [cells[b] for b in blocks])


def test_split_mode_skips_side_data():
    src = _mesh()
    out = _out()
    point_map, cell_maps = _piece()
    remap_piece_sets(
        src, out, point_map, cell_maps, keep_invalid=False, drop_empty_blocks=True
    )
    assert set(out.point_sets) == {"corner"}
    assert_array_equal(out.point_sets["corner"], [2])
    assert set(out.cell_sets) == {"sel"}
    assert [list(b) for b in out.cell_sets["sel"]] == [[0], [1]]


def test_partition_mode_copies_side_data_through():
    src = _mesh()
    out = _out()
    point_map, cell_maps = _piece()
    remap_piece_sets(
        src, out, point_map, cell_maps, keep_invalid=True, drop_empty_blocks=False
    )
    assert_array_equal(out.point_sets["tags"], [-7, 2])
    assert [list(b) for b in out.cell_sets["gmsh:bounding_entities"]] == [
        [-1, 5],
        [9],
    ]
    assert [list(b) for b in out.cell_sets["sel"]] == [[0], [1]]


def test_drop_empty_blocks_follows_the_piece():
    src = _mesh()
    out = _out(("triangle",))
    point_map = np.array([0, 1, 2, -1])
    cell_maps = [np.array([0, -1]), np.array([-1, -1, -1])]  # no line kept
    remap_piece_sets(
        src, out, point_map, cell_maps, keep_invalid=False, drop_empty_blocks=True
    )
    assert len(out.cell_sets["sel"]) == 1
    out2 = _out()
    remap_piece_sets(
        src, out2, point_map, cell_maps, keep_invalid=False, drop_empty_blocks=False
    )
    assert len(out2.cell_sets["sel"]) == 2
    assert len(out2.cell_sets["sel"][1]) == 0
