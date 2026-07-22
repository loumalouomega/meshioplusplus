"""Shared fixtures for the data-operation tests.

These operations need geometry with hand-computable answers rather than the
round-trip fixtures in ``helpers.py``, so they build their meshes locally --
the same approach ``test_stats.py`` and ``test_split.py`` take.
"""

import numpy as np

import meshioplusplus as mp

# Two unit squares side by side, split into a triangle block (2 cells) and a
# quad block (1 cell):
#
#   3---2---5
#   |  /|   |
#   | / |   |
#   0---1---4
#
POINTS = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 1.0, 0.0],
    ]
)

TRIANGLES = np.array([[0, 1, 2], [0, 2, 3]])
QUADS = np.array([[1, 4, 5, 2]])


def data_mesh():
    """The shared two-block fixture carrying point, cell and field data."""
    m = mp.Mesh(POINTS.copy(), [("triangle", TRIANGLES.copy()), ("quad", QUADS.copy())])
    # T = x + 10*y at each point.
    m.point_data["T"] = np.array([0.0, 1.0, 11.0, 10.0, 2.0, 12.0])
    m.point_data["v"] = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
        ]
    )
    m.cell_data["mat"] = [np.array([1.0, 2.0]), np.array([3.0])]
    m.cell_data["tag"] = [
        np.array([10, 20], dtype=np.int32),
        np.array([30], dtype=np.int32),
    ]
    m.field_data["meta"] = np.array([1.0, 2.0, 3.0])
    return m


def assert_same_geometry(a, b):
    """Every data operation must leave points and connectivity bit-identical."""
    assert np.array_equal(a.points, b.points)
    assert len(a.cells) == len(b.cells)
    for ca, cb in zip(a.cells, b.cells):
        assert ca.type == cb.type
        assert np.array_equal(np.asarray(ca.data), np.asarray(cb.data))
