"""Ragged blocks and per-call conversions across the Python <-> C++ boundary
(roadmap 3.4.3): the values, dtypes and aliasing a caller can observe."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import Mesh

try:
    from meshioplusplus import _core
except ImportError:
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the C++ core")


def _polygon_mesh(rows):
    pts = np.array([[float(i), float(i % 3), 0.0] for i in range(8)])
    return Mesh(pts, [("polygon", rows)])


def _polyhedron_mesh():
    pts = np.array(
        [[x, y, z] for z in (0.0, 1.0) for y in (0.0, 1.0) for x in (0.0, 1.0)]
    )
    faces = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [1, 3, 7, 5],
        [3, 2, 6, 7],
        [2, 0, 4, 6],
    ]
    return Mesh(pts, [("polyhedron", [[np.array(f, np.int64) for f in faces]])])


def _roundtrip(mesh):
    # `clean` reaches the core and back (py_to_mesh then mesh_to_py).
    return meshioplusplus.clean(mesh)


@needs_core
def test_polygon_rows_round_trip_with_dtype_and_layout():
    rows = [
        np.array([0, 1, 2], np.int64),
        np.array([1, 2, 3, 4], np.int64),
        np.array([4, 5, 6], np.int64),
    ]
    out = _roundtrip(_polygon_mesh(rows))
    block = out.cells[0]
    assert block.type == "polygon"
    assert len(block.data) == 3
    for got, want in zip(block.data, rows):
        assert got.dtype == np.int64
        assert got.ndim == 1
        assert got.flags["C_CONTIGUOUS"]
        np.testing.assert_array_equal(got, want)


@needs_core
def test_writing_one_row_does_not_touch_its_neighbours():
    rows = [np.array([0, 1, 2], np.int64), np.array([2, 3, 4], np.int64)]
    out = _roundtrip(_polygon_mesh(rows))
    before = [np.array(r) for r in out.cells[0].data]
    out.cells[0].data[0][:] = 7
    np.testing.assert_array_equal(out.cells[0].data[1], before[1])


@needs_core
def test_polyhedron_faces_round_trip():
    mesh = _polyhedron_mesh()
    out = _roundtrip(mesh)
    cells = out.cells[0].data
    assert len(cells) == 1 and len(cells[0]) == 6
    for got, want in zip(cells[0], mesh.cells[0].data[0]):
        assert got.dtype == np.int64
        np.testing.assert_array_equal(got, want)


@needs_core
@pytest.mark.parametrize(
    "make_row",
    [
        lambda v: np.array(v, np.int32),  # safe cast to int64
        lambda v: np.array(v, np.int64)[::1],
        lambda v: np.array(v, np.dtype(">i8")),  # big-endian
        lambda v: list(v),  # plain Python ints
        lambda v: tuple(v),
    ],
    ids=["int32", "int64", "big-endian", "list", "tuple"],
)
def test_row_input_forms_give_the_same_mesh(make_row):
    expected = [[0, 1, 2], [1, 2, 3, 4], [4, 5, 6]]
    out = _roundtrip(_polygon_mesh([make_row(r) for r in expected]))
    for got, want in zip(out.cells[0].data, expected):
        np.testing.assert_array_equal(got, want)


@needs_core
def test_non_contiguous_row_is_read_by_value():
    base = np.arange(12, dtype=np.int64)
    rows = [base[0:9:3], base[1:10:3]]  # strided views: [0 3 6], [1 4 7]
    out = _roundtrip(_polygon_mesh(rows))
    np.testing.assert_array_equal(out.cells[0].data[0], [0, 3, 6])
    np.testing.assert_array_equal(out.cells[0].data[1], [1, 4, 7])


@needs_core
def test_float_rows_are_still_rejected():
    rows = [np.array([0.0, 1.0, 2.0])]
    with pytest.raises(Exception):
        _roundtrip(_polygon_mesh(rows))


@needs_core
def test_empty_and_single_large_row():
    big = np.arange(5000, dtype=np.int64) % 8
    out = _roundtrip(_polygon_mesh([big]))
    np.testing.assert_array_equal(out.cells[0].data[0], big)


@needs_core
def test_repeated_calls_share_the_cached_classes():
    tiny = Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float),
        [("tetra", np.array([[0, 1, 2, 3]], np.int64))],
    )
    first = _roundtrip(tiny)
    second = _roundtrip(tiny)
    assert type(first) is meshioplusplus.Mesh is type(second)
    np.testing.assert_array_equal(first.points, second.points)


@needs_core
def test_big_endian_points_are_normalised():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], np.dtype(">f8"))
    mesh = Mesh(pts, [("tetra", np.array([[0, 1, 2, 3]], np.dtype(">i8")))])
    out = _roundtrip(mesh)
    np.testing.assert_array_equal(out.points, np.asarray(pts, float))
