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
    return meshioplusplus.clean(
        mesh,
        remove_orphans=False,
        drop_degenerate=False,
        drop_duplicate_cells=False,
    )


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
def test_single_large_row():
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


@needs_core
def test_rows_sharing_one_buffer_are_read_as_spans():
    """Roadmap 5.6: rows that are back-to-back views of one buffer are copied
    as one span; views that skip, overlap or repeat memory, separate arrays,
    int32 rows and lists, interleaved, still give each row its own values."""
    flat = np.arange(40, dtype=np.int64) % 8
    rows = [
        flat[0:3],
        flat[3:7],
        flat[7:10],  # three back to back
        flat[12:15],  # a gap
        flat[11:14],  # backwards, overlapping
        flat[11:14],  # the same view again
        np.array([1, 2, 3, 4], dtype=np.int32),
        [5, 6, 7],
        flat[20:25],
        flat[25:28],
        flat[30:36:2],  # strided: not C-contiguous
    ]
    out = _roundtrip(_polygon_mesh(rows))
    got = [np.asarray(r).tolist() for r in out.cells[0].data]
    assert got == [np.asarray(r).tolist() for r in rows]


# --- rows that carry their CSR arrays (roadmap 5.6) ------------------------

from meshioplusplus._ragged import CsrFaces, CsrRows  # noqa: E402


def _ragged_rows():
    return [
        np.array([0, 1, 2], np.int64),
        np.array([1, 2, 3, 4], np.int64),
        np.array([4, 5, 6], np.int64),
        np.array([2, 3, 5, 6, 7], np.int64),
    ]


def _as_lists(data):
    if data and isinstance(data[0], list):
        return [[f.tolist() for f in cell] for cell in data]
    return [r.tolist() for r in data]


@needs_core
def test_core_rows_carry_their_csr():
    out = _roundtrip(_polygon_mesh(_ragged_rows()))
    data = out.cells[0].data
    assert isinstance(data, CsrRows) and isinstance(data, list)
    flat, offsets, face_offsets = data.csr()
    assert face_offsets is None
    np.testing.assert_array_equal(offsets, [0, 3, 7, 10, 15])
    np.testing.assert_array_equal(flat, np.concatenate(_ragged_rows()))
    for r, row in enumerate(data):
        assert row.base is flat or np.shares_memory(row, flat)
        np.testing.assert_array_equal(row, flat[offsets[r] : offsets[r + 1]])


@needs_core
def test_carried_csr_comes_back_with_writes_into_rows():
    out = _roundtrip(_polygon_mesh(_ragged_rows()))
    out.cells[0].data[1][0] = 7  # a write into a row is a write into `flat`
    again = _roundtrip(out)
    assert isinstance(again.cells[0].data, CsrRows)
    want = [r.tolist() for r in _ragged_rows()]
    want[1][0] = 7
    assert _as_lists(again.cells[0].data) == want


_MUTATIONS = {
    "setitem": lambda d: d.__setitem__(0, np.array([7, 6, 5], np.int64)),
    "setslice": lambda d: d.__setitem__(slice(0, 1), [np.array([7, 6], np.int64)]),
    "delitem": lambda d: d.__delitem__(1),
    "iadd": lambda d: d.__iadd__([np.array([0, 7, 1], np.int64)]),
    "imul": lambda d: d.__imul__(2),
    "append": lambda d: d.append(np.array([0, 7, 1], np.int64)),
    "extend": lambda d: d.extend([np.array([0, 7, 1], np.int64)]),
    "insert": lambda d: d.insert(0, np.array([0, 7, 1], np.int64)),
    "pop": lambda d: d.pop(),
    "remove": lambda d: d.remove(d[0]),  # found by identity first
    "clear": lambda d: d.clear(),
    "sort": lambda d: d.sort(key=len),
    "reverse": lambda d: d.reverse(),
}


@needs_core
@pytest.mark.parametrize("name", sorted(_MUTATIONS))
def test_a_reshaped_block_drops_its_csr_and_round_trips(name):
    out = _roundtrip(_polygon_mesh(_ragged_rows()))
    data = out.cells[0].data
    _MUTATIONS[name](data)
    assert data.csr() is None
    want = _as_lists(data)
    if not want:
        return  # an empty polygon block is not a mesh the core keeps
    again = _roundtrip(Mesh(out.points, [("polygon", data)]))
    assert _as_lists(again.cells[0].data) == want


@needs_core
def test_polyhedron_face_lists_share_the_csr():
    out = _roundtrip(Mesh(_polyhedron_mesh().points, [_polyhedron_block(2)]))
    data = out.cells[0].data
    assert isinstance(data, CsrRows) and all(isinstance(c, CsrFaces) for c in data)
    flat, offsets, face_offsets = data.csr()
    np.testing.assert_array_equal(face_offsets, [0, 6, 12])
    assert len(offsets) == 13 and offsets[-1] == len(flat) == 48
    want = _as_lists(data)
    assert _as_lists(_roundtrip(out).cells[0].data) == want
    # Changing one cell's face list drops the whole block's arrays.
    data[1].pop()
    assert data.csr() is None
    want[1].pop()
    again = _roundtrip(Mesh(out.points, [("polyhedron", data)]))
    assert _as_lists(again.cells[0].data) == want


def _polyhedron_block(copies):
    faces = _polyhedron_mesh().cells[0].data[0]
    return ("polyhedron", [[f.copy() for f in faces] for _ in range(copies)])


@needs_core
def test_copies_and_pickles_keep_rows_on_their_own_buffer():
    import copy
    import pickle

    out = _roundtrip(_polygon_mesh(_ragged_rows()))
    data = out.cells[0].data
    for dup in (copy.deepcopy(data), pickle.loads(pickle.dumps(data))):
        assert isinstance(dup, CsrRows)
        flat = dup.csr()[0]
        assert not np.shares_memory(flat, data.csr()[0])
        dup[0][0] = 6
        assert flat[0] == 6 and data[0][0] == 0
    shallow = copy.copy(data)
    assert isinstance(shallow, CsrRows) and shallow.csr()[0] is data.csr()[0]
    data.append(np.array([1, 2, 3], np.int64))
    dropped = pickle.loads(pickle.dumps(data))
    assert isinstance(dropped, CsrRows) and dropped.csr() is None
    assert _as_lists(dropped) == _as_lists(data)


def test_uniform_rows_become_the_flat_buffer_reshaped():
    flat = np.arange(12, dtype=np.int64)
    rows = CsrRows.from_csr(flat, [0, 4, 8, 12])
    block = Mesh(np.zeros((12, 3)), [("polygon", rows)]).cells[0]
    assert isinstance(block.data, np.ndarray)
    np.testing.assert_array_equal(block.data, np.asarray(list(rows)))
    assert block.data.dtype == np.int64 and block.data.shape == (3, 4)


@needs_core
def test_inconsistent_carried_arrays_are_ignored():
    out = _roundtrip(_polygon_mesh(_ragged_rows()))
    data = out.cells[0].data
    flat, offsets, _ = data.csr()
    want = _as_lists(data)
    # Arrays that do not describe the list are not trusted: the rows are read.
    for forged in (
        (flat, offsets[:-1], None),
        (flat[:-1], offsets, None),
        (flat.astype(np.int32), offsets, None),
        (flat, offsets, offsets),
    ):
        data._state[0] = forged
        assert _as_lists(_roundtrip(out).cells[0].data) == want
