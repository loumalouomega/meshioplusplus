"""Tecplot SZL (``.szplt``) through TecIO: every fixture was written by TecIO
beside an ASCII twin holding the same data (``meshes/szplt/README.md``), and
the two must read identically.

The core engine needs a build with ``MESHIOPLUSPLUS_WITH_TECIO=ON``; the ctypes
twin needs ``MESHIOPLUSPLUS_TECIO_LIBRARY`` naming a shared TecIO. Each is
skipped when absent."""

import os
import pathlib
import sys

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._exceptions import ReadError
from meshioplusplus.szplt import _szplt

SZPLT = pathlib.Path(__file__).parent / "meshes" / "szplt"
NAMES = ["fe_mixed", "ordered", "transient", "triangles"]


def _core_read(path, time_step=0):
    return meshioplusplus.szplt.read(path, time_step=time_step)


@pytest.fixture(params=["core", "ctypes"])
def engine(request):
    if request.param == "core":
        if not getattr(_core, "__has_tecio__", False):
            pytest.skip("core built without TecIO")
        return _core_read, meshioplusplus.szplt.time_values
    if not os.environ.get("MESHIOPLUSPLUS_TECIO_LIBRARY"):
        pytest.skip("MESHIOPLUSPLUS_TECIO_LIBRARY names no shared TecIO")
    return _szplt.read, _szplt.time_values


def _assert_same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    assert sorted(a.point_data) == sorted(b.point_data)
    for key in a.point_data:
        np.testing.assert_array_equal(a.point_data[key], b.point_data[key])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for key in a.cell_data:
        for x, y in zip(a.cell_data[key], b.cell_data[key]):
            np.testing.assert_array_equal(x, y)
    assert sorted(r.name for r in a.regions) == sorted(r.name for r in b.regions)


@pytest.mark.parametrize("name", NAMES)
def test_reads_as_its_ascii_twin(engine, name):
    read, time_values = engine
    times = meshioplusplus.read_metadata(SZPLT / f"{name}.dat")["time_values"]
    assert time_values(SZPLT / f"{name}.szplt") == times
    for step in range(max(len(times), 1)):
        _assert_same(
            read(SZPLT / f"{name}.szplt", time_step=step),
            meshioplusplus.read(SZPLT / f"{name}.dat", time_step=step),
        )


def test_out_of_range_step(engine):
    read, _ = engine
    with pytest.raises(ReadError, match="3 steps"):
        read(SZPLT / "transient.szplt", time_step=3)


def test_dispatch_metadata_and_sniff():
    if not getattr(_core, "__has_tecio__", False):
        pytest.skip("core built without TecIO")
    mesh = meshioplusplus.read(SZPLT / "fe_mixed.szplt")
    assert [c.type for c in mesh.cells] == ["tetra", "hexahedron"]
    meta = meshioplusplus.read_metadata(SZPLT / "transient.szplt")
    assert meta["time_values"] == [0.0, 0.5, 1.0]
    from meshioplusplus._sniff import sniff_format

    assert sniff_format(SZPLT / "fe_mixed.szplt") == "szplt"


def test_without_tecio_names_the_remedies(monkeypatch):
    monkeypatch.delenv("MESHIOPLUSPLUS_TECIO_LIBRARY", raising=False)
    monkeypatch.setattr(_szplt, "_LIB", None)
    monkeypatch.setattr(_szplt.ctypes.util, "find_library", lambda name: None)
    with pytest.raises(ImportError, match=r"save the file as \.plt"):
        _szplt.read(SZPLT / "fe_mixed.szplt")
    assert "ctypes" in sys.modules


# ---------------------------------------------------------------------------
# The writer: the zones of the ASCII and `.plt` writers, through TecIO. The
# `.dat` twins of the fixtures are the oracle; each engine writes and reads.
# ---------------------------------------------------------------------------

from meshioplusplus.szplt import _write as py_write  # noqa: E402


@pytest.fixture(params=["core", "ctypes"])
def writer(request):
    if request.param == "core":
        if not getattr(_core, "__has_tecio__", False):
            pytest.skip("core built without TecIO")
        return lambda path, mesh: _core.szplt_write(str(path), mesh)
    if not os.environ.get("MESHIOPLUSPLUS_TECIO_LIBRARY"):
        pytest.skip("MESHIOPLUSPLUS_TECIO_LIBRARY names no shared TecIO")
    return py_write.write


def _any_reader():
    if getattr(_core, "__has_tecio__", False):
        return _core_read
    return _szplt.read


@pytest.mark.parametrize("name", ["fe_mixed", "triangles"])
def test_a_written_file_reads_as_the_ascii_twin(writer, name, tmp_path):
    twin = meshioplusplus.read(SZPLT / f"{name}.dat")
    out = tmp_path / f"{name}.szplt"
    writer(out, twin)
    again = _any_reader()(out)
    _assert_same(again, twin)


def test_wedges_and_pyramids_are_degenerate_bricks(writer, tmp_path):
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [0, 1, 1],
            [0.3, 0.3, 2],
        ],
        dtype=float,
    )
    mesh = meshioplusplus.Mesh(
        pts,
        [
            ("wedge", np.array([[0, 1, 2, 3, 4, 5]])),
            ("pyramid", np.array([[0, 1, 4, 3, 6]])),
        ],
    )
    mesh.point_data["p"] = np.arange(7.0)
    out = tmp_path / "w.szplt"
    writer(out, mesh)
    again = _any_reader()(out)
    # Degenerate bricks, as the `.plt` writer's: SZL has no wedge or pyramid.
    assert [c.type for c in again.cells] == ["hexahedron", "hexahedron"]
    rows = [sorted(set(row)) for c in again.cells for row in c.data.tolist()]
    assert [len(r) for r in rows] == [6, 5]
    assert len(again.points) == 11  # each block's own points (shared ones repeated)


def test_a_series_is_a_strand_of_zones(tmp_path):
    if not (
        getattr(_core, "__has_tecio__", False)
        or os.environ.get("MESHIOPLUSPLUS_TECIO_LIBRARY")
    ):
        pytest.skip("no TecIO")
    base = meshioplusplus.read(SZPLT / "fe_mixed.dat")
    steps = []
    for k in range(3):
        mesh = meshioplusplus.read(SZPLT / "fe_mixed.dat")
        mesh.point_data["P"] = mesh.point_data["P"] + k
        steps.append((0.5 * k, mesh))
    out = tmp_path / "series.szplt"
    meshioplusplus.write_sequence(out, iter(steps), file_format="szplt")
    assert meshioplusplus.szplt.time_values(out) == [0.0, 0.5, 1.0]
    for k, (_, mesh) in enumerate(steps):
        again = _any_reader()(out, time_step=k)
        np.testing.assert_array_equal(again.points, base.points)
        np.testing.assert_array_equal(again.point_data["P"], mesh.point_data["P"])
        assert [c.type for c in again.cells] == [c.type for c in base.cells]


def test_a_series_refuses_a_different_mesh(tmp_path):
    if not os.environ.get("MESHIOPLUSPLUS_TECIO_LIBRARY"):
        pytest.skip("MESHIOPLUSPLUS_TECIO_LIBRARY names no shared TecIO")
    mesh = meshioplusplus.read(SZPLT / "fe_mixed.dat")
    other = meshioplusplus.read(SZPLT / "fe_mixed.dat")
    other.point_data["extra"] = np.zeros(len(other.points))
    with pytest.raises(meshioplusplus.WriteError):
        with py_write.SeriesWriter(tmp_path / "bad.szplt") as w:
            w.write(0.0, mesh)
            w.write(1.0, other)


# ---------------------------------------------------------------------------
# Variable types and partitioned zones.
# ---------------------------------------------------------------------------


def _variable_types(path):
    """TecIO's type of every variable of the first zone (1 float, 2 double, 3 int32,
    4 int16, 5 byte), asked through the shared TecIO."""
    import ctypes

    lib = _szplt.library()
    handle = ctypes.c_void_p()
    assert lib.tecFileReaderOpen(str(path).encode(), ctypes.byref(handle)) == 0
    try:
        count = ctypes.c_int32()
        assert lib.tecDataSetGetNumVars(handle, ctypes.byref(count)) == 0
        out = []
        for var in range(1, count.value + 1):
            kind = ctypes.c_int32()
            assert lib.tecZoneVarGetType(handle, 1, var, ctypes.byref(kind)) == 0
            out.append(kind.value)
        return out
    finally:
        lib.tecFileReaderClose(ctypes.byref(handle))


def _typed_mesh():
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=np.float32
    )
    mesh = meshioplusplus.Mesh(pts, [("tetra", np.array([[0, 1, 2, 3], [1, 2, 3, 4]]))])
    mesh.point_data["f32"] = np.arange(5, dtype=np.float32) + 0.25
    mesh.point_data["i16"] = np.array([-3, 0, 7, 100, 30000], dtype=np.int16)
    mesh.point_data["i8"] = np.array([-3, 0, 7, 100, 127], dtype=np.int8)
    mesh.point_data["u8"] = np.array([0, 1, 2, 3, 255], dtype=np.uint8)
    mesh.point_data["i32"] = np.array([-1, 0, 1, 2, 2**30], dtype=np.int32)
    mesh.point_data["i64"] = np.arange(5, dtype=np.int64)
    mesh.cell_data["c32"] = [np.array([1.5, 2.5], dtype=np.float32)]
    return mesh


def test_variables_keep_the_type_that_holds_their_dtype(writer, tmp_path):
    if not os.environ.get("MESHIOPLUSPLUS_TECIO_LIBRARY"):
        pytest.skip("MESHIOPLUSPLUS_TECIO_LIBRARY names no shared TecIO")
    mesh = _typed_mesh()
    out = tmp_path / "typed.szplt"
    writer(out, mesh)
    names = ["X", "Y", "Z"] + sorted(mesh.point_data) + ["c32"]
    assert dict(zip(names, _variable_types(out))) == {
        "X": 1,
        "Y": 1,
        "Z": 1,
        "f32": 1,
        "i16": 4,
        "i8": 4,
        "u8": 5,
        "i32": 3,
        "i64": 2,
        "c32": 1,
    }
    again = _any_reader()(out)
    np.testing.assert_array_equal(again.points, mesh.points.astype(float))
    for name in mesh.point_data:
        np.testing.assert_array_equal(
            again.point_data[name], mesh.point_data[name].astype(float)
        )
    np.testing.assert_array_equal(again.cell_data["c32"][0], [1.5, 2.5])


def _partitioned_mesh(kind, parts=(0, 0, 1, 1)):
    # a row of cubes along x, two parts: the first two cubes, the last two
    pts, cells = [], []
    for i in range(5):
        for corner in ([0, 0, 0], [0, 1, 0], [0, 1, 1], [0, 0, 1]):
            pts.append([i, corner[1], corner[2]])
    for i in range(4):
        a, b = i * 4, (i + 1) * 4
        cells.append(list(range(a, a + 4)) + list(range(b, b + 4)))
    mesh = meshioplusplus.Mesh(
        np.array(pts, dtype=float), [("hexahedron", np.array(cells))]
    )
    mesh.point_data["xx"] = np.array(pts, dtype=float)[:, 0] * 10
    mesh.cell_data["partition:part"] = [np.array(parts)]
    mesh.cell_data["rho"] = [np.array([1.0, 2.0, 3.0, 4.0])]
    return mesh


def _cell_fingerprints(mesh):
    """Each cell as its sorted node coordinates with its data, in any order."""
    out = []
    for block in mesh.cells:
        for row in range(len(block.data)):
            pts = tuple(sorted(map(tuple, mesh.points[block.data[row]].tolist())))
            out.append((pts, float(mesh.cell_data["rho"][0][row])))
    return sorted(out)


@pytest.mark.parametrize("parts", [(0, 0, 1, 1), (1, 0, 1, 0), (0, 1, 2, 2)])
def test_cells_with_parts_are_written_as_partitions(writer, parts, tmp_path):
    mesh = _partitioned_mesh("hexahedron", parts)
    out = tmp_path / "parts.szplt"
    writer(out, mesh)
    again = _any_reader()(out)
    # the partitions come back as one zone: the nodes numbered partition by partition (each
    # one's owned nodes in order), the cells grouped by part, every cell's geometry intact
    assert len(again.points) == len(mesh.points)
    assert _cell_fingerprints(again) == _cell_fingerprints(mesh)
    np.testing.assert_array_equal(again.cell_data["partition:part"][0], sorted(parts))
    by_point = {
        tuple(p): v for p, v in zip(again.points.tolist(), again.point_data["xx"])
    }
    for p, v in zip(mesh.points.tolist(), mesh.point_data["xx"]):
        assert by_point[tuple(p)] == v


def test_a_part_without_a_node_of_its_own_is_written_whole(writer, tmp_path, capfd):
    # parts 0 1 2 1: the third part's nodes all belong to lower parts, which TecIO refuses
    mesh = _partitioned_mesh("hexahedron", (0, 1, 2, 1))
    out = tmp_path / "whole.szplt"
    writer(out, mesh)
    assert "owns no node of its own" in capfd.readouterr().err
    again = _any_reader()(out)
    np.testing.assert_array_equal(again.points, mesh.points)
    np.testing.assert_array_equal(again.cell_data["partition:part"][0], [0, 1, 2, 1])
    assert _cell_fingerprints(again) == _cell_fingerprints(mesh)


def test_the_partitions_of_two_engines_agree(tmp_path):
    if not (
        getattr(_core, "__has_tecio__", False)
        and os.environ.get("MESHIOPLUSPLUS_TECIO_LIBRARY")
    ):
        pytest.skip("needs both engines")
    mesh = _partitioned_mesh("hexahedron")
    a, b = tmp_path / "core.szplt", tmp_path / "twin.szplt"
    _core.szplt_write(str(a), mesh)
    py_write.write(b, mesh)
    x, y = _szplt.read(a), _szplt.read(b)
    np.testing.assert_array_equal(x.points, y.points)
    assert [c.data.tolist() for c in x.cells] == [c.data.tolist() for c in y.cells]
    assert a.read_bytes() == b.read_bytes()


def test_a_series_is_written_whole_even_with_parts(tmp_path):
    if not (
        getattr(_core, "__has_tecio__", False)
        or os.environ.get("MESHIOPLUSPLUS_TECIO_LIBRARY")
    ):
        pytest.skip("no TecIO")
    steps = [
        (0.0, _partitioned_mesh("hexahedron")),
        (1.0, _partitioned_mesh("hexahedron")),
    ]
    out = tmp_path / "series.szplt"
    meshioplusplus.write_sequence(out, iter(steps), file_format="szplt")
    again = _any_reader()(out, time_step=1)
    np.testing.assert_array_equal(again.points, steps[1][1].points)
