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
