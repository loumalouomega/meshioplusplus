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
