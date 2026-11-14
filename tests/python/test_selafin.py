"""TELEMAC Selafin (roadmap 2.4.3.4): files built here from the record layout."""

import struct

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _fallback
from meshioplusplus._fortran_records import fortran_record
from meshioplusplus.selafin import _selafin

L = (4, ">")


def _slf(path, *, double=False, steps=2, nplan=0, ndp=3, date=False):
    dt = ">f8" if double else ">f4"
    if ndp == 6:  # two planes of a triangle: one prism
        x = np.array([0.0, 1.0, 0.0, 0.0, 1.0, 0.0])
        y = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0])
        ikle = np.array([[1, 2, 3, 4, 5, 6]])
        names = ["ELEVATION Z", "VELOCITY U"]
    else:
        x = np.array([0.0, 1.0, 1.0, 0.0])
        y = np.array([0.0, 0.0, 1.0, 1.0])
        ikle = np.array([[1, 2, 3], [1, 3, 4]])
        names = ["WATER DEPTH", "VELOCITY U"]
    npoin = len(x)
    title = b"test case".ljust(72) + (b"SERAFIND" if double else b"SERAFIN ")
    out = [
        fortran_record(title, L),
        fortran_record(struct.pack(">2i", len(names), 0), L),
    ]
    for n, u in zip(names, ["M", "M/S"]):
        out.append(fortran_record(n.ljust(16).encode() + u.ljust(16).encode(), L))
    out.append(
        fortran_record(struct.pack(">10i", 1, 0, 0, 0, 0, 0, nplan, 0, 0, int(date)), L)
    )
    if date:
        out.append(fortran_record(struct.pack(">6i", 2024, 1, 2, 3, 4, 5), L))
    out.append(
        fortran_record(struct.pack(">4i", len(ikle), npoin, ikle.shape[1], 1), L)
    )
    out.append(fortran_record(ikle.astype(">i4").tobytes(), L))
    out.append(fortran_record(np.arange(npoin).astype(">i4").tobytes(), L))
    out.append(fortran_record(x.astype(dt).tobytes(), L))
    out.append(fortran_record(y.astype(dt).tobytes(), L))
    for s in range(steps):
        out.append(fortran_record(np.array([100.0 * s]).astype(dt).tobytes(), L))
        for v in range(len(names)):
            vals = np.arange(npoin) + 10 * v + 1000 * s
            if ndp == 6 and v == 0:
                vals = np.array([0.0, 0.0, 0.0, 2.0, 2.0, 2.0]) + s
            out.append(fortran_record(np.asarray(vals, float).astype(dt).tobytes(), L))
    path.write_bytes(b"".join(out))


@pytest.fixture(params=["core", "python"])
def reader(request, monkeypatch):
    if request.param == "core":
        monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
        monkeypatch.setattr(_fallback, "_strict", True)
        return mio.selafin.read
    return _selafin.read


@pytest.mark.parametrize("double", [False, True])
@pytest.mark.parametrize("date", [False, True])
def test_steps_and_precision(tmp_path, reader, double, date):
    path = tmp_path / "r2d.slf"
    _slf(path, double=double, date=date)
    mesh = reader(path)
    assert mesh.cells[0].type == "triangle"
    assert np.asarray(mesh.cells[0].data).tolist() == [[0, 1, 2], [0, 2, 3]]
    assert np.allclose(mesh.point_data["WATER DEPTH"], [0, 1, 2, 3])
    assert np.allclose(mesh.point_data["VELOCITY U"], [10, 11, 12, 13])
    assert bytes(mesh.field_data["selafin:unit:VELOCITY U"]) == b"M/S"
    assert bytes(mesh.field_data["selafin:title"]).startswith(b"test case")
    assert np.allclose(mesh.field_data["meshio:time"], [0.0])
    last = reader(path, time_step=-1)
    assert np.allclose(last.point_data["WATER DEPTH"], [1000, 1001, 1002, 1003])
    assert np.allclose(last.field_data["meshio:time"], [100.0])
    with pytest.raises(mio.ReadError, match="out of range"):
        reader(path, time_step=3)


def test_time_values_and_a_truncated_last_step(tmp_path):
    path = tmp_path / "r.slf"
    _slf(path, steps=3)
    data = path.read_bytes()
    path.write_bytes(data[:-10])  # a solver killed mid-write
    assert mio.selafin.time_values(path) == [0.0, 100.0]
    assert _selafin.time_values(path) == [0.0, 100.0]


def test_three_d_prisms_take_their_elevation(tmp_path, reader):
    path = tmp_path / "r3d.slf"
    _slf(path, ndp=6, nplan=2)
    mesh = reader(path, time_step=1)
    assert mesh.cells[0].type == "wedge"
    assert np.allclose(mesh.points[:, 2], [1, 1, 1, 3, 3, 3])


def test_write_round_trips_with_equal_bytes(tmp_path):
    path = tmp_path / "r.slf"
    _slf(path, double=True)
    mesh = mio.selafin.read(path, time_step=1)
    a, b = tmp_path / "a.slf", tmp_path / "b.slf"
    mio.selafin.write(a, mesh)
    _selafin.write(b, mesh)
    assert a.read_bytes() == b.read_bytes()
    back = mio.read(a)
    assert np.allclose(back.point_data["WATER DEPTH"], mesh.point_data["WATER DEPTH"])
    assert np.allclose(back.field_data["meshio:time"], [100.0])
    assert np.asarray(back.point_data["selafin:ipobo"]).tolist() == [0, 1, 2, 3]


def test_garbage_and_wrong_cells(tmp_path, reader):
    bad = tmp_path / "bad.slf"
    bad.write_bytes(b"\x00\x00\x00\x04abcd\x00\x00\x00\x04")
    with pytest.raises(mio.ReadError, match="selafin"):
        reader(bad)
    tets = mio.Mesh(np.eye(4, 3), [("tetra", np.array([[0, 1, 2, 3]]))])
    for write in (mio.selafin.write, _selafin.write):
        with pytest.raises(mio.WriteError, match="triangle or quad"):
            write(tmp_path / "t.slf", tets)
