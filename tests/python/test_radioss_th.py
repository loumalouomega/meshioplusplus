"""OpenRadioss time-history files (``<run>T01``): both engines on a T01 that
OpenRadioss wrote for a meshio++ deck, checked against OpenRadioss's own
``th_to_csv`` reading, and on synthetic files for the layout's other branches."""

import pathlib
import shutil
import struct

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.radioss_th import _th as py_th

COLUMN = pathlib.Path(__file__).parent / "meshes" / "radioss_th" / "column"
T01 = COLUMN / "columnT01"
REFERENCE = COLUMN / "columnT01_th_to_csv.csv"


@pytest.fixture(params=["core", "python"])
def engine(request):
    if request.param == "core":
        return meshioplusplus.radioss_th
    return py_th


def _record(payload):
    return struct.pack(">I", len(payload)) + payload + struct.pack(">I", len(payload))


def _ints(*values):
    return _record(struct.pack(f">{len(values)}i", *values))


def _floats(*values):
    return _record(struct.pack(f">{len(values)}f", *values))


def _text(text, width):
    return text.ljust(width).encode()


def th_file(version=3040, outputs=2, cut=True):
    """Globals 1, 2; parts 7 (IE, KE) and 8 (none); subset 0 (KE); TH groups:
    nodes 11 and 12 (three variables), rigid body 5 (none)."""
    title = 100 if version >= 4021 else 80 if version >= 3041 else 40
    out = _record(struct.pack(">i", version) + _text("run", 80))
    out += _record(_text("Mon Mar 30 16:39:04 2026 RADIOSS", 80))
    if version > 3050:
        out += _ints(1) + _ints(title) + _floats(1.0, 1000.0, 0.001)
    out += _ints(2, 1, 0, 1, 2, 2) + _ints(1, 2)
    out += _record(
        struct.pack(">i", 7) + _text("part7", title) + struct.pack(">4i", 0, 0, 0, 2)
    )
    out += _ints(1, 2)
    out += _record(
        struct.pack(">i", 8) + _text("part8", title) + struct.pack(">4i", 0, 0, 0, 0)
    )
    out += _record(struct.pack(">i", 1) + _text("steel", title))
    out += _record(struct.pack(">5i", 0, 0, 0, 2, 1) + _text("GLOBAL MODEL", title))
    out += _ints(7, 8) + _ints(2)
    out += _record(struct.pack(">5i", 3, 0, 0, 2, 3) + _text("th_nodes", title))
    out += _record(struct.pack(">i", 11) + _text("a", title))
    out += _record(struct.pack(">i", 12) + _text("b", title))
    out += _ints(1, 2, 3)
    out += _record(struct.pack(">5i", 4, 103, 0, 1, 0) + _text("TH RBODY", title))
    out += _record(struct.pack(">i", 5) + _text("body", title))
    for k in range(outputs):
        out += _floats(0.5 * k) + _floats(10 + k, 20 + k) + _floats(1 + k, 2 + k)
        out += (
            _floats(3 + k) + _floats(*[i + 100.0 * k for i in range(6)]) + _record(b"")
        )
    if cut:
        out += _floats(0.5 * outputs) + _floats(1.0)
    return out


def test_engines_agree():
    for step in range(len(py_th.time_values(T01))):
        a = meshioplusplus.radioss_th.read(T01, time_step=step)
        b = py_th.read(T01, time_step=step)
        assert sorted(a.field_data) == sorted(b.field_data)
        for name in a.field_data:
            assert a.field_data[name].dtype == b.field_data[name].dtype
            np.testing.assert_array_equal(a.field_data[name], b.field_data[name])


def test_real_run_matches_th_to_csv(engine):
    """Every output of OpenRadioss's T01, value for value, against th_to_csv
    (which prints 7 significant digits)."""
    reference = np.loadtxt(REFERENCE, delimiter=",")
    times = engine.time_values(T01)
    assert len(times) == len(reference) == 8
    np.testing.assert_allclose(times, reference[:, 0], rtol=1e-6)
    for step, row in enumerate(reference):
        mesh = engine.read(T01, time_step=step)
        fd = mesh.field_data
        assert len(mesh.points) == 0
        # the time, the 22 global variables (codes 1 to 22) in file order
        values = [fd["meshio:time"][0]]
        values += [
            fd[f"radioss_th:global:{py_th._global_name(c)}"][0] for c in range(1, 23)
        ]
        # part 1, then the TH groups in file order: bricks, then nodes
        tail = [fd[f"radioss_th:part:1:{v}"][0] for v in ("IE", "KE", "MASS")]
        tail += list(fd["radioss_th:brick:2"].ravel())
        tail += list(fd["radioss_th:node:1"].ravel())
        np.testing.assert_allclose(values + tail, row, rtol=1e-6, atol=1e-30)
    last = engine.read(T01, time_step=-1).field_data
    np.testing.assert_array_equal(last["radioss_th:node:1:ids"], [1, 9, 12])
    np.testing.assert_array_equal(
        last["radioss_th:node:1:variables"], [1, 2, 3, 4, 5, 6]
    )
    np.testing.assert_array_equal(last["radioss_th:brick:2:ids"], [1, 2])
    # node 1 is held
    assert not last["radioss_th:node:1"][0].any()


def test_synthetic_layout(engine, tmp_path):
    path = tmp_path / "runT01"
    path.write_bytes(th_file())
    assert engine.time_values(path) == [0.0, 0.5]  # the cut third output is dropped
    fd = engine.read(path, time_step=-1).field_data
    assert fd["radioss_th:global:internal_energy"][0] == 11.0
    assert fd["radioss_th:part:7:KE"][0] == 3.0
    assert fd["radioss_th:subset:0:KE"][0] == 4.0
    assert fd["radioss_th:node:3"].shape == (2, 3)
    assert fd["radioss_th:node:3"][1, 2] == 105.0
    assert fd["radioss_th:rbody:4"].shape == (1, 0)
    assert "radioss_th:unit_factors" not in fd
    newer = tmp_path / "newT01"
    newer.write_bytes(th_file(version=4021, cut=False))
    fd = engine.read(newer).field_data
    np.testing.assert_allclose(
        fd["radioss_th:unit_factors"], [1.0, 1000.0, 0.001], rtol=1e-7
    )
    assert fd["radioss_th:part:7:KE"][0] == 2.0


def test_narrowed(engine):
    mesh = engine.read(T01, arrays=["radioss_th:part:1:KE"])
    assert set(mesh.field_data) == {"meshio:time", "radioss_th:part:1:KE"}
    bare = engine.read(T01, points_only=True)
    assert set(bare.field_data) == {"meshio:time"}


def test_found_by_name_and_by_content(tmp_path):
    assert meshioplusplus.sniff_format(T01) == "radioss_th"
    renamed = tmp_path / "history.bin"
    shutil.copy(T01, renamed)
    assert meshioplusplus.sniff_format(renamed) == "radioss_th"
    assert _core.sniff_format(str(renamed)) == "radioss_th"
    mesh = meshioplusplus.read(T01, time_step=3)
    assert mesh.field_data["meshio:time"][0] == pytest.approx(
        meshioplusplus.radioss_th.time_values(T01)[3]
    )
    assert meshioplusplus.read_metadata(T01)["time_values"] == pytest.approx(
        meshioplusplus.radioss_th.time_values(T01)
    )
    steps = list(meshioplusplus.read_sequence(T01))
    assert len(steps) == 8


def test_refusals(engine, tmp_path):
    junk = tmp_path / "junkT01"
    junk.write_bytes(b"\x00" * 64)
    with pytest.raises(meshioplusplus.ReadError):
        engine.read(junk)
    with pytest.raises(meshioplusplus.ReadError, match="out of range"):
        engine.read(T01, time_step=8)
    cut = tmp_path / "cutT01"
    cut.write_bytes(th_file(outputs=0))
    with pytest.raises(meshioplusplus.ReadError, match="no complete output"):
        engine.read(cut)
