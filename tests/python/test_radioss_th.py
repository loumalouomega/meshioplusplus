"""OpenRadioss time-history files (``<run>T01``): both engines on a T01 that
OpenRadioss wrote for a meshio++ deck, checked against OpenRadioss's own
``th_to_csv`` reading, and on synthetic files for the layout's other branches.

``oracle_file`` is the hand-framed writer the reader's fixtures were first built
with. It stays as the independent oracle: it builds the layout record by record
(with the versions, materials and subset parts the library writer never emits),
so it checks the reader's branches and, with ``expected_bytes``, the library
writer's bytes without going through the framer the writer uses. Everything else
writes with the library writer (``write``, ``SeriesWriter``).
"""

import io
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


def oracle_file(version=3040, outputs=2, cut=True):
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
    path.write_bytes(oracle_file())
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
    newer.write_bytes(oracle_file(version=4021, cut=False))
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
    cut.write_bytes(oracle_file(outputs=0))
    with pytest.raises(meshioplusplus.ReadError, match="no complete output"):
        engine.read(cut)


# --------------------------------------------------------------------------- #
# The writer (roadmap: result files are written too)                          #
# --------------------------------------------------------------------------- #

WRITERS = {
    "core": lambda p, m: _core.radioss_th_write(str(p), m),
    "python": lambda p, m: py_th.write(p, m),
}
SERIES = {
    "core": lambda p: _core.RadiossThSeriesWriter(str(p)),
    "python": lambda p: py_th.SeriesWriter(p),
}


@pytest.fixture(params=sorted(WRITERS))
def writer(request):
    return request.param


def sample_mesh(k=0, factors=False):
    """The oracle file's content as field data: globals 1, 2; part 7 (IE, KE);
    subset 0 (KE); TH groups node 3 (ids 11, 12; variables 1, 2, 3) and rigid
    body 4 (one entity, no variables)."""
    fd = {
        "meshio:time": np.array([0.5 * k]),
        "radioss_th:global:internal_energy": np.array([10.0 + k]),
        "radioss_th:global:kinetic_energy": np.array([20.0 + k]),
        "radioss_th:part:7:IE": np.array([1.0 + k]),
        "radioss_th:part:7:KE": np.array([2.0 + k]),
        "radioss_th:subset:0:KE": np.array([3.0 + k]),
        "radioss_th:node:3": np.arange(6.0).reshape(2, 3) + 100.0 * k,
        "radioss_th:node:3:ids": np.array([11, 12]),
        "radioss_th:node:3:variables": np.array([1, 2, 3]),
        "radioss_th:rbody:4": np.zeros((1, 0)),
        "radioss_th:rbody:4:ids": np.array([5]),
    }
    if factors:
        fd["radioss_th:unit_factors"] = np.array([1.0, 1000.0, 0.001])
    return meshioplusplus.Mesh(
        np.zeros((0, 3)), [("vertex", np.zeros((0, 1), dtype=np.int64))], field_data=fd
    )


def expected_bytes(steps, factors=False):
    """What the writer must produce for ``sample_mesh``: the layout record by
    record, framed by ``_record`` rather than by the library's framer."""
    title = 100 if factors else 40
    tag = f"Written by meshio++ v{meshioplusplus.__version__}"
    out = _record(struct.pack(">i", 4021 if factors else 3040) + _text(tag, 80))
    out += _record(_text("meshio++", 80))
    if factors:
        out += _ints(1) + _ints(title) + _floats(1.0, 1000.0, 0.001)
    out += _ints(1, 0, 0, 1, 2, 2) + _ints(1, 2)
    out += _record(
        struct.pack(">i", 7) + _text("PART 7", title) + struct.pack(">4i", 0, 0, 0, 2)
    )
    out += _ints(1, 2)
    out += _record(struct.pack(">5i", 0, 0, 0, 0, 1) + _text("SUBSET 0", title))
    out += _ints(2)
    out += _record(struct.pack(">5i", 3, 0, 0, 2, 3) + _text("TH node 3", title))
    out += _record(struct.pack(">i", 11) + _text("node 11", title))
    out += _record(struct.pack(">i", 12) + _text("node 12", title))
    out += _ints(1, 2, 3)
    out += _record(struct.pack(">5i", 4, 103, 0, 1, 0) + _text("TH rbody 4", title))
    out += _record(struct.pack(">i", 5) + _text("rbody 5", title))
    for k in steps:
        out += _floats(0.5 * k) + _floats(10 + k, 20 + k) + _floats(1 + k, 2 + k)
        out += _floats(3 + k) + _floats(*[i + 100.0 * k for i in range(6)])
        out += _record(b"")
    return out


def _same_fields(a, b):
    assert sorted(a.field_data) == sorted(b.field_data)
    for name in a.field_data:
        assert a.field_data[name].dtype == b.field_data[name].dtype
        np.testing.assert_array_equal(a.field_data[name], b.field_data[name])


@pytest.mark.parametrize("factors", [False, True])
def test_the_writer_follows_the_hand_framed_layout(writer, factors, tmp_path):
    path = tmp_path / "runT01"
    WRITERS[writer](path, sample_mesh(1, factors))
    assert path.read_bytes() == expected_bytes([1], factors)


def test_the_oracle_and_the_writer_read_back_alike(engine, writer, tmp_path):
    """The hand-framed file (version 3040, a material, subset parts) and the
    written one hold the same variables."""
    oracle = tmp_path / "oracleT01"
    oracle.write_bytes(oracle_file(cut=False))
    written = tmp_path / "writtenT01"
    WRITERS[writer](written, sample_mesh(1))
    mine = engine.read(written)
    theirs = engine.read(oracle, time_step=1)
    mine.field_data["meshio:time"] = theirs.field_data["meshio:time"]
    _same_fields(mine, theirs)


def test_the_engines_write_the_same_bytes(tmp_path):
    for name, mesh in {
        "plain": sample_mesh(2),
        "factors": sample_mesh(0, factors=True),
        "empty": meshioplusplus.Mesh(np.zeros((0, 3)), []),
    }.items():
        paths = {}
        for engine, write in WRITERS.items():
            paths[engine] = tmp_path / f"{engine}_{name}T01"
            write(paths[engine], mesh)
        assert paths["core"].read_bytes() == paths["python"].read_bytes()


def test_a_real_run_survives_a_rewrite(engine, writer, tmp_path):
    """Every output of OpenRadioss's T01, read, written and read again."""
    for step in range(len(engine.time_values(T01))):
        original = engine.read(T01, time_step=step)
        path = tmp_path / f"copy{step}T01"
        WRITERS[writer](path, original)
        _same_fields(engine.read(path), original)


def test_unnamed_codes_and_unlisted_kinds_round_trip(engine, writer, tmp_path):
    mesh = meshioplusplus.Mesh(
        np.zeros((0, 3)),
        [],
        field_data={
            "radioss_th:global:var40": np.array([1.5]),
            "radioss_th:part:-3:var26": np.array([2.5]),
            "radioss_th:part:-3:PW": np.array([3.5]),
            "radioss_th:type12:9": np.array([[1.0, 2.0]]),
            "radioss_th:type12:9:ids": np.array([-4]),
            "radioss_th:type12:9:variables": np.array([40, 41]),
        },
    )
    path = tmp_path / "otherT01"
    WRITERS[writer](path, mesh)
    fd = engine.read(path).field_data
    assert fd["radioss_th:global:var40"][0] == 1.5
    assert fd["radioss_th:part:-3:var26"][0] == 2.5
    assert fd["radioss_th:part:-3:PW"][0] == 3.5
    np.testing.assert_array_equal(fd["radioss_th:type12:9"], [[1.0, 2.0]])
    np.testing.assert_array_equal(fd["radioss_th:type12:9:ids"], [-4])
    np.testing.assert_array_equal(fd["radioss_th:type12:9:variables"], [40, 41])


def test_ids_and_codes_default_to_one_based_sequences(engine, writer, tmp_path):
    mesh = meshioplusplus.Mesh(
        np.zeros((0, 3)),
        [],
        field_data={"radioss_th:shell:2": np.arange(6, dtype=np.float32).reshape(3, 2)},
    )
    path = tmp_path / "shellT01"
    WRITERS[writer](path, mesh)
    fd = engine.read(path).field_data
    np.testing.assert_array_equal(fd["radioss_th:shell:2"], np.arange(6).reshape(3, 2))
    np.testing.assert_array_equal(fd["radioss_th:shell:2:ids"], [1, 2, 3])
    np.testing.assert_array_equal(fd["radioss_th:shell:2:variables"], [1, 2])


def test_what_has_no_place_in_the_file_is_dropped_with_a_warning(
    writer, tmp_path, capfd
):
    mesh = sample_mesh()
    mesh.field_data["fd"] = np.array([1.5, 2.5])
    mesh.field_data["radioss_th:global:no_such_energy"] = np.array([1.0])
    mesh.field_data["radioss_th:part:5:IE:extra"] = np.array([1.0])
    path = tmp_path / "dropT01"
    WRITERS[writer](path, mesh)
    err = " ".join(capfd.readouterr().err.split())
    assert "not written" in err
    for name in (
        "fd",
        "radioss_th:global:no_such_energy",
        "radioss_th:part:5:IE:extra",
    ):
        assert name in err
    assert path.read_bytes() == expected_bytes([0])


@pytest.mark.parametrize(
    "name, value, message",
    [
        ("radioss_th:node:3", np.arange(3.0), "2-D"),
        ("radioss_th:node:3:ids", np.array([1, 2, 3]), "holds 3 values, expected 2"),
        ("radioss_th:node:3:ids", np.array([1.5, 2.0]), "not an integer"),
        ("radioss_th:node:3:ids", np.array([1, 2**31]), "not an integer"),
        ("radioss_th:node:3:variables", np.array([1, 2]), "holds 2 values, expected 3"),
        ("radioss_th:unit_factors", np.array([1.0, 2.0]), "three values"),
    ],
)
def test_malformed_groups_are_refused(writer, tmp_path, name, value, message):
    mesh = sample_mesh()
    mesh.field_data[name] = value
    with pytest.raises(meshioplusplus.WriteError, match=message):
        WRITERS[writer](tmp_path / "badT01", mesh)


def test_the_public_write_is_found_by_name_and_by_format(tmp_path):
    path = tmp_path / "crashT01"
    meshioplusplus.write(path, sample_mesh(1))
    assert path.read_bytes() == expected_bytes([1])
    other = tmp_path / "history.bin"
    meshioplusplus.write(other, sample_mesh(1), file_format="radioss_th")
    assert other.read_bytes() == path.read_bytes()
    assert meshioplusplus.radioss_th.time_values(path) == [0.5]
    assert "radioss_th" in meshioplusplus.formats()["writable"]
    buf = io.BytesIO()
    meshioplusplus.radioss_th.write(buf, sample_mesh(1))
    assert buf.getvalue() == path.read_bytes()


def _write_series(engine, path, steps):
    w = SERIES[engine](path)
    for time, mesh in steps:
        w.write(time, mesh)
    (w.finalize if engine == "core" else w.close)()


def test_a_series_reads_back_step_by_step(tmp_path):
    out = tmp_path / "runT01"
    times = [0.0, 0.5, 1.25]
    written = meshioplusplus.write_sequence(
        out, [(t, sample_mesh(k)) for k, t in enumerate(times)]
    )
    assert written == [str(out)]
    assert meshioplusplus.radioss_th.time_values(out) == times
    for k, t in enumerate(times):
        mesh = meshioplusplus.read(out, time_step=k)
        assert mesh.field_data["meshio:time"][0] == t
        assert mesh.field_data["radioss_th:part:7:KE"][0] == 2.0 + k
        assert mesh.field_data["radioss_th:node:3"][1, 2] == 5.0 + 100.0 * k
    assert [t for t, _ in meshioplusplus.read_sequence(out)] == times


def test_the_engines_write_the_same_series(tmp_path):
    steps = [(0.5 * k, sample_mesh(k, factors=True)) for k in range(3)]
    for engine in SERIES:
        (tmp_path / engine).mkdir()
        _write_series(engine, tmp_path / engine / "runT01", steps)
        # the descriptions once, then an output per step
        assert (tmp_path / engine / "runT01").read_bytes() == expected_bytes(
            [0, 1, 2], factors=True
        )


@pytest.mark.parametrize("engine", sorted(SERIES))
def test_a_one_step_series_is_a_plain_write(engine, tmp_path):
    _write_series(engine, tmp_path / "seriesT01", [(0.5, sample_mesh(1))])
    assert (tmp_path / "seriesT01").read_bytes() == expected_bytes([1])


@pytest.mark.parametrize("engine", sorted(SERIES))
def test_a_series_refuses_a_step_with_other_variables(engine, tmp_path):
    w = SERIES[engine](tmp_path / "badT01")
    w.write(0.0, sample_mesh(0))
    other = sample_mesh(1)
    del other.field_data["radioss_th:part:7:KE"]
    with pytest.raises(meshioplusplus.WriteError, match="step 1's variables differ"):
        w.write(0.5, other)
    w.write(1.0, sample_mesh(2))
    (w.finalize if engine == "core" else w.close)()
    assert meshioplusplus.radioss_th.time_values(tmp_path / "badT01") == [0.0, 1.0]


@pytest.mark.parametrize("engine", sorted(SERIES))
def test_a_series_takes_the_time_from_the_step(engine, tmp_path):
    _write_series(engine, tmp_path / "timeT01", [(7.0, sample_mesh(0))])
    assert meshioplusplus.radioss_th.time_values(tmp_path / "timeT01") == [7.0]


def test_a_series_goes_to_a_named_format_without_the_file_name(tmp_path):
    out = tmp_path / "history.bin"
    meshioplusplus.write_sequence(
        out, [(0.0, sample_mesh(0)), (1.0, sample_mesh(1))], file_format="radioss_th"
    )
    assert meshioplusplus.radioss_th.time_values(out) == [0.0, 1.0]


def test_a_cut_short_written_file_keeps_its_complete_outputs(engine, tmp_path):
    steps = [(0.5 * k, sample_mesh(k)) for k in range(2)]
    full = tmp_path / "fullT01"
    _write_series("python", full, steps)
    data = full.read_bytes()
    for name, tail in {
        "a record that is not finished": _floats(1.0)[:7],
        "a lone time": _floats(1.0),
        "a time and a record": _floats(1.0) + _floats(1.0, 2.0),
        "a broken closing marker": _floats(1.0)[:-1] + b"\x07",
    }.items():
        path = tmp_path / "cutT01"
        path.write_bytes(data + tail)
        assert engine.time_values(path) == [0.0, 0.5], name
        _same_fields(engine.read(path, time_step=1), engine.read(full, time_step=1))
    # cut inside the second output
    path = tmp_path / "cutT01"
    path.write_bytes(data[: len(data) - 20])
    assert engine.time_values(path) == [0.0]


@pytest.mark.parametrize("width", [4, 8])
@pytest.mark.parametrize("order", ["<", ">"])
def test_the_record_framer_writes_what_it_reads(width, order):
    from meshioplusplus._fortran_records import (
        fortran_record,
        fortran_records,
        sniff_fortran_records,
    )

    payloads = [b"abcd", b"", bytes(range(13))]
    layout = (width, order)
    data = b"".join(fortran_record(p, layout) for p in payloads)
    marker = struct.pack(order + ("i" if width == 4 else "q"), 4)
    assert data.startswith(marker + b"abcd" + marker)
    assert sniff_fortran_records(data) == layout
    spans = fortran_records(data, layout, "test")
    assert [data[o : o + n] for o, n in spans] == payloads
    # a truncated tail is an error unless the reader is lenient
    with pytest.raises(meshioplusplus.ReadError, match="Fortran record"):
        fortran_records(data[:-3], layout, "test")
    with pytest.raises(meshioplusplus.ReadError, match="runs past"):
        fortran_records(data[: -width - 8], layout, "test")
    assert fortran_records(data[:-3], layout, "test", lenient=True) == spans[:2]
    bad = data[:-1] + b"\x7f"
    with pytest.raises(meshioplusplus.ReadError, match="mismatched"):
        fortran_records(bad, layout, "test")
    assert fortran_records(bad, layout, "test", lenient=True) == spans[:2]
