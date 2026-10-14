"""The LS-DYNA binout writer: both engines write the same bytes, which the reader
takes back (and lasso-python, where installed, reads independently)."""

import struct

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.lsdyna_binout import _binout as py_binout
from meshioplusplus.lsdyna_binout import _write as py_write

pytestmark = pytest.mark.filterwarnings("ignore::UserWarning")


@pytest.fixture(params=["core", "python"])
def engine(request):
    return request.param


def write_with(engine, path, mesh):
    if engine == "core":
        _core.lsdyna_binout_write(str(path), mesh)
    else:
        py_write.write(str(path), mesh)


def write_series_with(engine, path, steps):
    if engine == "core":
        writer = _core.LsdynaBinoutSeriesWriter(str(path))
        for time, mesh in steps:
            writer.write(float(time), mesh)
        if writer.num_steps():
            writer.finalize()
    else:
        with py_write.SeriesWriter(str(path)) as writer:
            for time, mesh in steps:
                writer.write(time, mesh)


# -- the hand-built oracle: the layout of the reader's own test builder -------------------


class Lsda:
    """An LSDA file with 8-byte lengths and offsets, little-endian: DATA records,
    then one symbol table of CD and VARIABLE entries."""

    def __init__(self):
        self.data = bytearray(bytes([8, 8, 8, 1, 1, 1, 0, 0]))
        self.data += struct.pack("<QBQ", 17, 7, 0)
        self.entries = []

    def var(self, directory, name, code, values):
        fmt = {3: "i", 9: "f", 10: "d"}[code]
        raw = struct.pack(f"<{len(values)}{fmt}", *values)
        offset = len(self.data)
        self.data += struct.pack(
            "<QBBB", 8 + 3 + len(name) + len(raw), 3, code, len(name)
        )
        self.data += name.encode() + raw
        self.entries.append((directory, name, code, offset, len(values)))

    def finish(self):
        out = bytearray(self.data)
        out[17:25] = struct.pack("<Q", len(out))
        out += struct.pack("<QB", 9, 5)
        cwd = None
        for directory, name, code, offset, count in self.entries:
            if directory != cwd:
                out += struct.pack("<QB", 9 + len(directory), 2) + directory.encode()
                cwd = directory
            out += struct.pack("<QB", 9 + len(name) + 1 + 16, 4) + name.encode()
            out += struct.pack("<BQQ", code, offset, count)
        out += struct.pack("<QBQ", 17, 6, 0)
        return bytes(out)


def oracle_mesh():
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0]]), [("vertex", np.array([[0], [1]]))]
    )
    mesh.point_data["lsdyna:nid"] = np.array([11, 12])
    mesh.point_data["displacement"] = np.array([[0.5, 0, 0], [0.5, 0, 0]])
    mesh.field_data["meshio:time"] = np.array([1.0])
    mesh.field_data["binout:nodout:cycle"] = np.array([100])
    mesh.field_data["binout:glstat:time"] = np.array([0.5])
    mesh.field_data["binout:glstat:kinetic_energy"] = np.array([7.5])
    return mesh


def oracle_bytes():
    f = Lsda()
    f.var("/nodout/metadata", "ids", 3, [11, 12])
    f.var("/nodout/d000001", "time", 9, [1.0])
    f.var("/nodout/d000001", "x_coordinate", 9, [0.0, 1.0])
    f.var("/nodout/d000001", "y_coordinate", 9, [0.0, 0.0])
    f.var("/nodout/d000001", "z_coordinate", 9, [0.0, 0.0])
    f.var("/nodout/d000001", "x_displacement", 9, [0.5, 0.5])
    f.var("/nodout/d000001", "y_displacement", 9, [0.0, 0.0])
    f.var("/nodout/d000001", "z_displacement", 9, [0.0, 0.0])
    f.var("/nodout/d000001", "cycle", 3, [100])
    f.var("/glstat/d000001", "time", 9, [0.5])
    f.var("/glstat/d000001", "kinetic_energy", 9, [7.5])
    return f.finish()


def test_the_writer_follows_the_hand_built_layout(engine, tmp_path):
    path = tmp_path / "binout"
    write_with(engine, path, oracle_mesh())
    assert path.read_bytes() == oracle_bytes()


def test_the_oracle_and_the_writer_read_back_alike(engine, tmp_path):
    path = tmp_path / "binout"
    write_with(engine, path, oracle_mesh())
    oracle = tmp_path / "oracle"
    oracle.write_bytes(oracle_bytes())
    a, b = py_binout.read(path), py_binout.read(oracle)
    np.testing.assert_array_equal(a.points, b.points)
    assert sorted(a.field_data) == sorted(b.field_data)
    assert a.point_data["lsdyna:nid"].tolist() == [11, 12]
    assert a.field_data["binout:glstat:kinetic_energy"][0] == 7.5


# -- a mesh and its steps ----------------------------------------------------------------


def step_mesh(t, glstat=None, extra=True):
    points = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0]]) + [t, 0, 0]
    mesh = meshioplusplus.Mesh(points, [("vertex", np.arange(3)[:, None])])
    mesh.point_data["lsdyna:nid"] = np.array([11, 12, 13])
    mesh.point_data["displacement"] = np.arange(9, dtype=float).reshape(3, 3) * t
    mesh.point_data["rotation"] = np.ones((3, 3)) * t
    if extra:
        mesh.point_data["temperature"] = np.array([1.0, 2.0, 3.0]) * t
    mesh.field_data["meshio:time"] = np.array([t])
    mesh.field_data["binout:nodout:cycle"] = np.array([int(100 * t)])
    if glstat is not None:
        mesh.field_data["binout:glstat:time"] = np.array([glstat])
        mesh.field_data["binout:glstat:kinetic_energy"] = np.array([7.5 * glstat])
        mesh.field_data["binout:elout/beam:time"] = np.array([glstat])
        mesh.field_data["binout:elout/beam:axial"] = np.array([1.0, 2.0]) * glstat
        mesh.field_data["binout:elout/beam:ids"] = np.array([5, 6])
    return mesh


def test_a_mesh_comes_back_by_name(engine, tmp_path):
    path = tmp_path / "binout"
    write_with(engine, path, step_mesh(0.5, glstat=0.25))
    back = meshioplusplus.lsdyna_binout.read(path)
    np.testing.assert_allclose(back.points, step_mesh(0.5).points, rtol=1e-6)
    assert back.point_data["lsdyna:nid"].tolist() == [11, 12, 13]
    np.testing.assert_allclose(
        back.point_data["displacement"], step_mesh(0.5).point_data["displacement"]
    )
    np.testing.assert_allclose(back.point_data["rotation"], np.full((3, 3), 0.5))
    np.testing.assert_allclose(back.point_data["temperature"], [0.5, 1.0, 1.5])
    assert back.field_data["binout:nodout:cycle"].tolist() == [50]
    assert back.field_data["binout:glstat:kinetic_energy"][0] == 1.875
    np.testing.assert_array_equal(
        back.field_data["binout:elout/beam:axial"], [0.25, 0.5]
    )
    assert back.field_data["binout:elout/beam:ids"].tolist() == [5, 6]
    assert back.field_data["meshio:time"][0] == 0.5


def test_the_engines_write_the_same_bytes(tmp_path):
    mesh = step_mesh(0.5, glstat=0.25)
    mesh.field_data["binout:nodout:big"] = np.array([1.0e300, -1.0e300, np.nan])
    core, python = tmp_path / "core", tmp_path / "python"
    _core.lsdyna_binout_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_the_public_write_is_found_by_name_and_format(tmp_path):
    mesh = oracle_mesh()
    meshioplusplus.write(tmp_path / "binout", mesh)
    meshioplusplus.write(tmp_path / "run.bin", mesh, file_format="lsdyna_binout")
    assert (tmp_path / "binout").read_bytes() == oracle_bytes()
    assert (tmp_path / "run.bin").read_bytes() == oracle_bytes()
    assert meshioplusplus.read(tmp_path / "binout").point_data[
        "lsdyna:nid"
    ].tolist() == [11, 12]
    with open(tmp_path / "buffer", "wb") as fh:
        meshioplusplus.lsdyna_binout.write(fh, mesh)
    assert (tmp_path / "buffer").read_bytes() == oracle_bytes()


def test_a_point_free_mesh_is_a_database(engine, tmp_path):
    mesh = meshioplusplus.Mesh(np.zeros((0, 3)), [])
    mesh.field_data["binout:glstat:kinetic_energy"] = np.array([3.5])
    mesh.field_data["meshio:time"] = np.array([2.0])
    path = tmp_path / "binout"
    write_with(engine, path, mesh)
    assert py_binout.time_values(path) == [2.0]
    back = py_binout.read(path)
    assert len(back.points) == 0
    assert back.field_data["binout:glstat:kinetic_energy"][0] == 3.5


def test_what_has_no_place_in_the_file_is_dropped_with_a_warning(
    engine, tmp_path, capfd
):
    mesh = step_mesh(0.5)
    mesh.point_data["pressure_vec"] = np.zeros((3, 2))
    mesh.field_data["note"] = np.array([1.0])
    mesh.cell_data["c"] = [np.zeros(3)]
    path = tmp_path / "binout"
    write_with(engine, path, mesh)
    err = " ".join(capfd.readouterr().err.split())
    assert "not written" in err
    for name in ("pressure_vec", "note", "c"):
        assert name in err
    assert "pressure_vec" not in py_binout.read(path).point_data


# -- series ----------------------------------------------------------------------------


def series_steps():
    # glstat is written every 0.2, nodout every 0.1: a step carries the latest glstat.
    return [
        (0.0, step_mesh(0.0)),
        (0.1, step_mesh(0.1, glstat=0.0)),
        (0.2, step_mesh(0.2, glstat=0.0)),
        (0.3, step_mesh(0.3, glstat=0.25)),
    ]


def test_a_series_keeps_each_databases_own_cadence(engine, tmp_path):
    path = tmp_path / "binout"
    write_series_with(engine, path, series_steps())
    assert py_binout.time_values(path) == pytest.approx([0.0, 0.1, 0.2, 0.3])
    archive = py_binout._Archive(str(path))
    assert sorted(py_binout._steps(archive, ("nodout",))) == [
        f"d00000{k}" for k in (1, 2, 3, 4)
    ]
    assert py_binout._steps(archive, ("glstat",)) == ["d000001", "d000002"]
    last = py_binout.read(path, time_step=-1)
    assert last.field_data["binout:glstat:kinetic_energy"][0] == 1.875
    assert last.field_data["binout:glstat:time"][0] == 0.25
    first = py_binout.read(path, time_step=1)
    assert first.field_data["binout:glstat:kinetic_energy"][0] == 0.0


def test_the_engines_write_the_same_series(tmp_path):
    write_series_with("core", tmp_path / "core", series_steps())
    write_series_with("python", tmp_path / "python", series_steps())
    assert (tmp_path / "core").read_bytes() == (tmp_path / "python").read_bytes()


def test_write_sequence_uses_the_series_writer(tmp_path):
    path = tmp_path / "binout"
    meshioplusplus.write_sequence(path, series_steps(), file_format="lsdyna_binout")
    assert meshioplusplus.read_metadata(path)["time_values"] == pytest.approx(
        [0.0, 0.1, 0.2, 0.3]
    )
    assert len([t for t, _ in meshioplusplus.read_sequence(path)]) == 4


def test_a_rewrite_is_stable(tmp_path):
    first, second, third = (tmp_path / n for n in ("binout", "binout_b", "binout_c"))
    meshioplusplus.write_sequence(first, series_steps(), file_format="lsdyna_binout")
    meshioplusplus.write_sequence(
        second,
        meshioplusplus.read_sequence(first, file_format="lsdyna_binout"),
        file_format="lsdyna_binout",
    )
    meshioplusplus.write_sequence(
        third,
        meshioplusplus.read_sequence(second, file_format="lsdyna_binout"),
        file_format="lsdyna_binout",
    )
    assert second.read_bytes() == third.read_bytes()


def test_a_one_step_series_is_a_plain_write(engine, tmp_path):
    mesh = step_mesh(0.5, glstat=0.25)
    mesh.field_data["meshio:time"] = np.array([0.5])
    plain, series = tmp_path / "plain", tmp_path / "series"
    write_with(engine, plain, mesh)
    write_series_with(engine, series, [(0.5, mesh)])
    assert plain.read_bytes() == series.read_bytes()


@pytest.mark.parametrize(
    "make, message",
    [
        (lambda m: m.point_data.pop("temperature"), "differ from the first"),
        (
            lambda m: m.point_data.__setitem__("lsdyna:nid", np.array([1, 2, 3])),
            "differ from the first",
        ),
        (
            lambda m: m.field_data.__setitem__("binout:glstat:extra", np.array([1.0])),
            "of 'glstat' differ",
        ),
    ],
)
def test_a_series_refuses_a_step_with_other_nodes_or_variables(
    engine, tmp_path, make, message
):
    second = step_mesh(0.2, glstat=0.0)
    make(second)
    with pytest.raises(Exception, match=message):
        write_series_with(
            engine, tmp_path / "bad", [(0.1, step_mesh(0.1, glstat=0.0)), (0.2, second)]
        )


def test_a_series_refuses_another_point_count(engine, tmp_path):
    other = meshioplusplus.Mesh(np.zeros((2, 3)), [("vertex", np.arange(2)[:, None])])
    with pytest.raises(Exception, match="has 2 points, the first 3"):
        write_series_with(
            engine, tmp_path / "bad", [(0.0, step_mesh(0.0)), (0.1, other)]
        )


def test_a_database_may_appear_in_a_later_step(engine, tmp_path):
    path = tmp_path / "binout"
    write_series_with(
        engine, path, [(0.0, step_mesh(0.0)), (0.1, step_mesh(0.1, glstat=0.1))]
    )
    last = py_binout.read(path, time_step=-1)
    assert last.field_data["binout:glstat:kinetic_energy"][0] == pytest.approx(0.75)


# -- refusals --------------------------------------------------------------------------


def test_an_id_outside_32_bits_is_refused(engine, tmp_path):
    mesh = step_mesh(0.5)
    mesh.point_data["lsdyna:nid"] = np.array([1, 2, 2**40])
    with pytest.raises(Exception, match="32-bit"):
        write_with(engine, tmp_path / "x", mesh)


@pytest.mark.parametrize(
    "name, message",
    [
        ("binout:d000001:x", "not a name a database can have"),
        ("binout:metadata:x", "not a name a database can have"),
        ("binout:a//b:x", "not a name a binout can hold"),
    ],
)
def test_a_database_name_a_binout_cannot_hold_is_refused(
    engine, tmp_path, name, message
):
    mesh = step_mesh(0.5)
    mesh.field_data[name] = np.array([1.0])
    with pytest.raises(Exception, match=message):
        write_with(engine, tmp_path / "x", mesh)


def test_a_mesh_with_nothing_to_write_is_refused(engine, tmp_path):
    with pytest.raises(Exception, match="no database"):
        write_with(engine, tmp_path / "x", meshioplusplus.Mesh(np.zeros((0, 3)), []))


def test_a_variable_that_repeats_one_of_nodout_is_refused(engine, tmp_path):
    mesh = step_mesh(0.5)
    mesh.field_data["binout:nodout:temperature"] = np.array([1.0])
    with pytest.raises(Exception, match="repeats a variable of nodout|of its own"):
        write_with(engine, tmp_path / "x", mesh)


# -- an independent reader --------------------------------------------------------------


def test_lasso_reads_the_written_file(tmp_path):
    lasso = pytest.importorskip("lasso.dyna", exc_type=ImportError)
    path = tmp_path / "binout"
    write_series_with("core", path, series_steps())
    binout = lasso.Binout(str(path))
    assert set(binout.read()) >= {"nodout", "glstat", "elout"}
    np.testing.assert_array_equal(binout.read("nodout", "ids"), [11, 12, 13])
    np.testing.assert_allclose(
        binout.read("nodout", "time"), [0.0, 0.1, 0.2, 0.3], rtol=1e-6
    )
    x = binout.read("nodout", "x_coordinate")
    assert x.shape == (4, 3)
    np.testing.assert_allclose(x[3], [0.3, 1.3, 1.3], rtol=1e-6)
    np.testing.assert_allclose(binout.read("glstat", "kinetic_energy"), [0.0, 1.875])
    np.testing.assert_allclose(
        binout.read("elout", "beam", "axial"), [[0.0, 0.0], [0.25, 0.5]]
    )
