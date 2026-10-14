"""The Abaqus results-file writer: both engines write the same bytes, which a
record-by-record oracle spells out, and which the reader (and pybaqus, for the
ASCII layout, where installed) takes back."""

import pathlib
import sys

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.abaqus_fil import _abaqus_fil as py_fil
from meshioplusplus.abaqus_fil import _write as py_write

TOOLS = pathlib.Path(__file__).resolve().parents[2] / "tools"
MESHES = pathlib.Path(__file__).parent / "meshes" / "abaqus_fil"

pytestmark = pytest.mark.filterwarnings("ignore::UserWarning")


@pytest.fixture(params=["core", "python"])
def engine(request):
    return request.param


@pytest.fixture(params=[False, True], ids=["ascii", "binary"])
def binary(request):
    return request.param


def write_with(engine, path, mesh, binary=False):
    if engine == "core":
        _core.abaqus_fil_write(str(path), mesh, binary)
    else:
        py_write.write(str(path), mesh, binary=binary)


def write_series_with(engine, path, steps, binary=False):
    if engine == "core":
        writer = _core.AbaqusFilSeriesWriter(str(path), binary)
        for time, mesh in steps:
            writer.write(float(time), mesh)
        if writer.num_steps():
            writer.finalize()
    else:
        with py_write.SeriesWriter(str(path), binary=binary) as writer:
            for time, mesh in steps:
                writer.write(time, mesh)


# -- the oracle: one C3D8R, records spelled out with the fixture generator's encoders ----


def generator():
    sys.path.insert(0, str(TOOLS))
    try:
        import gen_abaqus_fil_fixtures as gen
    finally:
        sys.path.remove(str(TOOLS))
    return gen


UNIT = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ],
    dtype=float,
)


def oracle_mesh():
    mesh = meshioplusplus.Mesh(UNIT, [("hexahedron", np.arange(8)[None, :])])
    mesh.point_data["U"] = np.arange(24, dtype=float).reshape(8, 3) / 8
    mesh.cell_data["S"] = [np.arange(6, dtype=float)[None, :]]
    mesh.field_data["meshio:time"] = np.array([0.5])
    return mesh


def oracle_records():
    I = lambda v: ("I", v)  # noqa: E731,E741
    D = lambda v: ("D", float(v))  # noqa: E731
    A = lambda v: ("A", v.ljust(8)[:8])  # noqa: E731
    out = []

    def rec(key, *items):
        out.append([I(key)] + list(items))

    rec(1921, A("meshio++"), A(""), A(""), A(""), I(1), I(8), D(1.0))
    title = "meshio++ results file".ljust(80)
    rec(1922, *[A(title[8 * k : 8 * k + 8]) for k in range(10)])
    rec(1902, *[I(v) for v in [1, 2, 3] + [0] * 32])
    rec(1900, I(1), A("C3D8R"), *[I(n) for n in range(1, 9)])
    for n, xyz in enumerate(UNIT, start=1):
        rec(1901, I(n), *[D(x) for x in xyz])
    blank = "".ljust(80)
    rec(
        2000,
        D(0.5),
        D(0.5),
        D(0.0),
        D(0.0),
        I(1),
        I(1),
        I(1),
        I(0),
        D(0.0),
        D(0.0),
        D(1.0),
        *[A(blank[8 * k : 8 * k + 8]) for k in range(10)],
    )
    rec(1911, I(0), A(""), A(""))
    rec(1, I(1), I(0), I(0), I(1), A(""), I(3), I(3), I(0), I(0))
    rec(11, *[D(v) for v in range(6)])
    rec(1911, I(1), A(""))
    for n in range(1, 9):
        rec(101, I(n), *[D((n - 1) * 3 / 8 + c / 8) for c in range(3)])
    rec(2001)
    return out


def test_the_writer_follows_the_record_oracle(engine, binary, tmp_path):
    gen = generator()
    path = tmp_path / "oracle.fil"
    write_with(engine, path, oracle_mesh(), binary)
    recs = oracle_records()
    expected = gen.binary_bytes(recs, "<") if binary else gen.ascii_bytes(recs)
    assert path.read_bytes() == expected


def test_the_oracle_and_the_writer_read_back_alike(engine, binary, tmp_path):
    gen = generator()
    path = tmp_path / "written.fil"
    write_with(engine, path, oracle_mesh(), binary)
    oracle = tmp_path / "oracle.fil"
    recs = oracle_records()
    oracle.write_bytes(gen.binary_bytes(recs, "<") if binary else gen.ascii_bytes(recs))
    a, b = py_fil.read(path), py_fil.read(oracle)
    np.testing.assert_array_equal(a.points, b.points)
    np.testing.assert_array_equal(a.point_data["U"], b.point_data["U"])
    np.testing.assert_array_equal(a.cell_data["S"][0], b.cell_data["S"][0])
    assert a.field_data["meshio:time"].ravel()[0] == 0.5


# -- a real-shaped model, round tripped -------------------------------------------------


def same_mesh(a, b):
    np.testing.assert_allclose(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    assert sorted(a.point_data) == sorted(b.point_data)
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.point_data:
        np.testing.assert_allclose(
            a.point_data[name], b.point_data[name], rtol=1e-12, equal_nan=True
        )
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_allclose(x, y, rtol=1e-12, equal_nan=True)
    assert sorted((r.name, r.kind) for r in a.regions) == sorted(
        (r.name, r.kind) for r in b.regions
    )
    for name in a.field_data:
        np.testing.assert_allclose(
            a.field_data[name], b.field_data[name], equal_nan=True
        )


def test_a_real_shaped_model_survives_a_round_trip(engine, binary, tmp_path):
    source = MESHES / "model.fil"
    times = py_fil.time_values(source)
    steps = [(t, py_fil.read(source, time_step=i)) for i, t in enumerate(times)]
    path = tmp_path / "again.fil"
    write_series_with(engine, path, steps, binary)
    assert py_fil.time_values(path) == times
    for i, (_, mesh) in enumerate(steps):
        same_mesh(mesh, py_fil.read(path, time_step=i))
        np.testing.assert_array_equal(
            mesh.cell_data["abaqus:id"][0],
            py_fil.read(path, time_step=i).cell_data["abaqus:id"][0],
        )


def test_the_engines_write_the_same_bytes(binary, tmp_path):
    source = MESHES / "model.fil"
    steps = [
        (t, py_fil.read(source, time_step=i))
        for i, t in enumerate(py_fil.time_values(source))
    ]
    core, python = tmp_path / "core.fil", tmp_path / "python.fil"
    write_series_with("core", core, steps, binary)
    write_series_with("python", python, steps, binary)
    assert core.read_bytes() == python.read_bytes()
    one = tmp_path / "one_core.fil"
    other = tmp_path / "one_python.fil"
    write_with("core", one, steps[0][1], binary)
    write_with("python", other, steps[0][1], binary)
    assert one.read_bytes() == other.read_bytes()


def test_a_cell_block_with_no_abaqus_element_is_dropped_with_a_warning(
    engine, tmp_path, capfd
):
    mesh = oracle_mesh()
    mesh.cells.append(meshioplusplus.CellBlock("vertex", np.array([[0]])))
    mesh.cell_data["S"].append(np.zeros((1, 6)))
    path = tmp_path / "vertex.fil"
    write_with(engine, path, mesh)
    assert "vertex" in " ".join(capfd.readouterr().err.split())
    assert [c.type for c in py_fil.read(path).cells] == ["hexahedron"]


def test_the_public_write_is_found_by_extension_and_encoding(tmp_path):
    mesh = oracle_mesh()
    meshioplusplus.write(tmp_path / "a.fil", mesh)
    meshioplusplus.write(tmp_path / "b.fil", mesh, binary=True)
    meshioplusplus.write(tmp_path / "c.bin", mesh, file_format="abaqus_fil")
    assert (tmp_path / "a.fil").read_bytes() == (tmp_path / "c.bin").read_bytes()
    assert (tmp_path / "a.fil").read_bytes()[:1] == b"*"
    assert (tmp_path / "b.fil").read_bytes()[:4] != b"*I 1"
    for name in ("a.fil", "b.fil"):
        np.testing.assert_array_equal(
            meshioplusplus.read(tmp_path / name).point_data["U"], mesh.point_data["U"]
        )
    with open(tmp_path / "buffer", "wb") as fh:
        meshioplusplus.abaqus_fil.write(fh, mesh)
    assert (tmp_path / "buffer").read_bytes() == (tmp_path / "a.fil").read_bytes()


# -- every cell type, every location ------------------------------------------------------

TYPES = {
    "line": 2,
    "line3": 3,
    "triangle": 3,
    "quad": 4,
    "triangle6": 6,
    "quad8": 8,
    "quad9": 9,
    "tetra": 4,
    "tetra10": 10,
    "pyramid": 5,
    "pyramid13": 13,
    "wedge": 6,
    "wedge15": 15,
    "hexahedron": 8,
    "hexahedron20": 20,
}


def every_type():
    rng = np.random.default_rng(5)
    mesh = meshioplusplus.Mesh(
        rng.random((30, 3)),
        [(t, (np.arange(k)[None, :] * 2 + 1) % 30) for t, k in TYPES.items()],
    )
    return mesh


def test_every_cell_type_survives_a_round_trip(engine, binary, tmp_path):
    mesh = every_type()
    path = tmp_path / "types.fil"
    write_with(engine, path, mesh, binary)
    back = py_fil.read(path)
    assert [c.type for c in back.cells] == [c.type for c in mesh.cells]
    for a, b in zip(mesh.cells, back.cells):
        np.testing.assert_array_equal(a.data, b.data)
    np.testing.assert_allclose(back.points, mesh.points)


def results_mesh():
    points = np.array(
        [[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1]]
    )
    mesh = meshioplusplus.Mesh(
        points,
        [("hexahedron", np.arange(8)[None, :] % 6), ("quad", np.array([[0, 1, 2, 3]]))],
    )
    mesh.cells[0] = meshioplusplus.CellBlock("tetra", np.array([[0, 1, 3, 4]]))
    mesh.point_data["U"] = np.arange(18, dtype=float).reshape(6, 3) / 4
    mesh.point_data["abaqus:id"] = np.array([11, 12, 13, 14, 15, 16])
    mesh.point_data["RF"] = np.full((6, 3), np.nan)
    mesh.point_data["RF"][2] = [1.0, 2.0, 3.0]
    mesh.point_data["PEEQ"] = np.array([1.0, 2.0, np.nan, 4.0, 5.0, 6.0])
    mesh.cell_data["abaqus:id"] = [np.array([101]), np.array([102])]
    shell = np.full((1, 12), np.nan)
    shell[0, :3] = [1.0, 2.0, 3.0]
    mesh.cell_data["S"] = [np.arange(12, dtype=float).reshape(1, 12), shell]
    mesh.field_data["abaqus:layout:S"] = np.array([2, 6])
    mesh.cell_data["S@sp5"] = [np.full((1, 3), np.nan), np.array([[7.0, 8.0, 9.0]])]
    mesh.cell_data["MISES"] = [np.array([10.0]), np.array([np.nan])]
    mesh.cell_data["NFORC"] = [
        np.arange(12, dtype=float).reshape(1, 12),
        np.full((1, 12), np.nan),
    ]
    mesh.field_data["abaqus:layout:NFORC"] = np.array([4, 3])
    return mesh


def test_results_come_back_by_node_and_element(engine, binary, tmp_path):
    mesh = results_mesh()
    path = tmp_path / "results.fil"
    write_with(engine, path, mesh, binary)
    back = py_fil.read(path)
    assert back.point_data["abaqus:id"].tolist() == [11, 12, 13, 14, 15, 16]
    assert back.cell_data["abaqus:id"][0].tolist() == [101]
    np.testing.assert_allclose(back.point_data["U"], mesh.point_data["U"])
    np.testing.assert_allclose(back.point_data["RF"][2], [1, 2, 3])
    assert np.isnan(back.point_data["RF"][0]).all()
    # averaged at the nodes: the NaN point has no value
    np.testing.assert_allclose(
        back.point_data["PEEQ"], mesh.point_data["PEEQ"], equal_nan=True
    )
    # two integration points of six components, point-major
    np.testing.assert_allclose(back.cell_data["S"][0], np.arange(12.0)[None, :])
    assert back.field_data["abaqus:layout:S"].tolist() == [2, 6]
    # the shell's three components at section point 5
    np.testing.assert_allclose(back.cell_data["S@sp5"][1][0, :3], [7.0, 8.0, 9.0])
    # the centroid
    assert back.cell_data["MISES"][0][0] == 10.0
    assert np.isnan(back.cell_data["MISES"][1][0])
    # element nodes: NFORC's four nodes of the tetra, in the cell's node order
    np.testing.assert_allclose(back.cell_data["NFORC"][0], np.arange(12.0)[None, :])
    assert back.field_data["abaqus:layout:NFORC"].tolist() == [4, 3]


def test_rebar_and_numbered_keys_round_trip(engine, tmp_path):
    mesh = meshioplusplus.Mesh(UNIT, [("hexahedron", np.arange(8)[None, :])])
    mesh.cell_data["S@rebar:RB1"] = [np.array([[1.0, 2.0]])]
    mesh.field_data["abaqus:layout:S@rebar:RB1"] = np.array([2, 1])
    mesh.cell_data["key_900"] = [np.array([4.5])]
    mesh.point_data["key_999"] = np.arange(8, dtype=float)
    path = tmp_path / "rebar.fil"
    write_with(engine, path, mesh)
    back = py_fil.read(path)
    np.testing.assert_allclose(back.cell_data["S@rebar:RB1"][0], [[1.0, 2.0]])
    assert back.cell_data["key_900"][0][0] == 4.5
    np.testing.assert_allclose(back.point_data["key_999"], np.arange(8.0))


def test_a_mode_is_a_1980_record(engine, binary, tmp_path):
    mesh = meshioplusplus.Mesh(UNIT, [("hexahedron", np.arange(8)[None, :])])
    mesh.point_data["U"] = np.ones((8, 3))
    mesh.field_data["meshio:time"] = np.array([0.0])
    mesh.field_data["abaqus:mode"] = np.array([2])
    mesh.field_data["abaqus:eigenvalue"] = np.array([123.5])
    mesh.field_data["abaqus:generalized_mass"] = np.array([2.0])
    mesh.field_data["abaqus:participation_factor"] = np.array([0.1, 0.2])
    mesh.field_data["abaqus:effective_mass"] = np.array([0.01, 0.04])
    path = tmp_path / "mode.fil"
    write_with(engine, path, mesh, binary)
    back = py_fil.read(path)
    assert back.field_data["abaqus:mode"].ravel()[0] == 2
    assert back.field_data["abaqus:eigenvalue"].ravel()[0] == 123.5
    np.testing.assert_allclose(
        back.field_data["abaqus:participation_factor"], [0.1, 0.2]
    )
    np.testing.assert_allclose(back.field_data["abaqus:effective_mass"], [0.01, 0.04])


# -- sets ---------------------------------------------------------------------------------


def test_sets_come_back_as_regions_with_long_names_numbered(engine, binary, tmp_path):
    mesh = results_mesh()
    long_name = "ASSEMBLY_PART-1-1_CLAMPED_NODES"
    mesh.regions = [
        meshioplusplus.Region("FIXED", "point", np.array([0, 3])),
        meshioplusplus.Region(long_name, "point", np.arange(6)),
        meshioplusplus.Region("1", "cell", np.array([0, 1])),
        meshioplusplus.Region("EMPTY", "cell", np.zeros(0, dtype=np.int64)),
        meshioplusplus.Region("wall", "side", np.array([[0, 1]])),
    ]
    path = tmp_path / "sets.fil"
    write_with(engine, path, mesh, binary)
    back = py_fil.read(path)
    by_name = {(r.name, r.kind): r.entries.tolist() for r in back.regions}
    assert by_name[("FIXED", "point")] == [0, 3]
    assert by_name[(long_name, "point")] == [0, 1, 2, 3, 4, 5]
    assert by_name[("1", "cell")] == [0, 1]
    assert by_name[("EMPTY", "cell")] == []
    assert not any(r.kind == "side" for r in back.regions)


def test_what_has_no_place_in_the_file_is_dropped_with_a_warning(
    engine, tmp_path, capfd
):
    mesh = results_mesh()
    mesh.point_data["temperature_like"] = np.zeros(6)
    mesh.cell_data["unknown"] = [np.zeros(1), np.zeros(1)]
    mesh.field_data["note"] = np.array([1.0])
    mesh.regions = [meshioplusplus.Region("wall", "side", np.array([[0, 1]]))]
    path = tmp_path / "dropped.fil"
    write_with(engine, path, mesh)
    err = " ".join(capfd.readouterr().err.split())
    assert "not written" in err
    for name in ("temperature_like", "unknown", "note"):
        assert name in err
    assert "wall" in err


# -- NaN ----------------------------------------------------------------------------------


def test_nan_rows_are_left_out_and_trailing_nan_is_trimmed(engine, tmp_path):
    mesh = meshioplusplus.Mesh(UNIT, [("hexahedron", np.arange(8)[None, :])])
    u = np.arange(24, dtype=float).reshape(8, 3)
    u[1] = np.nan
    u[2, 2] = np.nan
    mesh.point_data["U"] = u
    path = tmp_path / "nan.fil"
    write_with(engine, path, mesh)
    back = py_fil.read(path).point_data["U"]
    assert np.isnan(back[1]).all()
    np.testing.assert_array_equal(back[2], [6.0, 7.0, np.nan])
    np.testing.assert_array_equal(back[3], u[3])


# -- series -------------------------------------------------------------------------------


def step_mesh(t, procedure=1):
    mesh = oracle_mesh()
    mesh.point_data["U"] = mesh.point_data["U"] * t
    mesh.cell_data["S"] = [mesh.cell_data["S"][0] * t]
    mesh.field_data["abaqus:step"] = np.array([1 if t < 2 else 2])
    mesh.field_data["abaqus:procedure"] = np.array([procedure])
    return mesh


def test_a_series_reads_back_increment_by_increment(engine, binary, tmp_path):
    path = tmp_path / "series.fil"
    write_series_with(
        engine, path, [(t, step_mesh(t)) for t in (0.5, 1.0, 2.5)], binary
    )
    assert py_fil.time_values(path) == [0.5, 1.0, 2.5]
    for i, t in enumerate((0.5, 1.0, 2.5)):
        back = py_fil.read(path, time_step=i)
        np.testing.assert_allclose(back.point_data["U"], step_mesh(t).point_data["U"])
        assert back.field_data["meshio:time"].ravel()[0] == t
        assert back.field_data["abaqus:step"].ravel()[0] == (1 if t < 2 else 2)
        assert back.field_data["abaqus:increment"].ravel()[0] == i + 1


def test_write_sequence_and_the_encoding_option(tmp_path):
    steps = [(t, step_mesh(t)) for t in (0.5, 1.0)]
    ascii_path, binary_path = tmp_path / "a.fil", tmp_path / "b.fil"
    meshioplusplus.write_sequence(ascii_path, steps)
    meshioplusplus.write_sequence(binary_path, steps, binary=True)
    assert (
        py_fil.time_values(ascii_path) == py_fil.time_values(binary_path) == [0.5, 1.0]
    )
    assert ascii_path.read_bytes()[:1] == b"*" and binary_path.read_bytes()[:1] != b"*"
    assert len(list(meshioplusplus.read_sequence(ascii_path))) == 2


def test_a_rewrite_is_stable(binary, tmp_path):
    first, second, third = (tmp_path / n for n in ("a.fil", "b.fil", "c.fil"))
    meshioplusplus.write_sequence(
        first, [(t, step_mesh(t)) for t in (0.5, 1.0)], binary=binary
    )
    meshioplusplus.write_sequence(
        second, meshioplusplus.read_sequence(first), binary=binary
    )
    meshioplusplus.write_sequence(
        third, meshioplusplus.read_sequence(second), binary=binary
    )
    assert second.read_bytes() == third.read_bytes()


def test_a_one_step_series_is_a_plain_write(engine, binary, tmp_path):
    mesh = oracle_mesh()
    plain, series = tmp_path / "plain.fil", tmp_path / "series.fil"
    write_with(engine, plain, mesh, binary)
    write_series_with(engine, series, [(0.5, mesh)], binary)
    assert plain.read_bytes() == series.read_bytes()


def test_a_series_refuses_other_points_or_cells(engine, tmp_path):
    other = meshioplusplus.Mesh(UNIT[:4], [("quad", np.arange(4)[None, :])])
    with pytest.raises(Exception, match="differ from the first"):
        write_series_with(
            engine, tmp_path / "bad.fil", [(0.5, oracle_mesh()), (1.0, other)]
        )


def test_a_mesh_with_no_results_is_just_the_model(engine, tmp_path):
    mesh = every_type()
    path = tmp_path / "model.fil"
    write_with(engine, path, mesh)
    assert py_fil.time_values(path) == []
    assert len(py_fil.read(path).cells) == len(TYPES)


# -- refusals -----------------------------------------------------------------------------


def test_labels_must_be_positive_and_unique(engine, tmp_path):
    mesh = oracle_mesh()
    mesh.point_data["abaqus:id"] = np.array([1, 2, 3, 4, 5, 6, 7, 7])
    with pytest.raises(Exception, match="not unique"):
        write_with(engine, tmp_path / "x.fil", mesh)
    mesh.point_data["abaqus:id"] = np.array([0, 1, 2, 3, 4, 5, 6, 7])
    with pytest.raises(Exception, match="positive"):
        write_with(engine, tmp_path / "x.fil", mesh)
    mesh.point_data["abaqus:id"] = np.array([1, 2, 3])
    with pytest.raises(Exception, match="one label per node"):
        write_with(engine, tmp_path / "x.fil", mesh)


def test_a_cell_or_set_naming_a_missing_point_is_refused(engine, tmp_path):
    mesh = oracle_mesh()
    mesh.cells[0].data[0, 0] = 99
    with pytest.raises(Exception, match="point the mesh does not have"):
        write_with(engine, tmp_path / "x.fil", mesh)
    mesh = oracle_mesh()
    mesh.regions = [meshioplusplus.Region("BAD", "point", np.array([99]))]
    with pytest.raises(Exception, match="point the mesh does not have"):
        write_with(engine, tmp_path / "x.fil", mesh)


def test_a_layout_that_does_not_match_is_refused(engine, tmp_path):
    mesh = oracle_mesh()
    mesh.field_data["abaqus:layout:S"] = np.array([2, 6])
    with pytest.raises(Exception, match="does not match"):
        write_with(engine, tmp_path / "x.fil", mesh)


def test_a_value_outside_what_an_item_holds_is_refused(engine, tmp_path):
    mesh = oracle_mesh()
    mesh.point_data["U"][0, 0] = 1e150
    with pytest.raises(Exception, match="range"):
        write_with(engine, tmp_path / "x.fil", mesh)
    mesh = oracle_mesh()
    mesh.point_data["U"][0, 0] = np.inf
    with pytest.raises(Exception, match="not finite"):
        write_with(engine, tmp_path / "x.fil", mesh)


# -- an independent reader --------------------------------------------------------------


def test_pybaqus_reads_the_written_file(tmp_path):
    pybaqus = pytest.importorskip("pybaqus", exc_type=ImportError)
    mesh = meshioplusplus.Mesh(
        np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [1, 1, 0],
                [0, 1, 0],
                [0, 0, 1],
                [1, 0, 1],
                [1, 1, 1],
                [0, 1, 1],
            ],
            dtype=float,
        ),
        [("hexahedron", np.arange(8)[None, :])],
    )
    mesh.point_data["abaqus:id"] = np.arange(101, 109)
    mesh.cell_data["abaqus:id"] = [np.array([5])]
    mesh.point_data["U"] = np.arange(24, dtype=float).reshape(8, 3) / 8
    mesh.cell_data["S"] = [np.arange(6, dtype=float)[None, :]]
    mesh.regions = [meshioplusplus.Region("FIXED", "point", np.array([0, 1]))]
    steps = [(0.5, mesh), (1.0, mesh)]
    path = tmp_path / "pybaqus.fil"
    meshioplusplus.write_sequence(path, steps, file_format="abaqus_fil")
    model = pybaqus.open_fil(str(path))
    assert sorted(model.nodes) == list(range(101, 109))
    assert list(model.elements) == [5]
    assert "FIXED" in model.node_sets
    np.testing.assert_allclose(
        model.get_nodal_vector_result("U", 1, 2), mesh.point_data["U"]
    )
    stress = np.column_stack(
        [np.ravel(model.get_element_result(f"S{c}", 1, 2)) for c in range(1, 7)]
    )
    np.testing.assert_allclose(stress, np.arange(6.0)[None, :])
