"""The Nastran OP2 writer: both engines write the same bytes, which a hand-framed
oracle spells out, and which the reader (and pyNastran, where installed) takes back."""

import struct

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._exceptions import WriteError
from meshioplusplus.nastran_op2 import _op2 as py_op2
from meshioplusplus.nastran_op2 import _write as py_write

pytestmark = pytest.mark.filterwarnings("ignore::UserWarning")


@pytest.fixture(params=["core", "python"])
def engine(request):
    return request.param


def write_with(engine, path, mesh):
    """Write through one engine's own entry point (never the other's)."""
    if engine == "core":
        _core.nastran_op2_write(str(path), mesh)
    else:
        py_write.write(str(path), mesh)


def write_series_with(engine, path, steps):
    if engine == "core":
        writer = _core.NastranOp2SeriesWriter(str(path))
        for time, mesh in steps:
            writer.write(float(time), mesh)
        if writer.num_steps():
            writer.finalize()
    else:
        with py_write.SeriesWriter(str(path)) as writer:
            for time, mesh in steps:
                writer.write(time, mesh)


# -- the hand-framed oracle: one quad, built word by word ----------------------------


def _marker(value):
    return struct.pack("<i", 4) + struct.pack("<i", value) + struct.pack("<i", 4)


def _block(payload):
    n = struct.pack("<i", len(payload))
    return n + payload + n


def _record(payload):
    return _marker(len(payload) // 4) + _block(payload)


def _i(*values):
    return struct.pack(f"<{len(values)}i", *values)


def _f(*values):
    return struct.pack(f"<{len(values)}f", *values)


def _name(text):
    return text.ljust(8).encode()


def _table(name, records):
    out = _record(_name(name)) + _marker(-1) + _record(_i(101, 0, 0, 0, 0, 0, 0))
    k = -2
    out += _marker(k) + _marker(1) + _marker(0) + _record(_name(name))
    for payload in records:
        k -= 1
        out += _marker(k) + _marker(1) + _marker(0) + _record(payload)
    k -= 1
    return out + _marker(k) + _marker(1) + _marker(0) + _marker(0)


def _result_header(tcode, etype, num_wide, s_code=0, thermal=0):
    words = [0x20202020] * 146
    for i in range(50):
        words[i] = 0
    words[0], words[1], words[2], words[3] = 11, tcode, etype, 1
    words[4], words[5] = 1, 0
    words[8], words[9], words[10], words[22] = 1, num_wide, s_code, thermal
    return struct.pack("<146I", *words)


END = _i(65535, 65535, 65535)


def oracle_bytes():
    """The file of ``oracle_mesh()``, framed by hand."""
    out = _marker(3) + _block(_i(1, 1, 0)) + _marker(7)
    out += _block(b"NASTRAN FORT TAPE ID CODE - ") + _record(_name("XXXXXXXX"))
    out += _marker(-1) + _marker(0)
    grids = b"".join(
        _i(n + 1, 0) + _f(*xyz) + _i(0, 0, 0)
        for n, xyz in enumerate([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
    )
    out += _table("GEOM1", [_i(4501, 45, 1) + grids, END])
    quad = _i(10, 7, 1, 2, 3, 4) + _i(0) * 8
    out += _table("GEOM2", [_i(2958, 51, 177) + quad, END])
    pshell = _i(7, 1) + _f(1.0) + _i(0) * 8
    out += _table("EPT", [_i(2302, 23, 283) + pshell, END])
    mat1 = _i(1) + _f(1.0, 0.0, 0.3) + _i(0) * 8
    out += _table("MPT", [_i(103, 1, 77) + mat1, END])
    disp = b"".join(
        _i(10 * n + 1, 1) + _f(*(n + 0.5 * c for c in range(6))) for n in range(1, 5)
    )
    out += _table("OUGV1", [_result_header(1, 0, 8), disp])
    stress = _i(10 * 10 + 1) + _f(*(100.0 + k for k in range(16)))
    out += _table("OES1X1", [_result_header(5, 33, 17, 1), stress])
    return out + _marker(0)


def oracle_mesh():
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    mesh = meshioplusplus.Mesh(points, [("quad", np.array([[0, 1, 2, 3]]))])
    steps = np.arange(1, 5)[:, None] + 0.5 * np.arange(6)[None, :]
    mesh.point_data["DISPLACEMENT"] = steps[:, :3].copy()
    mesh.point_data["DISPLACEMENT_ROT"] = steps[:, 3:].copy()
    mesh.cell_data["nastran:eid"] = [np.array([10])]
    mesh.cell_data["nastran:pid"] = [np.array([7])]
    names = ["FD1", "X1", "Y1", "TXY1", "ANGLE1", "MAJOR1", "MINOR1", "VON_MISES1"]
    names += ["FD2", "X2", "Y2", "TXY2", "ANGLE2", "MAJOR2", "MINOR2", "VON_MISES2"]
    for k, name in enumerate(names):
        mesh.cell_data["STRESS:" + name] = [np.array([100.0 + k])]
    return mesh


def test_the_writer_follows_the_hand_framed_layout(engine, tmp_path):
    path = tmp_path / "quad.op2"
    write_with(engine, path, oracle_mesh())
    assert path.read_bytes() == oracle_bytes()


def test_the_oracle_and_the_writer_read_back_alike(engine, tmp_path):
    path = tmp_path / "quad.op2"
    write_with(engine, path, oracle_mesh())
    oracle = tmp_path / "oracle.op2"
    oracle.write_bytes(oracle_bytes())
    a, b = py_op2.read(path), py_op2.read(oracle)
    np.testing.assert_array_equal(a.points, b.points)
    assert sorted(a.cell_data) == sorted(b.cell_data)
    assert a.cells[0].type == "quad"
    np.testing.assert_allclose(
        a.point_data["DISPLACEMENT"], oracle_mesh().point_data["DISPLACEMENT"]
    )
    np.testing.assert_array_equal(a.cell_data["STRESS:X1"][0], [101.0])
    assert [r.name for r in a.regions] == ["PSHELL_7"]


# -- a mesh with every card -----------------------------------------------------------

TYPES = {
    "vertex": 1,
    "line": 2,
    "triangle": 3,
    "triangle6": 6,
    "quad": 4,
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


def every_card(rng=None):
    rng = rng or np.random.default_rng(3)
    points = rng.random((30, 3))
    cells = [(t, (np.arange(k)[None, :] + 2) % 30) for t, k in TYPES.items()]
    mesh = meshioplusplus.Mesh(points, cells)
    pid = {"vertex": -1, "line": 4}
    mesh.cell_data["nastran:pid"] = [
        np.array([pid.get(t, 2 if t[:3] in ("tet", "pyr", "wed", "hex") else 1)])
        for t in TYPES
    ]
    return mesh


def test_every_card_survives_a_round_trip(engine, tmp_path):
    mesh = every_card()
    path = tmp_path / "cards.op2"
    write_with(engine, path, mesh)
    back = meshioplusplus.nastran_op2.read(path)
    assert {b.type for b in back.cells} == set(TYPES)
    by_type = {b.type: b.data for b in back.cells}
    for block in mesh.cells:
        np.testing.assert_array_equal(by_type[block.type], block.data)
    np.testing.assert_allclose(back.points, mesh.points, rtol=1e-6)
    assert sorted(r.name for r in back.regions) == ["PROD_4", "PSHELL_1", "PSOLID_2"]


def test_the_engines_write_the_same_bytes(tmp_path):
    mesh = every_card()
    mesh.point_data["DISPLACEMENT"] = np.random.default_rng(1).random((30, 3))
    mesh.point_data["VELOCITY"] = np.random.default_rng(2).random((30, 3))
    mesh.point_data["VELOCITY_ROT"] = np.random.default_rng(4).random((30, 3))
    mesh.point_data["TEMPERATURE"] = np.random.default_rng(5).random(30)
    mesh.point_data["SPC_FORCE"] = np.full((30, 3), np.nan)
    mesh.point_data["SPC_FORCE"][3] = [1.0, 2.0, 3.0]
    rng = np.random.default_rng(6)
    for group in ("STRESS", "STRAIN"):
        for m in ("X1", "VON_MISES1", "X2", "X", "TXY", "PRESSURE", "VON_MISES", "A"):
            mesh.cell_data[f"{group}:{m}"] = [rng.random(1) for _ in TYPES]
    mesh.cell_data["STRESS:X@corner"] = [rng.random((1, 8)) for _ in TYPES]
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_the_public_write_is_found_by_name_and_extension(tmp_path):
    mesh = oracle_mesh()
    meshioplusplus.write(tmp_path / "a.op2", mesh)
    meshioplusplus.write(tmp_path / "b.bin", mesh, file_format="nastran_op2")
    assert (tmp_path / "a.op2").read_bytes() == oracle_bytes()
    assert (tmp_path / "b.bin").read_bytes() == oracle_bytes()
    assert meshioplusplus.read(tmp_path / "b.bin").cells[0].type == "quad"
    with open(tmp_path / "c.op2", "wb") as fh:
        meshioplusplus.nastran_op2.write(fh, mesh)
    assert (tmp_path / "c.op2").read_bytes() == oracle_bytes()


# -- results ---------------------------------------------------------------------------


def results_mesh():
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1]], dtype=float
    )
    cells = [
        ("line", np.array([[4, 5]])),
        ("quad", np.array([[0, 1, 2, 3]])),
        ("tetra", np.array([[0, 1, 3, 4]])),
    ]
    mesh = meshioplusplus.Mesh(points, cells)
    mesh.cell_data["nastran:pid"] = [np.array([4]), np.array([1]), np.array([2])]
    mesh.point_data["DISPLACEMENT"] = np.arange(18, dtype=float).reshape(6, 3) / 4
    mesh.point_data["TEMPERATURE"] = np.arange(6, dtype=float) * 10
    mesh.cell_data["STRESS:A"] = [
        np.array([5.5]),
        np.array([np.nan]),
        np.array([np.nan]),
    ]
    mesh.cell_data["STRESS:X1"] = [
        np.array([np.nan]),
        np.array([7.0]),
        np.array([np.nan]),
    ]
    mesh.cell_data["STRESS:X"] = [
        np.array([np.nan]),
        np.array([np.nan]),
        np.array([9.0]),
    ]
    mesh.cell_data["STRAIN:X1"] = [
        np.array([np.nan]),
        np.array([0.5]),
        np.array([np.nan]),
    ]
    return mesh


def test_results_come_back_by_element_and_point(engine, tmp_path):
    path = tmp_path / "results.op2"
    write_with(engine, path, results_mesh())
    back = py_op2.read(path)
    types = [b.type for b in back.cells]
    eid = {t: back.cell_data["nastran:eid"][i][0] for i, t in enumerate(types)}
    assert eid == {"line": 1, "quad": 2, "tetra": 3}
    block = {t: i for i, t in enumerate(types)}
    assert back.cell_data["STRESS:A"][block["line"]][0] == 5.5
    assert back.cell_data["STRESS:X1"][block["quad"]][0] == 7.0
    assert back.cell_data["STRESS:X"][block["tetra"]][0] == 9.0
    assert back.cell_data["STRAIN:X1"][block["quad"]][0] == 0.5
    # a solid's corners repeat its centre value
    np.testing.assert_array_equal(
        back.cell_data["STRESS:X@corner"][block["tetra"]], [[9.0] * 4]
    )
    np.testing.assert_allclose(
        back.point_data["DISPLACEMENT"], results_mesh().point_data["DISPLACEMENT"]
    )
    np.testing.assert_array_equal(back.point_data["TEMPERATURE"], np.arange(6) * 10.0)


def test_nan_points_and_cells_are_left_out(engine, tmp_path):
    mesh = results_mesh()
    mesh.point_data["DISPLACEMENT"][2] = np.nan
    path = tmp_path / "nan.op2"
    write_with(engine, path, mesh)
    back = py_op2.read(path)
    assert np.isnan(back.point_data["DISPLACEMENT"][2]).all()
    assert not np.isnan(back.point_data["DISPLACEMENT"][3]).any()


def test_a_shear_or_octahedral_column_is_named_by_the_mesh(engine, tmp_path):
    mesh = results_mesh()
    del (
        mesh.cell_data["STRESS:X1"],
        mesh.cell_data["STRESS:X"],
        mesh.cell_data["STRAIN:X1"],
    )
    mesh.cell_data["STRESS:MAX_SHEAR1"] = [
        np.array([np.nan]),
        np.array([3.0]),
        np.array([np.nan]),
    ]
    mesh.cell_data["STRESS:OCT_SHEAR"] = [
        np.array([np.nan]),
        np.array([np.nan]),
        np.array([4.0]),
    ]
    path = tmp_path / "shear.op2"
    write_with(engine, path, mesh)
    back = py_op2.read(path)
    types = [b.type for b in back.cells]
    assert back.cell_data["STRESS:MAX_SHEAR1"][types.index("quad")][0] == 3.0
    assert back.cell_data["STRESS:OCT_SHEAR"][types.index("tetra")][0] == 4.0


def test_what_has_no_place_in_the_file_is_dropped_with_a_warning(
    engine, tmp_path, capfd
):
    mesh = results_mesh()
    mesh.point_data["pressure"] = np.zeros(6)
    mesh.field_data["note"] = np.array([1])
    mesh.cells.append(meshioplusplus.CellBlock("polygon", np.array([[0, 1, 2, 3, 4]])))
    path = tmp_path / "dropped.op2"
    write_with(engine, path, mesh)
    err = " ".join(capfd.readouterr().err.split())
    assert "not written" in err
    for name in ("pressure", "note", "polygon"):
        assert name in err
    assert "pressure" not in py_op2.read(path).point_data


# -- steps -----------------------------------------------------------------------------


def step_mesh(analysis, t, subcase=1):
    mesh = results_mesh()
    mesh.point_data["DISPLACEMENT"] = mesh.point_data["DISPLACEMENT"] * t
    mesh.cell_data["STRESS:A"][0] = np.array([5.5 * t])
    mesh.field_data["nastran:subcase"] = np.array([subcase])
    mesh.field_data["nastran:analysis"] = np.array([analysis])
    mesh.field_data["nastran:mode"] = np.array([int(t)])
    return mesh


@pytest.mark.parametrize("analysis", [1, 2, 5, 6, 8])
def test_a_series_reads_back_step_by_step(engine, tmp_path, analysis):
    path = tmp_path / "series.op2"
    times = [0.0] if analysis == 1 else [1.0, 2.0, 4.0]
    write_series_with(engine, path, [(t, step_mesh(analysis, t or 1.0)) for t in times])
    assert py_op2.time_values(path) == times
    for index, t in enumerate(times):
        back = py_op2.read(path, time_step=index)
        assert back.field_data["meshio:time"][0] == t
        assert back.field_data["nastran:analysis"][0] == analysis
        assert back.field_data["nastran:subcase"][0] == 1
        scale = t or 1.0
        np.testing.assert_allclose(
            back.point_data["DISPLACEMENT"][5], np.array([15, 16, 17]) / 4 * scale
        )
        line = [b.type for b in back.cells].index("line")
        assert back.cell_data["STRESS:A"][line][0] == 5.5 * scale


def test_the_engines_write_the_same_series(tmp_path):
    steps = [(t, step_mesh(6, t)) for t in (0.5, 1.5, 2.5)]
    write_series_with("core", tmp_path / "core.op2", steps)
    write_series_with("python", tmp_path / "python.op2", steps)
    assert (tmp_path / "core.op2").read_bytes() == (
        tmp_path / "python.op2"
    ).read_bytes()


def test_subcases_are_steps_of_their_own(engine, tmp_path):
    path = tmp_path / "subcases.op2"
    write_series_with(
        engine,
        path,
        [(0.0, step_mesh(1, 1.0, subcase=1)), (0.0, step_mesh(1, 2.0, subcase=2))],
    )
    assert py_op2.time_values(path) == [0.0, 0.0]
    assert py_op2.read(path, time_step=1).field_data["nastran:subcase"][0] == 2


def test_write_sequence_uses_the_series_writer(tmp_path):
    steps = [(t, step_mesh(6, t)) for t in (1.0, 2.0)]
    path = tmp_path / "seq.op2"
    meshioplusplus.write_sequence(path, steps, file_format="nastran_op2")
    assert py_op2.time_values(path) == [1.0, 2.0]
    assert len(meshioplusplus.read_metadata(path)["time_values"]) == 2
    assert [t for t, _ in meshioplusplus.read_sequence(path)] == [1.0, 2.0]


def test_a_rewrite_is_stable(tmp_path):
    # What a read adds (zero members, the corner arrays) is written back the same.
    first, second, third = (tmp_path / f"{n}.op2" for n in "abc")
    meshioplusplus.write_sequence(
        first, [(t, step_mesh(6, t)) for t in (1.0, 2.0, 3.0)]
    )
    meshioplusplus.write_sequence(second, meshioplusplus.read_sequence(first))
    meshioplusplus.write_sequence(third, meshioplusplus.read_sequence(second))
    assert second.read_bytes() == third.read_bytes()
    assert py_op2.time_values(third) == [1.0, 2.0, 3.0]


def test_a_one_step_series_is_a_plain_write(engine, tmp_path):
    mesh = step_mesh(6, 1.5)
    plain, series = tmp_path / "plain.op2", tmp_path / "series.op2"
    mesh.field_data["meshio:time"] = np.array([1.5])
    write_with(engine, plain, mesh)
    write_series_with(engine, series, [(1.5, mesh)])
    assert plain.read_bytes() == series.read_bytes()


@pytest.mark.parametrize(
    "make, message",
    [
        (lambda m: m.point_data.pop("DISPLACEMENT"), "result arrays differ"),
        (lambda m: m.cell_data.pop("STRESS:A"), "result arrays differ"),
        (
            lambda m: m.points.__setitem__(slice(None), 0.0) or m.cells.pop(),
            "points or cells",
        ),
    ],
)
def test_a_series_refuses_a_step_with_other_arrays_or_cells(
    engine, tmp_path, make, message
):
    second = step_mesh(6, 2.0)
    make(second)
    with pytest.raises(Exception, match=message):
        write_series_with(
            engine, tmp_path / "bad.op2", [(1.0, step_mesh(6, 1.0)), (2.0, second)]
        )


def test_a_series_refuses_two_steps_with_one_key(engine, tmp_path):
    steps = [(1.0, step_mesh(6, 1.0)), (1.0, step_mesh(6, 1.0))]
    with pytest.raises(Exception, match="of an earlier step"):
        write_series_with(engine, tmp_path / "twice.op2", steps)


def test_a_series_with_no_results_is_the_model(engine, tmp_path):
    path = tmp_path / "model.op2"
    write_series_with(engine, path, [(0.0, every_card()), (0.0, every_card())])
    assert py_op2.time_values(path) == []
    assert len(py_op2.read(path).cells) == len(TYPES)


# -- refusals --------------------------------------------------------------------------


def test_repeated_element_ids_are_refused(engine, tmp_path):
    mesh = results_mesh()
    mesh.cell_data["nastran:eid"] = [np.array([1]), np.array([1]), np.array([2])]
    with pytest.raises(Exception, match="not unique"):
        write_with(engine, tmp_path / "x.op2", mesh)


def test_a_property_id_shared_by_two_families_is_refused(engine, tmp_path):
    mesh = results_mesh()
    mesh.cell_data["nastran:pid"] = [np.array([4]), np.array([1]), np.array([1])]
    with pytest.raises(Exception, match="need ids of their own"):
        write_with(engine, tmp_path / "x.op2", mesh)


def test_a_cell_on_a_missing_point_is_refused(engine, tmp_path):
    mesh = results_mesh()
    mesh.cells[1].data[0, 0] = 99
    with pytest.raises(Exception, match="point the mesh does not have"):
        write_with(engine, tmp_path / "x.op2", mesh)


def test_an_array_of_the_wrong_size_is_refused(engine, tmp_path):
    mesh = results_mesh()
    mesh.point_data["DISPLACEMENT"] = np.zeros((6, 2))
    with pytest.raises(Exception, match="3 component"):
        write_with(engine, tmp_path / "x.op2", mesh)


def test_a_step_code_that_is_not_positive_is_refused(engine, tmp_path):
    mesh = results_mesh()
    mesh.field_data["nastran:subcase"] = np.array([0])
    with pytest.raises(Exception, match="positive"):
        write_with(engine, tmp_path / "x.op2", mesh)


def test_a_mesh_with_no_results_is_a_model(engine, tmp_path):
    mesh = every_card()
    path = tmp_path / "model.op2"
    write_with(engine, path, mesh)
    back = py_op2.read(path)
    assert back.field_data.get("meshio:time") is None
    assert py_op2.time_values(path) == []


def test_the_python_twin_raises_a_write_error(tmp_path):
    with pytest.raises(WriteError):
        py_write.write(str(tmp_path / "x.op2"), _bad_mesh())


def _bad_mesh():
    mesh = results_mesh()
    mesh.cell_data["nastran:eid"] = [np.array([0]), np.array([1]), np.array([2])]
    return mesh


# -- an independent reader ---------------------------------------------------------------


def test_pynastran_reads_the_written_file(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2", exc_type=ImportError)
    path = tmp_path / "pn.op2"
    meshioplusplus.write(path, step_mesh(6, 2.0))
    model = op2.OP2(debug=False)
    model.read_op2(str(path), build_dataframe=False)
    disp = model.displacements[1]
    np.testing.assert_allclose(
        disp.data[0, :, :3], step_mesh(6, 2.0).point_data["DISPLACEMENT"]
    )
    assert model.crod_stress[1].data[0, 0, 0] == 11.0
    assert model.cquad4_stress[1].data[0, 0, 1] == 7.0
