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


_GRID_IDS = [40, 7, 23, 5, 19, 11]


def _numbered_mesh():
    mesh = results_mesh()
    mesh.point_data["nastran:nid"] = np.array(_GRID_IDS, dtype=np.int64)
    return mesh


def test_grid_ids_other_than_one_to_n_round_trip(engine, tmp_path):
    """``point_data["nastran:nid"]`` numbers the GRIDs; the elements, the nodal
    tables and the solids' corner GRIDs use those ids. GEOM1 holds the GRIDs by id,
    so the points come back in id order."""
    path = tmp_path / "ids.op2"
    mesh = _numbered_mesh()
    write_with(engine, path, mesh)
    order = np.argsort(_GRID_IDS)
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        assert back.point_data["nastran:nid"].tolist() == sorted(_GRID_IDS)
        np.testing.assert_allclose(back.points, mesh.points[order], rtol=1e-6)
        by_type = {b.type: b for b in back.cells}
        for a in mesh.cells:  # the same GRIDs, by id
            b = by_type[a.type]
            np.testing.assert_array_equal(
                np.asarray(_GRID_IDS)[a.data],
                back.point_data["nastran:nid"][b.data],
            )
        np.testing.assert_allclose(
            back.point_data["DISPLACEMENT"], mesh.point_data["DISPLACEMENT"][order]
        )
        np.testing.assert_array_equal(
            back.point_data["TEMPERATURE"], mesh.point_data["TEMPERATURE"][order]
        )
        types = [b.type for b in back.cells]
        np.testing.assert_array_equal(
            back.cell_data["STRESS:X@corner"][types.index("tetra")], [[9.0] * 4]
        )


def test_the_engines_write_the_same_bytes_with_grid_ids(tmp_path):
    mesh = _numbered_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


@pytest.mark.parametrize("ids", [[1, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5]])
def test_bad_grid_ids_are_refused(engine, tmp_path, ids):
    mesh = results_mesh()
    mesh.point_data["nastran:nid"] = np.array(ids)
    with pytest.raises((WriteError, RuntimeError)):
        write_with(engine, tmp_path / "bad.op2", mesh)


def test_pynastran_reads_the_grid_ids(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2", exc_type=ImportError)
    path = tmp_path / "pn_ids.op2"
    _core.nastran_op2_write(str(path), _numbered_mesh())
    model = op2.OP2(debug=False)
    model.read_op2(str(path), build_dataframe=False)
    np.testing.assert_array_equal(
        model.displacements[1].node_gridtype[:, 0], sorted(_GRID_IDS)
    )


def _material_mesh():
    mesh = results_mesh()  # PROD 4, PSHELL 1, PSOLID 2
    mesh.field_data.update(
        {
            "nastran:PSHELL_1:T": np.array([0.25]),
            "nastran:PSHELL_1:MID1": np.array([7]),
            "nastran:PROD_4:A": np.array([2.5]),
            "nastran:PROD_4:J": np.array([0.5]),
            "nastran:MAT1_7:E": np.array([2.1e5]),
            "nastran:MAT1_7:NU": np.array([0.29]),
            "nastran:MAT1_7:RHO": np.array([7.8e-9]),
            "nastran:MAT1_1:E": np.array([70000.0]),
        }
    )
    return mesh


def test_property_and_material_values_round_trip(engine, tmp_path):
    """``nastran:<CARD>_<id>:<MEMBER>`` sets the property and material values, and
    the readers return them."""
    path = tmp_path / "values.op2"
    write_with(engine, path, _material_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        data = back.field_data
        assert data["nastran:PSHELL_1:T"][0] == pytest.approx(0.25)
        assert data["nastran:PSHELL_1:MID1"].tolist() == [7]
        assert data["nastran:PROD_4:A"][0] == pytest.approx(2.5)
        assert data["nastran:PROD_4:J"][0] == pytest.approx(0.5)
        assert data["nastran:PROD_4:MID"].tolist() == [1]
        assert data["nastran:PSOLID_2:MID"].tolist() == [1]
        assert data["nastran:MAT1_7:E"][0] == pytest.approx(2.1e5)
        assert data["nastran:MAT1_7:NU"][0] == pytest.approx(0.29)
        assert data["nastran:MAT1_7:RHO"][0] == pytest.approx(7.8e-9, rel=1e-6)
        assert data["nastran:MAT1_1:E"][0] == pytest.approx(70000.0)
        assert data["nastran:MAT1_1:NU"][0] == pytest.approx(0.3)


def test_the_engines_write_the_same_bytes_with_values(tmp_path):
    mesh = _material_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_default_values_are_the_old_placeholders(engine, tmp_path):
    path = tmp_path / "plain.op2"
    write_with(engine, path, results_mesh())
    data = meshioplusplus.nastran_op2.read(path).field_data
    assert data["nastran:PSHELL_1:T"][0] == 1.0
    assert data["nastran:MAT1_1:E"][0] == 1.0
    assert data["nastran:MAT1_1:NU"][0] == pytest.approx(0.3)


def _complex_mesh():
    mesh = results_mesh()
    rng = np.random.default_rng(11)
    mesh.point_data["EIGENVECTOR_real"] = rng.random((6, 3))
    mesh.point_data["EIGENVECTOR_imag"] = rng.random((6, 3))
    mesh.point_data["EIGENVECTOR_ROT_real"] = rng.random((6, 3))
    mesh.point_data["EIGENVECTOR_ROT_imag"] = rng.random((6, 3))
    mesh.point_data["DISPLACEMENT_real"] = rng.random((6, 3))
    mesh.point_data["DISPLACEMENT_imag"] = rng.random((6, 3))
    mesh.field_data["nastran:analysis"] = np.array([5])
    mesh.field_data["meshio:time"] = np.array([12.5])
    return mesh


def test_complex_nodal_tables_round_trip(engine, tmp_path):
    """``<name>_real``/``_imag`` (and ``_ROT_``) are a complex table: table code
    1000 + the real one, format code 2, 14 words a node."""
    path = tmp_path / "complex.op2"
    mesh = _complex_mesh()
    write_with(engine, path, mesh)
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        for name in (
            "EIGENVECTOR_real",
            "EIGENVECTOR_imag",
            "EIGENVECTOR_ROT_real",
            "EIGENVECTOR_ROT_imag",
            "DISPLACEMENT_real",
            "DISPLACEMENT_imag",
        ):
            np.testing.assert_allclose(
                back.point_data[name], mesh.point_data[name], rtol=1e-6, err_msg=name
            )
        # the real table of the same name is untouched (a real DISPLACEMENT)
        np.testing.assert_allclose(
            back.point_data["DISPLACEMENT"], mesh.point_data["DISPLACEMENT"], rtol=1e-6
        )


def test_the_engines_write_the_same_bytes_with_complex_tables(tmp_path):
    mesh = _complex_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_pynastran_reads_the_complex_table(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2", exc_type=ImportError)
    path = tmp_path / "pn_complex.op2"
    mesh = _complex_mesh()
    del mesh.point_data["EIGENVECTOR_real"], mesh.point_data["EIGENVECTOR_imag"]
    del mesh.point_data["EIGENVECTOR_ROT_real"], mesh.point_data["EIGENVECTOR_ROT_imag"]
    del mesh.point_data["DISPLACEMENT"]
    mesh.cell_data = {
        k: v for k, v in mesh.cell_data.items() if k.startswith("nastran")
    }
    _core.nastran_op2_write(str(path), mesh)
    model = op2.OP2(debug=False)
    model.read_op2(str(path), build_dataframe=False)
    table = model.displacements[1]
    np.testing.assert_allclose(
        table.data[0, :, :3].real, mesh.point_data["DISPLACEMENT_real"], rtol=1e-6
    )
    np.testing.assert_allclose(
        table.data[0, :, :3].imag, mesh.point_data["DISPLACEMENT_imag"], rtol=1e-6
    )


def _force_mesh():
    mesh = results_mesh()  # line (CROD 1), quad (CQUAD4 2), tetra (CTETRA 3)
    mesh.cell_data["ELEMENT_FORCE:AF"] = [
        np.array([12.5]),
        np.array([np.nan]),
        np.array([np.nan]),
    ]
    mesh.cell_data["ELEMENT_FORCE:TRQ"] = [
        np.array([-1.5]),
        np.array([np.nan]),
        np.array([np.nan]),
    ]
    for k, name in enumerate(["MX", "MY", "MXY", "BMX", "BMY", "BMXY", "TX", "TY"]):
        mesh.cell_data["ELEMENT_FORCE:" + name] = [
            np.array([np.nan]),
            np.array([10.0 + k]),
            np.array([np.nan]),
        ]
    return mesh


def test_element_forces_of_rods_and_plates_round_trip(engine, tmp_path):
    path = tmp_path / "forces.op2"
    write_with(engine, path, _force_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        types = [b.type for b in back.cells]
        rod, plate = types.index("line"), types.index("quad")
        assert back.cell_data["ELEMENT_FORCE:AF"][rod][0] == 12.5
        assert back.cell_data["ELEMENT_FORCE:TRQ"][rod][0] == -1.5
        for k, name in enumerate(["MX", "MY", "MXY", "BMX", "BMY", "BMXY", "TX", "TY"]):
            assert back.cell_data["ELEMENT_FORCE:" + name][plate][0] == 10.0 + k


def test_the_engines_write_the_same_bytes_with_element_forces(tmp_path):
    mesh = _force_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_pynastran_reads_the_element_forces(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2", exc_type=ImportError)
    path = tmp_path / "pn_forces.op2"
    _core.nastran_op2_write(str(path), _force_mesh())
    model = op2.OP2(debug=False)
    model.read_op2(str(path), build_dataframe=False)
    rod = model.crod_force[1]
    assert rod.data[0, 0, 0] == 12.5 and rod.data[0, 0, 1] == -1.5
    plate = model.cquad4_force[1]
    np.testing.assert_allclose(plate.data[0, 0, :8], 10.0 + np.arange(8))


@pytest.mark.parametrize(
    "points_only,arrays,expect_nid",
    [
        (True, None, False),
        (False, ["DISPLACEMENT"], False),
        (False, ["DISPLACEMENT", "nastran:nid"], True),
        (False, None, True),
    ],
)
def test_grid_id_array_follows_the_read_options(
    engine, tmp_path, points_only, arrays, expect_nid
):
    """``nastran:nid`` is point data: ``points_only`` drops it, and ``arrays`` keeps it
    only when it is named (both readers)."""
    path = tmp_path / "read_options.op2"
    write_with(engine, path, _numbered_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path, points_only=points_only, arrays=arrays)
        assert ("nastran:nid" in back.point_data) == expect_nid
        if expect_nid:
            np.testing.assert_array_equal(
                back.point_data["nastran:nid"], sorted(_GRID_IDS)
            )


# -- springs and dampers -----------------------------------------------------------------


def _spring_mesh():
    """Springs CELAS1 30 and 31 between GRIDs 1-2 and 2-3 (PELAS 5), one grounded at
    GRID 4 (its second component is 0), a CDAMP1 50 between 2 and 4 (PDAMP 6), and a
    card 99 (40) the writer does not know."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
    )
    mesh = meshioplusplus.Mesh(
        points,
        [
            ("line", np.array([[0, 1], [1, 2]])),
            ("vertex", np.array([[3]])),
            ("line", np.array([[1, 3]])),
            ("line", np.array([[0, 2]])),
        ],
    )
    mesh.cell_data["nastran:eid"] = [
        np.array([30, 31]),
        np.array([32]),
        np.array([50]),
        np.array([40]),
    ]
    mesh.cell_data["nastran:pid"] = [
        np.array([5, 5]),
        np.array([5]),
        np.array([6]),
        np.array([7]),
    ]
    mesh.cell_data["nastran:card"] = [
        np.array([6, 6]),
        np.array([6]),
        np.array([2]),
        np.array([99]),
    ]
    mesh.cell_data["nastran:c1"] = [
        np.array([1, 1]),
        np.array([2]),
        np.array([1]),
        np.array([0]),
    ]
    mesh.cell_data["nastran:c2"] = [
        np.array([1, 2]),
        np.array([0]),
        np.array([1]),
        np.array([0]),
    ]
    mesh.field_data["nastran:PELAS_5:K"] = np.array([1000.0])
    mesh.field_data["nastran:PELAS_5:GE"] = np.array([0.02])
    mesh.field_data["nastran:PDAMP_6:B"] = np.array([3.5])
    return mesh


def _elements(mesh):
    """``{eid: (nodes, pid, card, c1, c2)}`` of a mesh's element blocks."""
    out = {}
    columns = zip(
        mesh.cells,
        mesh.cell_data["nastran:eid"],
        mesh.cell_data["nastran:pid"],
        mesh.cell_data["nastran:card"],
        mesh.cell_data["nastran:c1"],
        mesh.cell_data["nastran:c2"],
    )
    for block, eids, pids, cards, c1, c2 in columns:
        for row in range(len(eids)):
            out[int(eids[row])] = (
                block.data[row].tolist(),
                int(pids[row]),
                int(cards[row]),
                int(c1[row]),
                int(c2[row]),
            )
    return out


def test_springs_and_dampers_round_trip(engine, tmp_path, capfd):
    path = tmp_path / "springs.op2"
    write_with(engine, path, _spring_mesh())
    capfd.readouterr()  # the CBAR's warning
    expect = {
        30: ([0, 1], 5, 6, 1, 1),
        31: ([1, 2], 5, 6, 1, 2),
        32: ([3], 5, 6, 2, 0),
        50: ([1, 3], 6, 2, 1, 1),
    }
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        assert _elements(back) == expect
        assert sorted(r.name for r in back.regions) == ["PDAMP_6", "PELAS_5"]
        assert back.field_data["nastran:PELAS_5:K"][0] == 1000.0
        assert back.field_data["nastran:PELAS_5:GE"][0] == pytest.approx(0.02)
        assert back.field_data["nastran:PELAS_5:S"][0] == 0.0
        assert back.field_data["nastran:PDAMP_6:B"][0] == 3.5


def test_the_engines_write_the_same_bytes_with_springs(tmp_path):
    mesh = _spring_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_cards_the_writer_does_not_hold_are_left_out_with_a_warning(
    engine, tmp_path, capfd
):
    mesh = _spring_mesh()
    # A CSHEAR is a quad: on a triangle its number names a card that is left out.
    mesh.cells.append(meshioplusplus.CellBlock("triangle", np.array([[0, 1, 2]])))
    mesh.cell_data["nastran:eid"].append(np.array([60]))
    mesh.cell_data["nastran:pid"].append(np.array([8]))
    mesh.cell_data["nastran:card"].append(np.array([31]))  # CSHEAR
    for name in ("nastran:c1", "nastran:c2"):
        mesh.cell_data[name].append(np.array([0]))
    path = tmp_path / "left_out.op2"
    write_with(engine, path, mesh)
    err = " ".join(capfd.readouterr().err.split())
    assert "element cards not written yet" in err
    assert "card 99" in err and "CSHEAR" in err
    back = py_op2.read(path)
    assert sorted(_elements(back)) == [30, 31, 32, 50]


def test_pynastran_reads_the_springs_and_dampers(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_springs.op2"
    _core.nastran_op2_write(str(path), _spring_mesh())
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    assert sorted(model.elements) == [30, 31, 32, 50]
    assert model.elements[30].type == "CELAS1" and model.elements[30].pid == 5
    assert model.elements[32].node_ids == [4, None]
    assert (model.elements[32].c1, model.elements[32].c2) == (2, 0)
    assert model.elements[50].type == "CDAMP1" and model.elements[50].node_ids == [2, 4]
    assert model.properties[5].type == "PELAS"
    assert model.properties[5].k == 1000.0
    assert model.properties[5].ge == pytest.approx(0.02)
    assert model.properties[6].type == "PDAMP" and model.properties[6].b == 3.5


# -- CBAR and PBAR -----------------------------------------------------------------------


def _bar_mesh():
    """Two CBARs of PBAR 3: the first oriented by X = (0, 0, 1), the second by GRID 3 as
    its G0, with pins (1, 2) and offsets at both ends."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    )
    mesh = meshioplusplus.Mesh(points, [("line", np.array([[0, 1], [1, 2]]))])
    mesh.cell_data["nastran:eid"] = [np.array([24, 25])]
    mesh.cell_data["nastran:pid"] = [np.array([3, 3])]
    mesh.cell_data["nastran:card"] = [np.array([24, 24])]
    mesh.cell_data["nastran:x"] = [np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0]])]
    mesh.cell_data["nastran:g0"] = [np.array([0, 3])]
    mesh.cell_data["nastran:pin"] = [np.array([[0, 0], [1, 2]])]
    mesh.cell_data["nastran:w"] = [
        np.array([[0.1, 0.2, 0.3, 0.4, 0.5, 0.6], [0.0] * 6])
    ]
    mesh.field_data["nastran:PBAR_3:A"] = np.array([2.5])
    mesh.field_data["nastran:PBAR_3:I1"] = np.array([0.25])
    mesh.field_data["nastran:PBAR_3:J"] = np.array([0.5])
    mesh.field_data["nastran:PBAR_3:MID"] = np.array([1], dtype=np.int64)
    return mesh


def test_cbars_and_their_pbars_round_trip(engine, tmp_path):
    path = tmp_path / "cbar.op2"
    write_with(engine, path, _bar_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        assert np.concatenate(back.cell_data["nastran:eid"]).tolist() == [24, 25]
        assert np.concatenate(back.cell_data["nastran:card"]).tolist() == [24, 24]
        np.testing.assert_array_equal(
            np.concatenate(back.cell_data["nastran:g0"]), [0, 3]
        )
        np.testing.assert_array_equal(
            np.concatenate(back.cell_data["nastran:pin"]), [[0, 0], [1, 2]]
        )
        np.testing.assert_allclose(
            np.concatenate(back.cell_data["nastran:x"]), [[0, 0, 1], [0, 0, 0]]
        )
        w = np.concatenate(back.cell_data["nastran:w"])
        np.testing.assert_allclose(w[0], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6], rtol=1e-6)
        assert not w[1].any()
        assert back.field_data["nastran:PBAR_3:A"][0] == 2.5
        assert back.field_data["nastran:PBAR_3:I1"][0] == 0.25
        assert back.field_data["nastran:PBAR_3:J"][0] == 0.5
        assert back.field_data["nastran:PBAR_3:MID"][0] == 1
        assert [r.name for r in back.regions] == ["PBAR_3"]


def test_the_engines_write_the_same_bytes_with_cbars(tmp_path):
    mesh = _bar_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_cbar_needs_its_orientation(engine, tmp_path):
    mesh = _bar_mesh()
    for name in ("nastran:x", "nastran:g0"):
        del mesh.cell_data[name]
    with pytest.raises(Exception, match="orientation"):
        write_with(engine, tmp_path / "x.op2", mesh)


def test_pynastran_reads_the_cbars_and_pbars(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_cbar.op2"
    _core.nastran_op2_write(str(path), _bar_mesh())
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    first, second = model.elements[24], model.elements[25]
    assert first.type == second.type == "CBAR" and first.pid == 3
    assert first.node_ids == [1, 2] and second.node_ids == [2, 3]
    np.testing.assert_allclose(first.x, [0, 0, 1])
    assert first.g0 is None
    assert second.g0 == 3 and second.x is None
    assert (second.pa, second.pb) == (1, 2)
    np.testing.assert_allclose(first.wa, [0.1, 0.2, 0.3], rtol=1e-6)
    np.testing.assert_allclose(first.wb, [0.4, 0.5, 0.6], rtol=1e-6)
    prop = model.properties[3]
    assert prop.type == "PBAR" and prop.mid == 1
    assert prop.A == 2.5 and prop.i1 == 0.25 and prop.j == 0.5


# -- PCOMP composites --------------------------------------------------------------------


def _composite_mesh(sym=0):
    """A quad (CQUAD4 7) whose property 9 is a three-ply PCOMP, its plies' materials 1
    and 2: Z0 0.5, NSM 0.1, SB 2, failure theory 3 (TSAI), TREF 20, GE 0.02."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    mesh = meshioplusplus.Mesh(points, [("quad", np.array([[0, 1, 2, 3]]))])
    mesh.cell_data["nastran:eid"] = [np.array([7])]
    mesh.cell_data["nastran:pid"] = [np.array([9])]
    mesh.field_data["nastran:PCOMP_9:Z0"] = np.array([0.5])
    mesh.field_data["nastran:PCOMP_9:NSM"] = np.array([0.1])
    mesh.field_data["nastran:PCOMP_9:SB"] = np.array([2.0])
    mesh.field_data["nastran:PCOMP_9:FT"] = np.array([3], dtype=np.int64)
    mesh.field_data["nastran:PCOMP_9:TREF"] = np.array([20.0])
    mesh.field_data["nastran:PCOMP_9:GE"] = np.array([0.02])
    mesh.field_data["nastran:PCOMP_9:SYM"] = np.array([sym], dtype=np.int64)
    mesh.field_data["nastran:PCOMP_9:PLY"] = np.array(
        [[1, 0.1, 0.0, 1], [2, 0.2, 45.0, 0], [1, 0.1, -45.0, 1]], dtype=float
    )
    return mesh


def test_a_pcomp_round_trips_with_its_plies(engine, tmp_path):
    path = tmp_path / "pcomp.op2"
    write_with(engine, path, _composite_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        assert [r.name for r in back.regions] == ["PCOMP_9"]
        np.testing.assert_allclose(
            back.field_data["nastran:PCOMP_9:PLY"],
            [[1, 0.1, 0.0, 1], [2, 0.2, 45.0, 0], [1, 0.1, -45.0, 1]],
            rtol=1e-6,
        )
        assert back.field_data["nastran:PCOMP_9:Z0"][0] == 0.5
        assert back.field_data["nastran:PCOMP_9:SB"][0] == 2.0
        assert back.field_data["nastran:PCOMP_9:FT"][0] == 3
        assert back.field_data["nastran:PCOMP_9:TREF"][0] == 20.0
        assert back.field_data["nastran:PCOMP_9:SYM"][0] == 0
        assert back.field_data["nastran:PCOMP_9:GE"][0] == pytest.approx(0.02)
        assert "nastran:PSHELL_9:T" not in back.field_data


def test_a_symmetric_pcomp_keeps_its_layup_flag(engine, tmp_path):
    path = tmp_path / "pcomp_sym.op2"
    write_with(engine, path, _composite_mesh(sym=1))
    back = py_op2.read(path)
    assert back.field_data["nastran:PCOMP_9:SYM"][0] == 1
    assert back.field_data["nastran:PCOMP_9:PLY"].shape == (3, 4)


def test_the_engines_write_the_same_bytes_with_a_pcomp(tmp_path):
    mesh = _composite_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_pynastran_reads_the_pcomp(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_pcomp.op2"
    _core.nastran_op2_write(str(path), _composite_mesh())
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    prop = model.properties[9]
    assert prop.type == "PCOMP" and prop.nplies == 3
    assert prop.material_ids == [1, 2, 1]
    np.testing.assert_allclose(prop.thicknesses, [0.1, 0.2, 0.1], rtol=1e-6)
    np.testing.assert_allclose(prop.thetas, [0.0, 45.0, -45.0])
    assert prop.souts == ["YES", "NO", "YES"]
    assert (prop.z0, prop.sb, prop.ft) == (0.5, 2.0, "TSAI")
    assert prop.tref == 20.0 and prop.lam is None
    assert sorted(model.materials) == [1, 2]


# -- CELAS2 and CDAMP2 -------------------------------------------------------------------


def _inline_spring_mesh():
    """A CELAS2 30 between GRIDs 1 and 2 (K 1000, GE 0.02), a CDAMP2 31 between 2 and 3
    (B 2.5) and a CELAS2 32 grounded at GRID 4 (K 55). None has a property."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    )
    mesh = meshioplusplus.Mesh(
        points, [("line", np.array([[0, 1], [1, 2]])), ("vertex", np.array([[3]]))]
    )
    mesh.cell_data["nastran:eid"] = [np.array([30, 31]), np.array([32])]
    mesh.cell_data["nastran:pid"] = [np.array([-1, -1]), np.array([-1])]
    mesh.cell_data["nastran:card"] = [np.array([7, 3]), np.array([7])]
    mesh.cell_data["nastran:c1"] = [np.array([1, 2]), np.array([3])]
    mesh.cell_data["nastran:c2"] = [np.array([1, 2]), np.array([0])]
    mesh.cell_data["nastran:k"] = [np.array([1000.0, 0.0]), np.array([55.0])]
    mesh.cell_data["nastran:ge"] = [np.array([0.02, 0.0]), np.array([0.0])]
    mesh.cell_data["nastran:s"] = [np.array([0.0, 0.0]), np.array([0.0])]
    mesh.cell_data["nastran:b"] = [np.array([0.0, 2.5]), np.array([0.0])]
    return mesh


def _rows_by_eid(mesh):
    rows = {}
    for block, eids in enumerate(mesh.cell_data["nastran:eid"]):
        for row, eid in enumerate(eids.tolist()):
            rows[eid] = {
                name: mesh.cell_data[name][block][row]
                for name in (
                    "nastran:card",
                    "nastran:k",
                    "nastran:ge",
                    "nastran:b",
                    "nastran:c1",
                    "nastran:c2",
                )
            }
    return rows


def test_inline_springs_and_dampers_round_trip(engine, tmp_path):
    path = tmp_path / "inline.op2"
    write_with(engine, path, _inline_spring_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        rows = _rows_by_eid(read(path))
        assert sorted(rows) == [30, 31, 32]
        assert rows[30]["nastran:card"] == 7 and rows[30]["nastran:k"] == 1000.0
        assert rows[30]["nastran:ge"] == pytest.approx(0.02)
        assert (rows[30]["nastran:c1"], rows[30]["nastran:c2"]) == (1, 1)
        assert rows[31]["nastran:card"] == 3 and rows[31]["nastran:b"] == 2.5
        assert (rows[31]["nastran:c1"], rows[31]["nastran:c2"]) == (2, 2)
        assert rows[32]["nastran:card"] == 7 and rows[32]["nastran:k"] == 55.0
        assert (rows[32]["nastran:c1"], rows[32]["nastran:c2"]) == (3, 0)


def test_the_engines_write_the_same_bytes_with_inline_springs(tmp_path):
    mesh = _inline_spring_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_celas2_needs_its_stiffness(engine, tmp_path):
    mesh = _inline_spring_mesh()
    del mesh.cell_data["nastran:k"]
    with pytest.raises(Exception, match="nastran:k"):
        write_with(engine, tmp_path / "x.op2", mesh)


def test_pynastran_reads_the_inline_springs_and_dampers(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_inline.op2"
    _core.nastran_op2_write(str(path), _inline_spring_mesh())
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    first, damper, grounded = model.elements[30], model.elements[31], model.elements[32]
    assert first.type == "CELAS2" and first.k == 1000.0
    assert first.ge == pytest.approx(0.02) and first.s == 0.0
    assert (first.c1, first.c2) == (1, 1) and first.node_ids == [1, 2]
    assert damper.type == "CDAMP2" and damper.b == 2.5
    assert damper.node_ids == [2, 3] and (damper.c1, damper.c2) == (2, 2)
    assert grounded.type == "CELAS2" and grounded.k == 55.0
    assert grounded.node_ids == [4, None] and (grounded.c1, grounded.c2) == (3, 0)


# -- coordinate systems -------------------------------------------------------------------


def _cord_mesh():
    """Four GRIDs: one in rectangular system 3 (origin (1, 2, 3), turned a quarter turn
    about z), one in cylindrical system 4 and one in spherical system 5, both at the
    origin, and one in basic coordinates. A CELAS1 spring (PELAS 5) joins GRIDs 1 and 4.
    """
    c30, s30 = np.cos(np.radians(30.0)), np.sin(np.radians(30.0))
    t60, f45 = np.radians(60.0), np.radians(45.0)
    points = np.array(
        [
            [1.0, 3.0, 3.0],  # system 3, local (1, 0, 0)
            [2.0 * c30, 2.0 * s30, 1.0],  # system 4, local (r 2, theta 30, z 1)
            [  # system 5, local (r 2, theta 60, phi 45)
                2.0 * np.sin(t60) * np.cos(f45),
                2.0 * np.sin(t60) * np.sin(f45),
                2.0 * np.cos(t60),
            ],
            [0.5, 0.5, 0.5],
        ]
    )
    mesh = meshioplusplus.Mesh(points, [("line", np.array([[0, 3]]))])
    mesh.cell_data["nastran:eid"] = [np.array([30])]
    mesh.cell_data["nastran:pid"] = [np.array([5])]
    mesh.cell_data["nastran:card"] = [np.array([6])]  # CELAS1
    mesh.cell_data["nastran:c1"] = [np.array([1])]
    mesh.cell_data["nastran:c2"] = [np.array([1])]
    mesh.field_data["nastran:PELAS_5:K"] = np.array([100.0])
    mesh.point_data["nastran:cp"] = np.array([3, 4, 5, 0])
    systems = {
        3: (1, [1.0, 2.0, 3.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        4: (2, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
        5: (3, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
    }
    for cid, (ctype, origin, x, y, z) in systems.items():
        prefix = f"nastran:CORD2_{cid}:"
        mesh.field_data[prefix + "TYPE"] = np.array([ctype], dtype=np.int64)
        mesh.field_data[prefix + "A"] = np.array(origin)
        mesh.field_data[prefix + "X"] = np.array(x)
        mesh.field_data[prefix + "Y"] = np.array(y)
        mesh.field_data[prefix + "Z"] = np.array(z)
    return mesh


def test_grids_in_coordinate_systems_round_trip(engine, tmp_path):
    path = tmp_path / "cord.op2"
    mesh = _cord_mesh()
    write_with(engine, path, mesh)
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        np.testing.assert_allclose(back.points, mesh.points, atol=1e-5)
        np.testing.assert_array_equal(back.point_data["nastran:cp"], [3, 4, 5, 0])
        for cid in (3, 4, 5):
            prefix = f"nastran:CORD2_{cid}:"
            np.testing.assert_array_equal(
                back.field_data[prefix + "TYPE"], mesh.field_data[prefix + "TYPE"]
            )
            for member in "AXYZ":
                np.testing.assert_allclose(
                    back.field_data[prefix + member],
                    mesh.field_data[prefix + member],
                    atol=1e-5,
                )


def test_the_engines_write_the_same_bytes_with_coordinate_systems(tmp_path):
    mesh = _cord_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_grid_in_an_undefined_system_is_refused(engine, tmp_path):
    mesh = _cord_mesh()
    del mesh.field_data["nastran:CORD2_4:Y"]
    with pytest.raises(Exception, match="coordinate system 4"):
        write_with(engine, tmp_path / "x.op2", mesh)


def test_an_output_system_is_written_as_the_grids_cd(engine, tmp_path):
    mesh = _cord_mesh()
    mesh.point_data["nastran:cd"] = np.array([0, 3, 0, 0])
    path = tmp_path / "cd.op2"
    write_with(engine, path, mesh)
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        np.testing.assert_array_equal(read(path).point_data["nastran:cd"], [0, 3, 0, 0])


def test_pynastran_reads_the_coordinate_systems(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_cord.op2"
    mesh = _cord_mesh()
    _core.nastran_op2_write(str(path), mesh)
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    model.cross_reference()
    assert [model.coords[c].type for c in (3, 4, 5)] == ["CORD2R", "CORD2C", "CORD2S"]
    for gid, expected in zip((1, 2, 3, 4), mesh.points):
        np.testing.assert_allclose(model.nodes[gid].get_position(), expected, atol=1e-5)
    assert model.elements[30].node_ids == [1, 4]


# -- CBEAM and PBEAM beams ----------------------------------------------------------------


def _beam_mesh():
    """Two CBEAMs of PBEAM 9 along the x axis: the first oriented by X = (0, 0, 1), the
    second by X = (0, 1, 0) with pins (1, 2) and offsets at both ends. PBEAM 9 has two
    segments, so three stations: the middle one lies between the ends (in slot 1)."""
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    mesh = meshioplusplus.Mesh(points, [("line", np.array([[0, 1], [1, 2]]))])
    mesh.cell_data["nastran:eid"] = [np.array([40, 41])]
    mesh.cell_data["nastran:pid"] = [np.array([9, 9])]
    mesh.cell_data["nastran:card"] = [np.array([54, 54])]
    mesh.cell_data["nastran:x"] = [np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])]
    mesh.cell_data["nastran:g0"] = [np.array([0, 0])]
    mesh.cell_data["nastran:pin"] = [np.array([[0, 0], [1, 2]])]
    mesh.cell_data["nastran:w"] = [
        np.array([[0.0] * 6, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]])
    ]
    mesh.field_data["nastran:PBEAM_9:MID"] = np.array([1], dtype=np.int64)
    mesh.field_data["nastran:PBEAM_9:NSEG"] = np.array([2], dtype=np.int64)
    mesh.field_data["nastran:PBEAM_9:CCF"] = np.array([1], dtype=np.int64)
    mesh.field_data["nastran:PBEAM_9:SO"] = np.array([1.0, 1.0, 1.0])
    mesh.field_data["nastran:PBEAM_9:XXB"] = np.array([0.0, 0.5, 1.0])
    mesh.field_data["nastran:PBEAM_9:A"] = np.array([0.25, 0.3, 0.25])
    mesh.field_data["nastran:PBEAM_9:I1"] = np.array([0.005, 0.006, 0.005])
    mesh.field_data["nastran:PBEAM_9:I2"] = np.array([0.005, 0.006, 0.005])
    mesh.field_data["nastran:PBEAM_9:J"] = np.array([0.01, 0.012, 0.01])
    mesh.field_data["nastran:PBEAM_9:C1"] = np.array([[0.5, 0.5]] * 3)
    mesh.field_data["nastran:PBEAM_9:D2"] = np.array([[-0.5, 0.5]] * 3)
    return mesh


def test_cbeams_and_their_pbeams_round_trip(engine, tmp_path):
    path = tmp_path / "cbeam.op2"
    write_with(engine, path, _beam_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        assert _elements(back) == {
            40: ([0, 1], 9, 54, 0, 0),
            41: ([1, 2], 9, 54, 0, 0),
        }
        np.testing.assert_allclose(
            np.concatenate(back.cell_data["nastran:x"]), [[0, 0, 1], [0, 1, 0]]
        )
        np.testing.assert_array_equal(
            np.concatenate(back.cell_data["nastran:pin"]), [[0, 0], [1, 2]]
        )
        w = np.concatenate(back.cell_data["nastran:w"])
        np.testing.assert_allclose(w[1], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6], rtol=1e-6)
        assert not w[0].any()
        assert back.field_data["nastran:PBEAM_9:MID"][0] == 1
        assert back.field_data["nastran:PBEAM_9:NSEG"][0] == 2
        assert back.field_data["nastran:PBEAM_9:CCF"][0] == 1
        np.testing.assert_allclose(back.field_data["nastran:PBEAM_9:XXB"], [0, 0.5, 1])
        np.testing.assert_allclose(
            back.field_data["nastran:PBEAM_9:A"], [0.25, 0.3, 0.25], rtol=1e-6
        )
        np.testing.assert_allclose(
            back.field_data["nastran:PBEAM_9:J"], [0.01, 0.012, 0.01], rtol=1e-6
        )
        assert back.field_data["nastran:PBEAM_9:C1"].shape == (3, 2)
        np.testing.assert_allclose(
            back.field_data["nastran:PBEAM_9:D2"], [[-0.5, 0.5]] * 3, rtol=1e-6
        )
        assert [r.name for r in back.regions] == ["PBEAM_9"]


def test_the_engines_write_the_same_bytes_with_beams(tmp_path):
    mesh = _beam_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_cbeam_needs_its_orientation(engine, tmp_path):
    mesh = _beam_mesh()
    for name in ("nastran:x", "nastran:g0"):
        del mesh.cell_data[name]
    with pytest.raises(Exception, match="orientation"):
        write_with(engine, tmp_path / "x.op2", mesh)


def _g0_beam_mesh():
    """The beam mesh with a fourth point, GRID 4, as the G0 of its second CBEAM."""
    mesh = _beam_mesh()
    mesh.points = np.vstack([mesh.points, [0.0, 1.0, 0.0]])
    mesh.cell_data["nastran:x"] = [np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0]])]
    mesh.cell_data["nastran:g0"] = [np.array([0, 4])]
    return mesh


def test_a_cbeam_oriented_by_a_grid_round_trips(engine, tmp_path):
    path = tmp_path / "cbeam_g0.op2"
    write_with(engine, path, _g0_beam_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        np.testing.assert_array_equal(
            np.concatenate(back.cell_data["nastran:g0"]), [0, 4]
        )
        np.testing.assert_allclose(
            np.concatenate(back.cell_data["nastran:x"]), [[0, 0, 1], [0, 0, 0]]
        )


def test_the_engines_write_the_same_bytes_with_g0_beams(tmp_path):
    mesh = _g0_beam_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_cbeam_with_only_a_g0_is_written(engine, tmp_path):
    mesh = _g0_beam_mesh()
    del mesh.cell_data["nastran:x"]
    mesh.cell_data["nastran:g0"] = [np.array([4, 4])]
    write_with(engine, tmp_path / "only.op2", mesh)


def test_pynastran_reads_the_g0_cbeams(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_cbeam_g0.op2"
    mesh = _g0_beam_mesh()
    # pyNastran takes a CBEAM's orientation flag from the low bits of word 10 (PA), not
    # word 9 where the record holds it, so it reads a G0 only when PA is 2.
    mesh.cell_data["nastran:pin"] = [np.array([[0, 0], [2, 0]])]
    _core.nastran_op2_write(str(path), mesh)
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    first, second = model.elements[40], model.elements[41]
    np.testing.assert_allclose(first.x, [0, 0, 1])
    assert first.g0 is None
    assert second.g0 == 4 and second.x is None
    assert second.node_ids == [2, 3] and second.pa == 2


def test_a_pbeam_needs_a_segment(engine, tmp_path):
    mesh = _beam_mesh()
    mesh.field_data["nastran:PBEAM_9:NSEG"] = np.array([0], dtype=np.int64)
    with pytest.raises(Exception, match="segments"):
        write_with(engine, tmp_path / "seg.op2", mesh)


def test_pynastran_reads_the_cbeams_and_pbeams(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_cbeam.op2"
    _core.nastran_op2_write(str(path), _beam_mesh())
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    first, second = model.elements[40], model.elements[41]
    assert first.type == second.type == "CBEAM" and first.pid == 9
    assert first.node_ids == [1, 2] and second.node_ids == [2, 3]
    np.testing.assert_allclose(first.x, [0, 0, 1])
    np.testing.assert_allclose(second.x, [0, 1, 0])
    assert (second.pa, second.pb) == (1, 2)
    np.testing.assert_allclose(second.wa, [0.1, 0.2, 0.3], rtol=1e-6)
    np.testing.assert_allclose(second.wb, [0.4, 0.5, 0.6], rtol=1e-6)
    prop = model.properties[9]
    assert prop.type == "PBEAM" and prop.mid == 1
    np.testing.assert_allclose(prop.xxb, [0, 0.5, 1])
    np.testing.assert_allclose(prop.A, [0.25, 0.3, 0.25], rtol=1e-6)
    np.testing.assert_allclose(prop.i1, [0.005, 0.006, 0.005], rtol=1e-5)


# -- CBUSH bushings and PBUSH properties --------------------------------------------------


def _bush_mesh():
    """Two CBUSHes of PBUSH 7 between GRIDs 1-2 and 2-3: the first oriented by X = (1, 0, 0),
    the second by GRID 3 as its G0. The stiffnesses and damping are given per direction, and
    GE by one value, as the 18-word NX layout holds it."""
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    mesh = meshioplusplus.Mesh(points, [("line", np.array([[0, 1], [1, 2]]))])
    mesh.cell_data["nastran:eid"] = [np.array([50, 51])]
    mesh.cell_data["nastran:pid"] = [np.array([7, 7])]
    mesh.cell_data["nastran:card"] = [np.array([26, 26])]
    mesh.cell_data["nastran:x"] = [np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])]
    mesh.cell_data["nastran:g0"] = [np.array([0, 3])]
    mesh.field_data["nastran:PBUSH_7:K"] = np.array([1e6, 2e6, 3e6, 4e5, 5e5, 6e5])
    mesh.field_data["nastran:PBUSH_7:B"] = np.array([0.1, 0.2, 0.3, 0.0, 0.0, 0.0])
    mesh.field_data["nastran:PBUSH_7:GE"] = np.array([0.05])
    mesh.field_data["nastran:PBUSH_7:SA"] = np.array([0.5])
    mesh.field_data["nastran:PBUSH_7:ST"] = np.array([0.6])
    mesh.field_data["nastran:PBUSH_7:EA"] = np.array([0.7])
    mesh.field_data["nastran:PBUSH_7:ET"] = np.array([0.8])
    return mesh


def test_cbushes_and_their_pbushes_round_trip(engine, tmp_path):
    path = tmp_path / "cbush.op2"
    write_with(engine, path, _bush_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        assert _elements(back) == {
            50: ([0, 1], 7, 26, 0, 0),
            51: ([1, 2], 7, 26, 0, 0),
        }
        np.testing.assert_allclose(
            np.concatenate(back.cell_data["nastran:x"]), [[1, 0, 0], [0, 0, 0]]
        )
        np.testing.assert_array_equal(
            np.concatenate(back.cell_data["nastran:g0"]), [0, 3]
        )
        np.testing.assert_allclose(
            back.field_data["nastran:PBUSH_7:K"],
            [1e6, 2e6, 3e6, 4e5, 5e5, 6e5],
            rtol=1e-6,
        )
        np.testing.assert_allclose(
            back.field_data["nastran:PBUSH_7:B"], [0.1, 0.2, 0.3, 0, 0, 0], rtol=1e-6
        )
        # The written record holds the 23-word layout, so GE comes back for each direction.
        np.testing.assert_allclose(
            back.field_data["nastran:PBUSH_7:GE"], [0.05] * 6, rtol=1e-6
        )
        for member, value in (("SA", 0.5), ("ST", 0.6), ("EA", 0.7), ("ET", 0.8)):
            assert back.field_data[f"nastran:PBUSH_7:{member}"][0] == pytest.approx(
                value, rel=1e-6
            )
        assert [r.name for r in back.regions] == ["PBUSH_7"]


def test_the_engines_write_the_same_bytes_with_bushes(tmp_path):
    mesh = _bush_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_cbush_needs_its_orientation(engine, tmp_path):
    mesh = _bush_mesh()
    for name in ("nastran:x", "nastran:g0"):
        del mesh.cell_data[name]
    with pytest.raises(Exception, match="orientation"):
        write_with(engine, tmp_path / "x.op2", mesh)


def test_pynastran_reads_the_cbushes_and_pbushes(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_cbush.op2"
    _core.nastran_op2_write(str(path), _bush_mesh())
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    first, second = model.elements[50], model.elements[51]
    assert first.type == second.type == "CBUSH" and first.pid == 7
    assert first.node_ids == [1, 2] and second.node_ids == [2, 3]
    np.testing.assert_allclose(first.x, [1, 0, 0])
    assert second.g0 == 3 and second.x == [None, None, None]
    prop = model.properties[7]
    assert prop.type == "PBUSH"
    np.testing.assert_allclose(prop.Ki, [1e6, 2e6, 3e6, 4e5, 5e5, 6e5], rtol=1e-6)
    np.testing.assert_allclose(prop.Bi, [0.1, 0.2, 0.3, 0, 0, 0], rtol=1e-6)
    np.testing.assert_allclose(prop.GEi, [0.05] * 6, rtol=1e-6)
    assert (prop.sa, prop.st, prop.ea, prop.et) == pytest.approx((0.5, 0.6, 0.7, 0.8))


def _default_bush_mesh():
    """The bush mesh with no PBUSH values, so the writer gives each one its default."""
    mesh = _bush_mesh()
    for member in ("K", "B", "GE", "SA", "ST", "EA", "ET"):
        del mesh.field_data[f"nastran:PBUSH_7:{member}"]
    return mesh


def test_a_pbush_without_values_writes_the_defaults(engine, tmp_path):
    path = tmp_path / "defaults.op2"
    write_with(engine, path, _default_bush_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        np.testing.assert_allclose(back.field_data["nastran:PBUSH_7:K"], [1.0] * 6)
        np.testing.assert_allclose(back.field_data["nastran:PBUSH_7:B"], [0.0] * 6)
        np.testing.assert_allclose(back.field_data["nastran:PBUSH_7:GE"], [0.0] * 6)
        for member in ("SA", "ST", "EA", "ET"):
            assert back.field_data[f"nastran:PBUSH_7:{member}"][0] == pytest.approx(1.0)


def test_the_engines_write_the_same_bytes_with_default_bushes(tmp_path):
    mesh = _default_bush_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


# -- CVISC viscous dampers and PVISC properties -------------------------------------------


def _visc_mesh(values=True):
    """Two CVISCs of PVISC 8 between GRIDs 1-2 and 2-3."""
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    mesh = meshioplusplus.Mesh(points, [("line", np.array([[0, 1], [1, 2]]))])
    mesh.cell_data["nastran:eid"] = [np.array([60, 61])]
    mesh.cell_data["nastran:pid"] = [np.array([8, 8])]
    mesh.cell_data["nastran:card"] = [np.array([39, 39])]
    if values:
        mesh.field_data["nastran:PVISC_8:CE"] = np.array([2.5])
        mesh.field_data["nastran:PVISC_8:CR"] = np.array([0.75])
    return mesh


def test_cviscs_and_their_pviscs_round_trip(engine, tmp_path):
    path = tmp_path / "cvisc.op2"
    write_with(engine, path, _visc_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        assert back.cells[0].type == "line"
        np.testing.assert_array_equal(back.cells[0].data, [[0, 1], [1, 2]])
        np.testing.assert_array_equal(back.cell_data["nastran:eid"][0], [60, 61])
        np.testing.assert_array_equal(back.cell_data["nastran:pid"][0], [8, 8])
        np.testing.assert_array_equal(back.cell_data["nastran:card"][0], [39, 39])
        np.testing.assert_allclose(back.field_data["nastran:PVISC_8:CE"], [2.5])
        np.testing.assert_allclose(back.field_data["nastran:PVISC_8:CR"], [0.75])


def test_a_pvisc_without_values_writes_the_defaults(engine, tmp_path):
    path = tmp_path / "defaults.op2"
    write_with(engine, path, _visc_mesh(values=False))
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        np.testing.assert_allclose(back.field_data["nastran:PVISC_8:CE"], [1.0])
        np.testing.assert_allclose(back.field_data["nastran:PVISC_8:CR"], [0.0])


@pytest.mark.parametrize("values", [True, False])
def test_the_engines_write_the_same_bytes_with_viscous_dampers(tmp_path, values):
    mesh = _visc_mesh(values)
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_pynastran_reads_the_cviscs_and_pviscs(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_cvisc.op2"
    _core.nastran_op2_write(str(path), _visc_mesh())
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    first, second = model.elements[60], model.elements[61]
    assert first.type == second.type == "CVISC" and first.pid == 8
    assert first.node_ids == [1, 2] and second.node_ids == [2, 3]
    prop = model.properties[8]
    assert prop.type == "PVISC"
    assert (prop.ce, prop.cr) == pytest.approx((2.5, 0.75))


# -- PLOTEL plot elements ---------------------------------------------------------------------


def _plotel_mesh():
    """Two PLOTELs between GRIDs 1-2 and 2-3 (they have no property)."""
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    mesh = meshioplusplus.Mesh(points, [("line", np.array([[0, 1], [1, 2]]))])
    mesh.cell_data["nastran:eid"] = [np.array([70, 71])]
    mesh.cell_data["nastran:card"] = [np.array([52, 52])]
    return mesh


def test_plotels_round_trip(engine, tmp_path):
    path = tmp_path / "plotel.op2"
    write_with(engine, path, _plotel_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        assert back.cells[0].type == "line"
        np.testing.assert_array_equal(back.cells[0].data, [[0, 1], [1, 2]])
        np.testing.assert_array_equal(back.cell_data["nastran:eid"][0], [70, 71])
        np.testing.assert_array_equal(back.cell_data["nastran:pid"][0], [-1, -1])
        np.testing.assert_array_equal(back.cell_data["nastran:card"][0], [52, 52])


def test_the_engines_write_the_same_bytes_with_plotels(tmp_path):
    mesh = _plotel_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_plotel_cannot_sit_in_a_vertex_block(engine, tmp_path):
    mesh = _plotel_mesh()
    mesh.cells[0] = meshioplusplus.CellBlock("vertex", np.array([[0], [1]]))
    with pytest.raises(Exception, match="cannot hold PLOTEL"):
        write_with(engine, tmp_path / "bad.op2", mesh)


def test_pynastran_reads_the_plotels(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_plotel.op2"
    _core.nastran_op2_write(str(path), _plotel_mesh())
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    assert sorted(model.plotels) == [70, 71]
    assert model.plotels[70].node_ids == [1, 2]
    assert model.plotels[71].node_ids == [2, 3]


# -- CSHEAR shear panels and PSHEAR properties -----------------------------------------------


def _shear_mesh(values=True):
    """A CQUAD4 (PSHELL 1) and a CSHEAR (PSHEAR 4) in one quad block, GRIDs 1 to 6."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
        + [[2.0, 0.0, 0.0], [2.0, 1.0, 0.0]]
    )
    mesh = meshioplusplus.Mesh(
        points, [("quad", np.array([[0, 1, 2, 3], [1, 4, 5, 2]]))]
    )
    mesh.cell_data["nastran:eid"] = [np.array([80, 81])]
    mesh.cell_data["nastran:pid"] = [np.array([1, 4])]
    mesh.cell_data["nastran:card"] = [np.array([51, 31])]
    if values:
        mesh.field_data["nastran:PSHEAR_4:MID"] = np.array([3], dtype=np.int64)
        mesh.field_data["nastran:PSHEAR_4:T"] = np.array([0.25])
        mesh.field_data["nastran:PSHEAR_4:NSM"] = np.array([0.5])
        mesh.field_data["nastran:PSHEAR_4:F1"] = np.array([0.75])
        mesh.field_data["nastran:PSHEAR_4:F2"] = np.array([1.5])
    return mesh


def _by_eid(mesh):
    """``{eid: (nodes, pid, card)}`` of a mesh read from an OP2."""
    out = {}
    for block, eids, pids, cards in zip(
        mesh.cells,
        mesh.cell_data["nastran:eid"],
        mesh.cell_data["nastran:pid"],
        mesh.cell_data["nastran:card"],
    ):
        for row, eid in enumerate(eids):
            out[int(eid)] = (block.data[row].tolist(), int(pids[row]), int(cards[row]))
    return out


def test_cshears_and_pshears_round_trip(engine, tmp_path):
    path = tmp_path / "shear.op2"
    write_with(engine, path, _shear_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        assert _by_eid(back) == {
            80: ([0, 1, 2, 3], 1, 51),
            81: ([1, 4, 5, 2], 4, 31),
        }
        assert back.field_data["nastran:PSHEAR_4:MID"][0] == 3
        for member, value in (("T", 0.25), ("NSM", 0.5), ("F1", 0.75), ("F2", 1.5)):
            np.testing.assert_allclose(
                back.field_data[f"nastran:PSHEAR_4:{member}"], [value]
            )
        # the PSHEAR's material is written as a MAT1
        assert "nastran:MAT1_3:E" in back.field_data


def test_a_pshear_without_values_writes_the_defaults(engine, tmp_path):
    path = tmp_path / "defaults.op2"
    write_with(engine, path, _shear_mesh(values=False))
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        assert back.field_data["nastran:PSHEAR_4:MID"][0] == 1
        np.testing.assert_allclose(back.field_data["nastran:PSHEAR_4:T"], [1.0])
        np.testing.assert_allclose(back.field_data["nastran:PSHEAR_4:F1"], [0.0])


@pytest.mark.parametrize("values", [True, False])
def test_the_engines_write_the_same_bytes_with_shear_panels(tmp_path, values):
    mesh = _shear_mesh(values)
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_pshear_cannot_share_an_id_with_a_pshell(engine, tmp_path):
    mesh = _shear_mesh()
    mesh.cell_data["nastran:pid"] = [np.array([4, 4])]
    with pytest.raises(Exception, match="PSHELL and PSHEAR|PSHEAR and PSHELL"):
        write_with(engine, tmp_path / "shared.op2", mesh)


def test_pynastran_reads_the_cshears_and_pshears(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_shear.op2"
    _core.nastran_op2_write(str(path), _shear_mesh())
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    assert model.elements[80].type == "CQUAD4"
    shear = model.elements[81]
    assert shear.type == "CSHEAR" and shear.pid == 4
    assert shear.node_ids == [2, 5, 6, 3]
    prop = model.properties[4]
    assert prop.type == "PSHEAR" and prop.mid == 3
    assert (prop.t, prop.nsm, prop.f1, prop.f2) == pytest.approx((0.25, 0.5, 0.75, 1.5))


# -- MAT2 and MAT8 materials -----------------------------------------------------------------

_MAT2_5 = {
    "G11": 4.0e4,
    "G12": 1.0e4,
    "G13": 0.5e4,
    "G22": 3.0e4,
    "G23": 0.25e4,
    "G33": 2.0e4,
    "RHO": 7.8e-9,
    "A1": 1.0e-5,
    "A2": 2.0e-5,
    "A3": 3.0e-5,
    "TREF": 20.0,
    "GE": 0.02,
    "ST": 100.0,
    "SC": 200.0,
    "SS": 50.0,
}
_MAT8_2 = {
    "E1": 1.5e5,
    "E2": 9.0e3,
    "NU12": 0.28,
    "G12": 5.0e3,
    "G1Z": 4.0e3,
    "G2Z": 3.0e3,
    "RHO": 1.6e-9,
    "A1": 1.0e-6,
    "A2": 3.0e-5,
    "TREF": 25.0,
    "XT": 1500.0,
    "XC": 1200.0,
    "YT": 50.0,
    "YC": 250.0,
    "S": 70.0,
    "GE": 0.01,
    "F12": -0.5,
    "STRN": 1.0,
}


def _anisotropic_mesh(mat1=True, mat2=True, mat8=True):
    """A PCOMP quad (CQUAD4 7, PCOMP 9, plies of materials 1, 2 and 1) and a PSHELL
    triangle (CTRIA3 8, PSHELL 1 on material 5), with material 2 a MAT8 and material 5 a
    MAT2 when asked, and material 1 a MAT1 of its own."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
        + [[2.0, 0.0, 0.0]]
    )
    mesh = meshioplusplus.Mesh(
        points,
        [("quad", np.array([[0, 1, 2, 3]])), ("triangle", np.array([[1, 4, 2]]))],
    )
    mesh.cell_data["nastran:eid"] = [np.array([7]), np.array([8])]
    mesh.cell_data["nastran:pid"] = [np.array([9]), np.array([1])]
    mesh.field_data["nastran:PCOMP_9:PLY"] = np.array(
        [[1, 0.1, 0.0, 1], [2, 0.2, 45.0, 0], [1, 0.1, -45.0, 1]], dtype=float
    )
    mesh.field_data["nastran:PSHELL_1:MID1"] = np.array([5], dtype=np.int64)
    if mat1:
        mesh.field_data["nastran:MAT1_1:E"] = np.array([70000.0])
    if mat2:
        for member, value in _MAT2_5.items():
            mesh.field_data[f"nastran:MAT2_5:{member}"] = np.array([value])
        mesh.field_data["nastran:MAT2_5:MCSID"] = np.array([3], dtype=np.int64)
    if mat8:
        for member, value in _MAT8_2.items():
            mesh.field_data[f"nastran:MAT8_2:{member}"] = np.array([value])
    return mesh


def test_mat2_and_mat8_materials_round_trip(engine, tmp_path):
    path = tmp_path / "anisotropic.op2"
    write_with(engine, path, _anisotropic_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        for member, value in _MAT2_5.items():
            np.testing.assert_allclose(
                back.field_data[f"nastran:MAT2_5:{member}"], [value], rtol=1e-6
            )
        assert back.field_data["nastran:MAT2_5:MCSID"].tolist() == [3]
        for member, value in _MAT8_2.items():
            np.testing.assert_allclose(
                back.field_data[f"nastran:MAT8_2:{member}"], [value], rtol=1e-6
            )
        # a material is the one card it was given as: no MAT1 twin of the others
        assert "nastran:MAT1_1:E" in back.field_data
        assert "nastran:MAT1_2:E" not in back.field_data
        assert "nastran:MAT1_5:E" not in back.field_data


def test_materials_without_values_stay_mat1(engine, tmp_path):
    path = tmp_path / "plain.op2"
    write_with(engine, path, _anisotropic_mesh(mat1=False, mat2=False, mat8=False))
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        for mid in (1, 2, 5):
            assert back.field_data[f"nastran:MAT1_{mid}:E"][0] == 1.0
        assert not any(":MAT2_" in name or ":MAT8_" in name for name in back.field_data)


@pytest.mark.parametrize(
    "flags", [(True, True, True), (False, True, False), (False, False, True)]
)
def test_the_engines_write_the_same_bytes_with_mat2_and_mat8(tmp_path, flags):
    mesh = _anisotropic_mesh(*flags)
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_material_id_given_as_two_cards_is_refused(engine, tmp_path):
    mesh = _anisotropic_mesh()
    mesh.field_data["nastran:MAT1_2:E"] = np.array([1.0e5])  # 2 is also a MAT8
    with pytest.raises(Exception, match="material id 2 is given as both MAT1 and MAT8"):
        write_with(engine, tmp_path / "conflict.op2", mesh)


def test_pynastran_reads_the_mat2_and_mat8(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_anisotropic.op2"
    _core.nastran_op2_write(str(path), _anisotropic_mesh())
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    assert sorted(model.materials) == [1, 2, 5]
    ply = model.materials[2]
    assert ply.type == "MAT8"
    assert (ply.e11, ply.e22, ply.nu12) == pytest.approx((1.5e5, 9.0e3, 0.28))
    assert (ply.g12, ply.Xt, ply.F12) == pytest.approx((5.0e3, 1500.0, -0.5))
    shell = model.materials[5]
    assert shell.type == "MAT2"
    assert (shell.G11, shell.G12, shell.G22, shell.G33) == pytest.approx(
        (4.0e4, 1.0e4, 3.0e4, 2.0e4)
    )
    assert shell.mcsid == 3
    assert model.materials[1].type == "MAT1"
    assert model.properties[1].mid1 == 5


# -- CGAP gaps and PGAP properties ------------------------------------------------------------


def _gap_mesh(values=True):
    """Two CGAPs of PGAP 12 between GRIDs 1-2 and 2-3: the first oriented by X = (0, 0, 1),
    the second by GRID 4 as its G0."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    )
    mesh = meshioplusplus.Mesh(points, [("line", np.array([[0, 1], [1, 2]]))])
    mesh.cell_data["nastran:eid"] = [np.array([80, 81])]
    mesh.cell_data["nastran:pid"] = [np.array([12, 12])]
    mesh.cell_data["nastran:card"] = [np.array([19, 19])]
    mesh.cell_data["nastran:x"] = [np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0]])]
    mesh.cell_data["nastran:g0"] = [np.array([0, 4])]
    if values:
        for member, value in (
            ("U0", 0.01),
            ("F0", 2.0),
            ("KA", 5.0e6),
            ("KB", 0.5),
            ("KT", 1.5e3),
            ("MU1", 0.3),
            ("MU2", 0.2),
            ("TMAX", 4.0),
            ("MAR", 50.0),
            ("TRMIN", 0.01),
        ):
            mesh.field_data[f"nastran:PGAP_12:{member}"] = np.array([value])
    return mesh


def test_cgaps_and_their_pgaps_round_trip(engine, tmp_path):
    path = tmp_path / "cgap.op2"
    write_with(engine, path, _gap_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        assert back.cells[0].type == "line"
        np.testing.assert_array_equal(back.cells[0].data, [[0, 1], [1, 2]])
        np.testing.assert_array_equal(back.cell_data["nastran:eid"][0], [80, 81])
        np.testing.assert_array_equal(back.cell_data["nastran:pid"][0], [12, 12])
        np.testing.assert_array_equal(back.cell_data["nastran:card"][0], [19, 19])
        np.testing.assert_array_equal(back.cell_data["nastran:g0"][0], [0, 4])
        np.testing.assert_allclose(
            back.cell_data["nastran:x"][0], [[0, 0, 1], [0, 0, 0]]
        )
        assert [r.name for r in back.regions] == ["PGAP_12"]
        for member, value in (("U0", 0.01), ("KA", 5.0e6), ("MU1", 0.3), ("MAR", 50.0)):
            np.testing.assert_allclose(
                back.field_data[f"nastran:PGAP_12:{member}"], [value], rtol=1e-6
            )


def test_a_pgap_without_values_writes_the_defaults(engine, tmp_path):
    path = tmp_path / "defaults.op2"
    write_with(engine, path, _gap_mesh(values=False))
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        np.testing.assert_allclose(back.field_data["nastran:PGAP_12:KA"], [1.0e8])
        np.testing.assert_allclose(back.field_data["nastran:PGAP_12:MAR"], [100.0])


@pytest.mark.parametrize("values", [True, False])
def test_the_engines_write_the_same_bytes_with_gaps(tmp_path, values):
    mesh = _gap_mesh(values)
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_cgap_needs_its_orientation(engine, tmp_path):
    mesh = _gap_mesh()
    for name in ("nastran:x", "nastran:g0"):
        del mesh.cell_data[name]
    with pytest.raises(Exception, match="orientation"):
        write_with(engine, tmp_path / "x.op2", mesh)


def test_pynastran_reads_the_cgaps_and_pgaps(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_cgap.op2"
    _core.nastran_op2_write(str(path), _gap_mesh())
    model = op2.OP2Geom(make_geom=True, debug=False)
    model.read_op2(str(path))
    first, second = model.elements[80], model.elements[81]
    assert first.type == second.type == "CGAP" and first.pid == 12
    assert first.node_ids == [1, 2] and second.node_ids == [2, 3]
    np.testing.assert_allclose(first.x, [0, 0, 1])
    assert first.g0 is None and second.g0 == 4
    prop = model.properties[12]
    assert prop.type == "PGAP"
    assert (prop.u0, prop.f0, prop.ka, prop.kb, prop.kt) == pytest.approx(
        (0.01, 2.0, 5.0e6, 0.5, 1.5e3), rel=1e-6
    )
    assert (prop.mu1, prop.mu2, prop.tmax, prop.mar, prop.trmin) == pytest.approx(
        (0.3, 0.2, 4.0, 50.0, 0.01), rel=1e-6
    )


# -- element results of every card, in every layout the reader decodes --------------------
#
# The writer inverts the reader's layout table, so for every card, element type, result
# group and complex/random variant: fill the arrays of the first layout the reader has for
# it, write through both engines (the same bytes), and read every value back.


def _rod_mesh(code):
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
        [("line", np.array([[0, 1], [1, 2]]))],
    )
    mesh.cell_data["nastran:eid"] = [np.array([5, 6])]
    mesh.cell_data["nastran:pid"] = [np.array([4, 4])]
    mesh.cell_data["nastran:card"] = [np.array([code, code])]
    return mesh


_RESULT_MESHES = {
    "CROD": lambda: _rod_mesh(30),
    "CTUBE": lambda: _rod_mesh(37),
    "CONROD": lambda: _rod_mesh(16),
    "CQUAD4": _shear_mesh,
    "CTRIA3": every_card,
    "CQUAD8": every_card,
    "CTRIA6": every_card,
    "CTETRA": every_card,
    "CHEXA": every_card,
    "CPENTA": every_card,
    "CPYRAM": every_card,
    "CBAR": _bar_mesh,
    "CBEAM": _beam_mesh,
    "CBUSH": _bush_mesh,
    "CSHEAR": _shear_mesh,
    "CELAS1": _spring_mesh,
    "CDAMP1": _spring_mesh,
    "CELAS2": _inline_spring_mesh,
    "CDAMP2": _inline_spring_mesh,
    "CVISC": _visc_mesh,
}


def _result_cases():
    cases = []
    for group, family, _base, _table, _thermal in py_write._RESULT_GROUPS:
        variants = [("", False, False), ("", True, False)]
        if family in ("stress", "force"):
            variants += [(s, False, True) for s in py_write._RANDOM_KINDS]
        for suffix, complex_, random in variants:
            for card, alternatives in py_write._RESULT_CARDS:
                for etype, strict in alternatives:
                    layouts = py_write._candidates(
                        family, etype, complex_, random and family == "stress"
                    )
                    if layouts:
                        label = f"{group}-{card}-{etype}-{suffix or 'c' * complex_ or 'real'}"
                        cases.append(
                            pytest.param(
                                (group, family, card, etype, strict),
                                (suffix, complex_, random),
                                id=label,
                            )
                        )
    return cases


def _result_mesh(card_info, variant):
    """A mesh whose ``card`` elements carry random values for every member of the
    first layout of the element type; the values by element id, per array name."""
    group, family, card, etype, _strict = card_info
    suffix, complex_, random = variant
    mesh = _RESULT_MESHES[card]()
    model = py_write._Model(mesh)
    eid, _pid, parts, order = model.element_rows(card)
    lay, _wide, _s_code = py_write._candidates(
        family, etype, complex_, random and family == "stress"
    )[0]
    loc = py_write._location(lay)
    rng = np.random.default_rng(etype)
    if lay["kind"] == "blocks" and lay["stations"]:
        width = lay["blocks"]
    elif lay["kind"] == "blocks":
        width = lay["blocks"] - 1
    else:
        width = 3
    shapes = {}
    for member in lay["members"]:
        for name, _word in py_write._names_of(member, suffix):
            if lay["kind"] == "row" or (
                lay["kind"] == "blocks" and not lay["stations"]
            ):
                shapes[name] = ()
            if loc is not None:
                shapes[name + loc] = (width,)
    # the values by element id, laid into the blocks' cell data
    positions = {}
    at = 0
    for blk, rows in parts:
        positions[blk["index"]] = (np.asarray(rows), at)
        at += len(rows)
    sorted_at = np.empty(len(eid), dtype=np.int64)
    sorted_at[order] = np.arange(len(eid))
    expected = {}
    for key, shape in shapes.items():
        values = rng.random((len(eid),) + shape) + 1.0
        values = values.astype(np.float32).astype(np.float64)
        expected[key] = dict(zip(eid.tolist(), values))
        blocks = []
        for bi, block in enumerate(mesh.cells):
            array = np.full((len(block.data),) + shape, np.nan)
            if bi in positions:
                rows, start = positions[bi]
                for j, r in enumerate(rows.tolist()):
                    array[r] = values[sorted_at[start + j]]
            blocks.append(array)
        mesh.cell_data[f"{group}:{key}"] = blocks
    return mesh, expected


@pytest.mark.parametrize("card_info, variant", _result_cases())
def test_element_results_of_every_layout(card_info, variant, tmp_path):
    group = card_info[0]
    mesh, expected = _result_mesh(card_info, variant)
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(core)
        got = {}
        for bi, eids in enumerate(back.cell_data["nastran:eid"]):
            for name in expected:
                array = back.cell_data.get(f"{group}:{name}")
                if array is None:
                    continue
                for row, eid in enumerate(eids.tolist()):
                    value = np.asarray(array[bi][row])
                    if not np.isnan(value).all():
                        got.setdefault(name, {})[eid] = value
        for name, wanted in expected.items():
            assert name in got, f"{name} was not read back"
            for eid, value in wanted.items():
                assert eid in got[name], f"{name} of element {eid} was not read back"
                width = np.shape(value)[-1] if np.ndim(value) else None
                have = got[name][eid]
                if width is not None:
                    have = have[:width]
                np.testing.assert_allclose(have, value, rtol=1e-6)


@pytest.mark.parametrize(
    "card_info, variant, where, shape",
    [
        (("STRESS", "stress", "CROD", 1, False), ("", True, False), "stress.crod_stress", (1, 2, 2)),
        (("STRESS", "stress", "CROD", 1, False), ("_RMS", False, True), "rms.crod_stress", (1, 2, 2)),
        (("STRESS", "stress", "CROD", 1, False), ("_NO", False, True), "no.crod_stress", (1, 2, 2)),
        (("ELEMENT_FORCE", "force", "CROD", 1, False), ("_RMS", False, True), "rms.crod_force", (1, 2, 2)),
        (("STRESS", "stress", "CBAR", 34, False), ("", False, False), "stress.cbar_stress", (1, 2, 15)),
        (("STRESS", "stress", "CBAR", 100, True), ("", False, False), "stress.cbar_stress_10nodes", (1, 6, 9)),
        (("STRESS", "stress", "CQUAD4", 95, True), ("", False, False), "stress.cquad4_composite_stress", (1, 3, 9)),
        (("STRESS", "stress", "CQUAD4", 144, True), ("", False, False), "stress.cquad4_stress", (1, 10, 8)),
        (("STRESS", "stress", "CQUAD8", 96, True), ("", False, False), "stress.cquad8_composite_stress", (1, 3, 9)),
        (("STRESS", "stress", "CBEAM", 2, False), ("", False, False), "stress.cbeam_stress", (1, 22, 8)),
        (("ELEMENT_FORCE", "force", "CBUSH", 102, False), ("", False, False), "force.cbush_force", (1, 2, 6)),
        (("ELEMENT_FORCE", "force", "CSHEAR", 4, False), ("", False, False), "force.cshear_force", (1, 1, 16)),
        (("ELEMENT_FORCE", "force", "CELAS1", 11, False), ("", False, False), "force.celas1_force", (1, 3, 1)),
    ],
)  # fmt: skip
def test_pynastran_reads_the_element_results(
    card_info, variant, where, shape, tmp_path
):
    op2 = pytest.importorskip("pyNastran.op2.op2", exc_type=ImportError)
    mesh, _expected = _result_mesh(card_info, variant)
    if variant[1] or variant[2]:  # complex and random tables are a frequency response's
        mesh.field_data["nastran:analysis"] = np.array([5])
        mesh.field_data["meshio:time"] = np.array([1.0])
    path = tmp_path / "pn_results.op2"
    _core.nastran_op2_write(str(path), mesh)
    model = op2.OP2(debug=False)
    model.read_op2(str(path), build_dataframe=False)
    group, name = where.split(".")
    found = getattr(getattr(model.op2_results, group), name)
    assert {k: v.data.shape for k, v in found.items()} == {1: shape}
    assert np.isfinite(np.abs(next(iter(found.values())).data)).all()


def test_the_engines_write_the_same_bytes_with_every_kind_of_result(tmp_path):
    mesh = _rod_mesh(30)
    for name in (
        "STRESS:A",
        "STRESS:A_real",
        "STRESS:A_imag",
        "STRESS:A_RMS",
        "STRESS:A_PSD",
        "STRAIN:A",
        "ELEMENT_FORCE:AF",
        "ELEMENT_FORCE:AF_NO",
        "ELEMENT_FORCE:AF_real",
        "ELEMENT_FORCE:AF_imag",
        "HEAT_FLUX:XFLUX",
    ):
        mesh.cell_data[name] = [np.array([1.5, 2.5]) + len(name)]
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()
    back = meshioplusplus.nastran_op2.read(core)
    assert "STRESS:A_PSD" in back.cell_data and "HEAT_FLUX:XFLUX" in back.cell_data


def _energy_and_grid_force_mesh():
    """Two CROD/CQUAD4 cells with strain and kinetic energies, grid point forces at
    their nodes and the totals at every point."""
    mesh = _shear_mesh()  # a CQUAD4 and a CSHEAR
    mesh.point_data["DISPLACEMENT"] = np.zeros((len(mesh.points), 3))
    mesh.cell_data["ENERGY:ENERGY"] = [np.array([1.5, 2.5])]
    mesh.cell_data["ENERGY:PCT"] = [np.array([37.5, 62.5])]
    mesh.cell_data["KINETIC_ENERGY:ENERGY"] = [np.array([0.5, 0.25])]
    for k, member in enumerate(("F1", "F2", "F3", "M1", "M2", "M3")):
        forces = np.arange(8, dtype=float).reshape(2, 4) + k
        forces[1, 2] = np.nan  # no force at one node
        mesh.cell_data["GRID_FORCE:" + member] = [forces]
        mesh.point_data[f"GRID_FORCE:TOTALS:{member}"] = (
            np.arange(len(mesh.points), dtype=float) + k
        )
    return mesh


def test_energies_and_grid_point_forces_round_trip(engine, tmp_path):
    mesh = _energy_and_grid_force_mesh()
    path = tmp_path / "energy.op2"
    write_with(engine, path, mesh)
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        eids = {}
        for block, ids in enumerate(back.cell_data["nastran:eid"]):
            for row, eid in enumerate(ids.tolist()):
                eids[eid] = (block, row)
        energy = {
            eid: back.cell_data["ENERGY:ENERGY"][b][r] for eid, (b, r) in eids.items()
        }
        assert sorted(energy.values()) == [1.5, 2.5]
        kinetic = sorted(
            back.cell_data["KINETIC_ENERGY:ENERGY"][b][r] for b, r in eids.values()
        )
        assert kinetic == [0.25, 0.5]
        total = sum(
            np.nansum(back.cell_data["GRID_FORCE:F1"][b][r]) for b, r in eids.values()
        )
        assert total == pytest.approx(
            0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 - 6
        )  # one node unset
        np.testing.assert_allclose(
            back.point_data["GRID_FORCE:TOTALS:M3"].ravel(),
            np.arange(len(mesh.points)) + 5.0,
        )


def test_the_engines_write_the_same_bytes_with_energies_and_grid_forces(tmp_path):
    mesh = _energy_and_grid_force_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


# -- thermal and anisotropic solid materials ----------------------------------------------

_MAT4_6 = {
    "K": 45.0,
    "CP": 460.0,
    "RHO": 7800.0,
    "H": 12.0,
    "MU": 0.5,
    "HGEN": 2.0,
    "REFENTH": 3.0,
    "TCH": 4.0,
    "TDELTA": 5.0,
    "QLAT": 6.0,
}
_MAT5_7 = {
    "K1": 10.0,
    "K2": 1.0,
    "K3": 2.0,
    "K4": 20.0,
    "K5": 3.0,
    "K6": 30.0,
    "CP": 400.0,
    "RHO": 7000.0,
    "HGEN": 1.5,
}
_MAT9_8 = {
    f"G{i}{j}": 1000.0 * (10 * i + j) for i in range(1, 7) for j in range(i, 7)
} | {
    "RHO": 2.5,
    "A1": 1e-6,
    "A2": 2e-6,
    "A3": 3e-6,
    "A4": 4e-6,
    "A5": 5e-6,
    "A6": 6e-6,
    "TREF": 20.0,
    "GE": 0.01,
}


def _thermal_mesh():
    """A tetra whose property 2 is on MAT9 8, plus a MAT4 6 and a MAT5 7."""
    mesh = oracle_mesh()
    mesh.field_data["nastran:PSHELL_7:MID1"] = np.array([8], dtype=np.int64)
    for card, ident, values in (
        ("MAT4", 6, _MAT4_6),
        ("MAT5", 7, _MAT5_7),
        ("MAT9", 8, _MAT9_8),
    ):
        for member, value in values.items():
            mesh.field_data[f"nastran:{card}_{ident}:{member}"] = np.array([value])
    return mesh


def test_thermal_and_solid_materials_round_trip(engine, tmp_path):
    path = tmp_path / "thermal.op2"
    write_with(engine, path, _thermal_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        for card, ident, values in (
            ("MAT4", 6, _MAT4_6),
            ("MAT5", 7, _MAT5_7),
            ("MAT9", 8, _MAT9_8),
        ):
            for member, value in values.items():
                np.testing.assert_allclose(
                    back.field_data[f"nastran:{card}_{ident}:{member}"],
                    [value],
                    rtol=1e-6,
                )
        assert "nastran:MAT1_8:E" not in back.field_data  # MAT9 8 is not also a MAT1


def test_the_engines_write_the_same_bytes_with_thermal_materials(tmp_path):
    mesh = _thermal_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_pynastran_reads_the_thermal_and_solid_materials(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_thermal.op2"
    _core.nastran_op2_write(str(path), _thermal_mesh())
    model = op2.read_op2_geom(str(path), debug=False)
    assert model.thermal_materials[6].type == "MAT4"
    assert model.thermal_materials[6].k == pytest.approx(45.0)
    assert model.thermal_materials[7].type == "MAT5"
    assert model.thermal_materials[7].kyy == pytest.approx(20.0)
    mat9 = model.materials[8]
    assert mat9.type == "MAT9"
    assert mat9.rho == pytest.approx(2.5)
    assert mat9.G11 == pytest.approx(11000.0) and mat9.G66 == pytest.approx(66000.0)


# -- point masses: CMASS1 with its PMASS, CMASS2 and CONM1 --------------------------------

_CONM1_TERMS = np.arange(1.0, 22.0)


def _mass_mesh():
    """A CMASS1 (PMASS 5) between GRIDs 1 and 2 and one grounded at GRID 3, a CMASS2
    between GRIDs 2 and 3 and a CONM1 at GRID 4."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
    )
    mesh = meshioplusplus.Mesh(
        points,
        [
            ("line", np.array([[0, 1]])),
            ("vertex", np.array([[2]])),
            ("line", np.array([[1, 2]])),
            ("vertex", np.array([[3]])),
        ],
    )
    mesh.cell_data["nastran:eid"] = [
        np.array([20]),
        np.array([21]),
        np.array([22]),
        np.array([23]),
    ]
    mesh.cell_data["nastran:pid"] = [
        np.array([5]),
        np.array([5]),
        np.array([-1]),
        np.array([-1]),
    ]
    mesh.cell_data["nastran:card"] = [
        np.array([10]),
        np.array([10]),
        np.array([11]),
        np.array([14]),
    ]
    mesh.cell_data["nastran:c1"] = [
        np.array([1]),
        np.array([2]),
        np.array([3]),
        np.array([0]),
    ]
    mesh.cell_data["nastran:c2"] = [
        np.array([1]),
        np.array([0]),
        np.array([3]),
        np.array([0]),
    ]
    mesh.cell_data["nastran:m"] = [
        np.array([0.0]),
        np.array([0.0]),
        np.array([2.5]),
        np.array([0.0]),
    ]
    mesh.cell_data["nastran:conm1"] = [
        np.zeros((1, 21)),
        np.zeros((1, 21)),
        np.zeros((1, 21)),
        _CONM1_TERMS[None, :],
    ]
    mesh.field_data["nastran:PMASS_5:M"] = np.array([0.75])
    return mesh


def test_point_masses_round_trip(engine, tmp_path):
    path = tmp_path / "masses.op2"
    write_with(engine, path, _mass_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        by_eid = {}
        for block, eids in zip(back.cells, back.cell_data["nastran:eid"]):
            for row, eid in enumerate(eids.tolist()):
                by_eid[eid] = (block.type, block.data[row].tolist())
        assert by_eid[20] == ("line", [0, 1])
        assert by_eid[21] == ("vertex", [2])
        assert by_eid[22] == ("line", [1, 2])
        assert by_eid[23] == ("vertex", [3])
        cards = {}
        for ids, card in zip(
            back.cell_data["nastran:eid"], back.cell_data["nastran:card"]
        ):
            cards.update(zip(ids.tolist(), card.tolist()))
        assert cards == {20: 10, 21: 10, 22: 11, 23: 14}
        m = {}
        for ids, values in zip(
            back.cell_data["nastran:eid"], back.cell_data["nastran:m"]
        ):
            m.update(zip(ids.tolist(), values.tolist()))
        assert m[22] == 2.5
        conm1 = {}
        for ids, values in zip(
            back.cell_data["nastran:eid"], back.cell_data["nastran:conm1"]
        ):
            conm1.update(zip(ids.tolist(), values))
        np.testing.assert_allclose(conm1[23], _CONM1_TERMS)
        assert back.field_data["nastran:PMASS_5:M"][0] == pytest.approx(0.75)
        assert [r.name for r in back.regions] == ["PMASS_5"]


def test_the_engines_write_the_same_bytes_with_point_masses(tmp_path):
    mesh = _mass_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_cmass2_and_a_conm1_need_their_values(engine, tmp_path):
    for name, card in (("nastran:m", "CMASS2"), ("nastran:conm1", "CONM1")):
        mesh = _mass_mesh()
        del mesh.cell_data[name]
        with pytest.raises(Exception, match=f"a {card} needs its value"):
            write_with(engine, tmp_path / "x.op2", mesh)


def test_pynastran_reads_the_point_masses(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_masses.op2"
    _core.nastran_op2_write(str(path), _mass_mesh())
    model = op2.read_op2_geom(str(path), debug=False)
    assert model.masses[20].type == "CMASS1" and model.masses[20].nodes == [1, 2]
    assert model.masses[20].pid == 5 and model.masses[20].c1 == 1
    assert model.masses[21].nodes == [3, None]
    assert model.masses[22].type == "CMASS2" and model.masses[22].mass == 2.5
    assert model.masses[22].nodes == [2, 3] and model.masses[22].c1 == 3
    conm1 = model.masses[23]
    assert conm1.type == "CONM1" and conm1.nid == 4
    assert conm1.mass_matrix[0, 0] == 1.0 and conm1.mass_matrix[5, 5] == 21.0
    assert model.properties_mass[5].mass == pytest.approx(0.75)


# -- PBARL and PBEAML library sections, PCOMPG composites ---------------------------------


def _text(value):
    return np.frombuffer(value.encode(), dtype=np.uint8)


_BEAML_DIMS = np.array(
    [[10.0, 4.0, 2.0, 1.0, 1.5, 0.5], [10.0, 4.0, 2.5, 1.2, 1.5, 0.5]]
)


def _section_mesh():
    """A CBAR on a PBARL (11), a CBEAM on a PBEAML (12) and a quad on a PCOMPG (13)."""
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ]
    )
    mesh = meshioplusplus.Mesh(
        points,
        [
            ("line", np.array([[0, 1]])),
            ("line", np.array([[1, 2]])),
            ("quad", np.array([[0, 1, 4, 3]])),
        ],
    )
    mesh.cell_data["nastran:eid"] = [np.array([40]), np.array([41]), np.array([42])]
    mesh.cell_data["nastran:pid"] = [np.array([11]), np.array([12]), np.array([13])]
    mesh.cell_data["nastran:card"] = [np.array([24]), np.array([54]), np.array([51])]
    mesh.cell_data["nastran:x"] = [
        np.array([[0.0, 0.0, 1.0]]),
        np.array([[0.0, 0.0, 1.0]]),
        np.zeros((1, 3)),
    ]
    fields = mesh.field_data
    fields["nastran:PBARL_11:MID"] = np.array([3], dtype=np.int64)
    fields["nastran:PBARL_11:GROUP"] = _text("MSCBML0")
    fields["nastran:PBARL_11:TYPE"] = _text("BAR")
    fields["nastran:PBARL_11:DIM"] = np.array([2.0, 3.0])
    fields["nastran:PBARL_11:NSM"] = np.array([0.5])
    fields["nastran:PBEAML_12:MID"] = np.array([3], dtype=np.int64)
    fields["nastran:PBEAML_12:TYPE"] = _text("I")
    fields["nastran:PBEAML_12:SO"] = np.array([0.0, 0.0])
    fields["nastran:PBEAML_12:XXB"] = np.array([0.0, 1.0])
    fields["nastran:PBEAML_12:DIM"] = _BEAML_DIMS
    fields["nastran:PBEAML_12:NSM"] = np.array([0.1, 0.2])
    fields["nastran:PCOMPG_13:LAMOPT"] = np.array([0], dtype=np.int64)
    fields["nastran:PCOMPG_13:Z0"] = np.array([0.5])
    fields["nastran:PCOMPG_13:NSM"] = np.array([0.1])
    fields["nastran:PCOMPG_13:SB"] = np.array([2.0])
    fields["nastran:PCOMPG_13:FT"] = np.array([3], dtype=np.int64)
    fields["nastran:PCOMPG_13:TREF"] = np.array([20.0])
    fields["nastran:PCOMPG_13:GE"] = np.array([0.02])
    fields["nastran:PCOMPG_13:PLY"] = np.array(
        [[1, 4, 0.1, 0.0, 1], [2, 4, 0.2, 45.0, 0], [3, 4, 0.1, -45.0, 1]]
    )
    return mesh


def test_library_sections_and_global_plies_round_trip(engine, tmp_path):
    path = tmp_path / "sections.op2"
    write_with(engine, path, _section_mesh())
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        fields = back.field_data
        assert bytes(fields["nastran:PBARL_11:TYPE"]) == b"BAR"
        assert bytes(fields["nastran:PBARL_11:GROUP"]) == b"MSCBML0"
        assert fields["nastran:PBARL_11:MID"][0] == 3
        np.testing.assert_allclose(fields["nastran:PBARL_11:DIM"], [2.0, 3.0])
        np.testing.assert_allclose(fields["nastran:PBARL_11:NSM"], [0.5])
        assert bytes(fields["nastran:PBEAML_12:TYPE"]) == b"I"
        np.testing.assert_allclose(
            fields["nastran:PBEAML_12:DIM"], _BEAML_DIMS, rtol=1e-6
        )
        np.testing.assert_allclose(fields["nastran:PBEAML_12:XXB"], [0.0, 1.0])
        np.testing.assert_allclose(
            fields["nastran:PBEAML_12:NSM"], [0.1, 0.2], rtol=1e-6
        )
        np.testing.assert_allclose(
            fields["nastran:PCOMPG_13:PLY"],
            np.asarray(_section_mesh().field_data["nastran:PCOMPG_13:PLY"]),
            rtol=1e-6,
        )
        assert fields["nastran:PCOMPG_13:FT"][0] == 3
        assert fields["nastran:PCOMPG_13:Z0"][0] == pytest.approx(0.5)
        assert sorted(r.name for r in back.regions) == [
            "PBARL_11",
            "PBEAML_12",
            "PCOMPG_13",
        ]
        # no PBAR, PBEAM or PCOMP twin, and the plies' material is a MAT1 of its own
        assert not any(
            ":PBAR_" in n or ":PBEAM_" in n or ":PCOMP_" in n for n in fields
        )
        assert "nastran:MAT1_4:E" in fields and "nastran:MAT1_3:E" in fields


def test_the_engines_write_the_same_bytes_with_sections(tmp_path):
    mesh = _section_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_a_section_needs_its_dimensions(engine, tmp_path):
    mesh = _section_mesh()
    mesh.field_data["nastran:PBARL_11:DIM"] = np.array([2.0])
    with pytest.raises(Exception, match="PBARL_11:DIM' must have 2 value"):
        write_with(engine, tmp_path / "x.op2", mesh)
    mesh = _section_mesh()
    mesh.field_data["nastran:PBEAML_12:TYPE"] = _text("NOSUCH")
    with pytest.raises(Exception, match="section type 'NOSUCH'"):
        write_with(engine, tmp_path / "y.op2", mesh)


def test_pynastran_reads_the_sections_and_global_plies(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_sections.op2"
    _core.nastran_op2_write(str(path), _section_mesh())
    model = op2.read_op2_geom(str(path), debug=False)
    bar = model.properties[11]
    assert bar.type == "PBARL" and bar.beam_type == "BAR" and bar.mid == 3
    np.testing.assert_allclose(bar.dim, [2.0, 3.0])
    beam = model.properties[12]
    assert beam.type == "PBEAML" and beam.beam_type == "I"
    np.testing.assert_allclose(beam.dim, _BEAML_DIMS, rtol=1e-6)
    comp = model.properties[13]
    assert comp.type == "PCOMPG"
    assert comp.global_ply_ids == [1, 2, 3]
    np.testing.assert_allclose(comp.thicknesses, [0.1, 0.2, 0.1], rtol=1e-6)
    np.testing.assert_allclose(comp.thetas, [0.0, 45.0, -45.0])


# -- output coordinate systems (CD) -----------------------------------------------------


def _output_system_mesh():
    """A triangle whose GRIDs 1 and 2 are in cylindrical system 5 and GRID 3 in spherical
    system 6, with a displacement, a complex displacement and totals of grid point forces.
    """
    points = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [1.0, 1.0, 1.0]])
    mesh = meshioplusplus.Mesh(points, [("triangle", np.array([[0, 1, 2]]))])
    mesh.point_data["nastran:cd"] = np.array([5, 5, 6], dtype=np.int64)
    axes = {"A": [0.0, 0.0, 0.0], "X": [1.0, 0.0, 0.0], "Y": [0.0, 1.0, 0.0]}
    for cid, kind in ((5, 2), (6, 3)):
        mesh.field_data[f"nastran:CORD2_{cid}:TYPE"] = np.array([kind], dtype=np.int64)
        mesh.field_data[f"nastran:CORD2_{cid}:Z"] = np.array([0.0, 0.0, 1.0])
        for member, vector in axes.items():
            mesh.field_data[f"nastran:CORD2_{cid}:{member}"] = np.array(vector)
    mesh.point_data["DISPLACEMENT"] = np.array(
        [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.3, -0.2, 0.5]]
    )
    mesh.point_data["DISPLACEMENT_ROT"] = np.array(
        [[0.0, 1.0, 0.0], [0.0, 0.0, 2.0], [0.1, 0.2, 0.3]]
    )
    mesh.point_data["VELOCITY_real"] = mesh.point_data["DISPLACEMENT"] * 2.0
    mesh.point_data["VELOCITY_imag"] = mesh.point_data["DISPLACEMENT"] * -1.0
    mesh.cell_data["nastran:eid"] = [np.array([7])]
    for k, member in enumerate(("F1", "F2", "F3", "M1", "M2", "M3")):
        mesh.point_data[f"GRID_FORCE:TOTALS:{member}"] = np.array([1.0, 2.0, 3.0]) * (
            k + 1
        )
    return mesh


def test_results_are_written_in_the_output_systems(engine, tmp_path):
    path = tmp_path / "cd.op2"
    mesh = _output_system_mesh()
    mesh.field_data["nastran:analysis"] = np.array([5], dtype=np.int64)
    mesh.field_data["meshio:time"] = np.array([1.0])
    write_with(engine, path, mesh)
    for read in (meshioplusplus.nastran_op2.read, py_op2.read):
        back = read(path)
        np.testing.assert_array_equal(back.point_data["nastran:cd"], [5, 5, 6])
        for name in (
            "DISPLACEMENT",
            "DISPLACEMENT_ROT",
            "VELOCITY_real",
            "VELOCITY_imag",
        ):
            np.testing.assert_allclose(
                back.point_data[name], mesh.point_data[name], atol=2e-6, err_msg=name
            )
        for member in ("F1", "F2", "F3", "M1", "M2", "M3"):
            name = f"GRID_FORCE:TOTALS:{member}"
            np.testing.assert_allclose(
                back.point_data[name].ravel(), mesh.point_data[name], atol=2e-6
            )


def test_the_engines_write_the_same_bytes_with_output_systems(tmp_path):
    mesh = _output_system_mesh()
    core, python = tmp_path / "core.op2", tmp_path / "python.op2"
    _core.nastran_op2_write(str(core), mesh)
    py_write.write(str(python), mesh)
    assert core.read_bytes() == python.read_bytes()


def test_pynastran_reads_the_output_systems(tmp_path):
    op2 = pytest.importorskip("pyNastran.op2.op2_geom", exc_type=ImportError)
    path = tmp_path / "pn_cd.op2"
    mesh = _output_system_mesh()
    for name in ("VELOCITY_real", "VELOCITY_imag"):
        del mesh.point_data[name]
    _core.nastran_op2_write(str(path), mesh)
    model = op2.read_op2_geom(str(path), debug=False)
    assert [model.nodes[n].cd for n in (1, 2, 3)] == [5, 5, 6]
    assert model.coords[5].type == "CORD2C" and model.coords[6].type == "CORD2S"
    shown = model.displacements[1].data[
        0
    ]  # the output components, per node in id order
    # GRID 2 at (0, 2, 0): basic +x is -theta; GRID 1 at (1, 0, 0): +x is +r
    np.testing.assert_allclose(shown[0, :3], [1.0, 0.0, 0.0], atol=1e-6)
    np.testing.assert_allclose(shown[1, :3], [0.0, -1.0, 0.0], atol=1e-6)


def test_an_output_system_needs_its_definition(engine, tmp_path):
    mesh = _output_system_mesh()
    del mesh.field_data["nastran:CORD2_6:TYPE"]
    with pytest.raises(Exception, match="coordinate system 6"):
        write_with(engine, tmp_path / "x.op2", mesh)
