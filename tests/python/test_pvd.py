"""ParaView collection `.pvd` (v15.0.0).

A time-indexed list of VTK XML files. `timestep` selects the step (`time_step=`)
and, within a step, `part` selects the piece (`piece=`). Behavioural tests run
against both engines: the C++ core (strict-core on, so a silent fallback would
fail loudly) and the pure-Python reference.
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import ReadError, WriteError, _core
from meshioplusplus._fallback import set_strict_core
from meshioplusplus.pvd import SeriesWriter, _pvd

XY = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
TRI = [[0, 1, 2]]


def _step(k, npts=3):
    """A one-triangle mesh shifted by k, carrying u = k."""
    pts = np.array(XY, dtype=np.float64) + k
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array(TRI))])
    mesh.point_data["u"] = np.full(3, float(k))
    return mesh


def _num_cells(mesh):
    return sum(len(cb.data) for cb in mesh.cells)


class _Engine:
    def __init__(self, mod):
        self.read = mod.read
        self.write = mod.write


@pytest.fixture(params=["core", "python"])
def engine(request):
    if request.param == "core":
        if not hasattr(_core, "pvd_read"):
            pytest.skip("this build has no pvd core")
        set_strict_core(True)
        yield _Engine(meshioplusplus.pvd)
        set_strict_core(None)
    else:
        yield _Engine(_pvd)


def _write_vtu(path, points=XY, tris=TRI):
    pts = " ".join(str(v) for row in points for v in row)
    conn = " ".join(str(v) for row in tris for v in row)
    offs = " ".join(str(3 * (i + 1)) for i in range(len(tris)))
    types = " ".join("5" for _ in tris)
    path.write_text(
        '<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" '
        'byte_order="LittleEndian">\n<UnstructuredGrid>\n'
        f'<Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(tris)}">\n<Points>'
        f'<DataArray type="Float64" NumberOfComponents="3" format="ascii">{pts}'
        "</DataArray></Points>\n<Cells>"
        f'<DataArray type="Int64" Name="connectivity" format="ascii">{conn}</DataArray>'
        f'<DataArray type="Int64" Name="offsets" format="ascii">{offs}</DataArray>'
        f'<DataArray type="UInt8" Name="types" format="ascii">{types}</DataArray>'
        "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n"
    )


def _write_pvd(path, entries):
    """`entries`: dicts of DataSet attributes (timestep/part/group/name/file)."""
    rows = []
    for e in entries:
        attrs = " ".join(f'{k}="{v}"' for k, v in e.items())
        rows.append(f"<DataSet {attrs}/>")
    path.write_text(
        '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" '
        'byte_order="LittleEndian">\n<Collection>\n'
        + "\n".join(rows)
        + "\n</Collection>\n</VTKFile>\n"
    )


def _series(tmp_path, times=(0.0, 0.5, 2.0), name="t.pvd"):
    path = tmp_path / name
    meshioplusplus.write_sequence(path, ((t, _step(i)) for i, t in enumerate(times)))
    return path


# --- write -------------------------------------------------------------------


def test_plain_write_is_a_one_step_collection(engine, tmp_path):
    mesh = _step(1)
    mesh.field_data["meshio:time"] = np.array([0.25])
    meshioplusplus.pvd.write(tmp_path / "c.pvd", mesh)
    text = (tmp_path / "c.pvd").read_text()
    assert '<VTKFile type="Collection"' in text
    assert text.count("<DataSet ") == 1
    assert 'timestep="0.25"' in text and 'file="c/c_0000.vtu"' in text
    assert (tmp_path / "c" / "c_0000.vtu").is_file()
    back = engine.read(tmp_path / "c.pvd")
    assert float(np.ravel(back.field_data["meshio:time"])[0]) == 0.25
    np.testing.assert_allclose(back.points, mesh.points)


def test_a_mesh_without_a_time_is_step_zero(tmp_path):
    meshioplusplus.pvd.write(tmp_path / "c.pvd", _step(0))
    assert 'timestep="0' in (tmp_path / "c.pvd").read_text()


def test_the_write_sequence_layout(tmp_path):
    path = _series(tmp_path)
    text = path.read_text()
    assert text.count("<DataSet ") == 3
    for i, t in enumerate(("0", "0.5", "2")):
        assert f'timestep="{t}' in text
        assert (
            f'file="t/t_{i:04d}.vtu"' in text
            and (tmp_path / "t" / f"t_{i:04d}.vtu").is_file()
        )


def test_a_run_that_is_killed_leaves_an_openable_collection(tmp_path):
    """The index is rewritten after every step, not only on close."""
    writer = SeriesWriter(tmp_path / "k.pvd")
    writer.write(0.0, _step(0))
    writer.write(1.0, _step(1))
    # no close(): simulate the process dying here
    assert (tmp_path / "k.pvd").read_text().count("<DataSet ") == 2
    assert meshioplusplus.read_metadata(tmp_path / "k.pvd")["time_values"] == [0.0, 1.0]


def test_a_series_with_no_step_is_refused(tmp_path):
    with pytest.raises(WriteError, match="at least one step"):
        with SeriesWriter(tmp_path / "e.pvd"):
            pass


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_a_non_finite_time_is_refused(bad, tmp_path):
    with pytest.raises(WriteError, match="must be finite"):
        SeriesWriter(tmp_path / "n.pvd").write(bad, _step(0))


def test_a_path_with_a_space_round_trips(tmp_path):
    path = _series(tmp_path, name="my run.pvd")
    assert (tmp_path / "my run" / "my run_0000.vtu").is_file()
    assert meshioplusplus.read(path, time_step=1).points[0][0] == 1.0


# --- read: the time axis -----------------------------------------------------


def test_time_step_selection(engine, tmp_path):
    path = _series(tmp_path)
    assert engine.read(path).points[0][0] == 0.0
    assert engine.read(path, time_step=1).points[0][0] == 1.0
    last = engine.read(path, time_step=-1)
    assert last.points[0][0] == 2.0
    assert float(np.ravel(last.field_data["meshio:time"])[0]) == 2.0
    for k in (3, -4):
        with pytest.raises(ReadError, match=r"time step .* is out of range.*3 steps"):
            engine.read(path, time_step=k)


def test_steps_are_the_distinct_times_in_ascending_order(engine, tmp_path):
    """Entries need not be listed in time order, nor uniformly spaced."""
    for i in range(3):
        _write_vtu(tmp_path / f"p{i}.vtu", [[v + i for v in row] for row in XY])
    _write_pvd(
        tmp_path / "u.pvd",
        [
            {"timestep": 2.5, "file": "p0.vtu"},
            {"timestep": 0.1, "file": "p1.vtu"},
            {"timestep": 1, "file": "p2.vtu"},
        ],
    )
    assert engine.read(tmp_path / "u.pvd", time_step=0).points[0][0] == 1.0  # t=0.1
    assert engine.read(tmp_path / "u.pvd", time_step=1).points[0][0] == 2.0  # t=1
    assert engine.read(tmp_path / "u.pvd", time_step=2).points[0][0] == 0.0  # t=2.5
    assert meshioplusplus.read_metadata(tmp_path / "u.pvd")["time_values"] == [
        0.1,
        1.0,
        2.5,
    ]


def test_an_entry_without_a_timestep_is_step_zero(engine, tmp_path):
    _write_vtu(tmp_path / "a.vtu")
    _write_pvd(tmp_path / "n.pvd", [{"file": "a.vtu"}])
    back = engine.read(tmp_path / "n.pvd")
    assert float(np.ravel(back.field_data["meshio:time"])[0]) == 0.0
    assert _num_cells(back) == 1


def test_read_sequence_closes_on_time(tmp_path):
    path = _series(tmp_path, times=(0.0, 0.5, 2.0))
    got = [(t, float(m.points[0][0])) for t, m in meshioplusplus.read_sequence(path)]
    assert got == [(0.0, 0.0), (0.5, 1.0), (2.0, 2.0)]
    ts = meshioplusplus.TimeSeries(path)
    assert len(ts) == 3 and ts.times == [0.0, 0.5, 2.0]


def test_fan_out_one_file_per_step(tmp_path):
    path = _series(tmp_path)
    written = meshioplusplus.write_sequence(
        str(tmp_path / "out_{step}.vtu"), meshioplusplus.read_sequence(path)
    )
    assert len(written) == 3
    assert meshioplusplus.read(written[2]).points[0][0] == 2.0


def test_metadata_reports_every_time_without_opening_a_later_step(tmp_path):
    """time_values come off the index alone, so step 1.. need not be readable."""
    path = _series(tmp_path)
    for i in (1, 2):
        (tmp_path / "t" / f"t_{i:04d}.vtu").unlink()
    meta = meshioplusplus.read_metadata(path)
    assert meta["time_values"] == [0.0, 0.5, 2.0]
    assert meta["format"] == "pvd"
    if hasattr(_core, "pvd_read"):
        assert meta["fell_back_to_full_read"] is False


def test_num_steps_uses_the_index(tmp_path):
    from meshioplusplus import _sequence

    assert _sequence.num_steps(_series(tmp_path)) == 3


# --- read: the part axis -----------------------------------------------------


def _two_steps_three_parts(tmp_path):
    for t in range(2):
        for p in range(3):
            _write_vtu(
                tmp_path / f"s{t}p{p}.vtu", [[v + 10 * t + p for v in r] for r in XY]
            )
    entries = [
        {"timestep": t, "part": p, "file": f"s{t}p{p}.vtu"}
        for t in range(2)
        for p in range(3)
    ]
    _write_pvd(tmp_path / "tp.pvd", entries)
    return tmp_path / "tp.pvd"


def test_part_and_timestep_are_orthogonal(engine, tmp_path):
    path = _two_steps_three_parts(tmp_path)
    step0 = engine.read(path)
    assert _num_cells(step0) == 3 and len(step0.points) == 9
    assert [r.name for r in step0.regions] == ["part_0", "part_1", "part_2"]
    step1 = engine.read(path, time_step=1)
    assert step1.points[0][0] == 10.0 and len(step1.points) == 9
    one = engine.read(path, time_step=1, piece=2)
    assert one.points[0][0] == 12.0 and _num_cells(one) == 1
    assert one.regions == [] or len(one.regions) == 0
    assert engine.read(path, time_step=1, piece=-1).points[0][0] == 12.0
    with pytest.raises(ReadError, match=r"piece 3 is out of range.*3 pieces"):
        engine.read(path, piece=3)


def test_parts_are_ordered_by_part_then_document_order(engine, tmp_path):
    for i in range(3):
        _write_vtu(tmp_path / f"p{i}.vtu", [[v + i for v in r] for r in XY])
    _write_pvd(
        tmp_path / "o.pvd",
        [
            {"timestep": 0, "part": 2, "file": "p2.vtu"},
            {"timestep": 0, "part": 0, "file": "p0.vtu"},
            {"timestep": 0, "part": 1, "file": "p1.vtu"},
        ],
    )
    assert [r.name for r in engine.read(tmp_path / "o.pvd").regions] == [
        "part_0",
        "part_1",
        "part_2",
    ]
    assert engine.read(tmp_path / "o.pvd", piece=0).points[0][0] == 0.0


def test_region_names_come_from_name_then_group(engine, tmp_path):
    for i in range(3):
        _write_vtu(tmp_path / f"p{i}.vtu", [[v + i for v in r] for r in XY])
    _write_pvd(
        tmp_path / "g.pvd",
        [
            {"timestep": 0, "part": 0, "name": "wing", "file": "p0.vtu"},
            {"timestep": 0, "part": 1, "group": "fluid", "file": "p1.vtu"},
            {"timestep": 0, "part": 2, "file": "p2.vtu"},
        ],
    )
    names = sorted(r.name for r in engine.read(tmp_path / "g.pvd").regions)
    assert names == ["fluid/part_1", "part_2", "wing"]


# --- composition: .pvd -> .pvtu -> .vtu --------------------------------------


def _partitioned_run(tmp_path, nsteps=3, nparts=3, ghost_layers=0):
    """`.pvd` -> `.pvtu` -> `.vtu` over a real partitioned transient run."""
    if ghost_layers and not hasattr(_core, "partition"):
        pytest.skip("ghost layers need the C++ core")
    n = 6
    xs, ys = np.meshgrid(np.arange(n + 1.0), np.arange(n + 1.0))
    base = np.c_[xs.ravel(), ys.ravel(), np.zeros((n + 1) ** 2)]
    tris = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            b, c = a + 1, a + n + 1
            tris += [[a, b, c + 1], [a, c + 1, c]]
    entries = []
    for t in range(nsteps):
        mesh = meshioplusplus.Mesh(base + [0, 0, t], [("triangle", np.array(tris))])
        mesh.point_data["u"] = np.full(len(base), float(t))
        pieces = meshioplusplus.partition(mesh, nparts, ghost_layers=ghost_layers)
        meshioplusplus.pvtu.write_pieces(tmp_path / f"step{t}.pvtu", pieces)
        entries.append({"timestep": 0.5 * t, "part": 0, "file": f"step{t}.pvtu"})
    _write_pvd(tmp_path / "run.pvd", entries)
    return tmp_path / "run.pvd", base, tris


def test_pvd_of_pvtu_of_vtu_nests(engine, tmp_path):
    path, base, tris = _partitioned_run(tmp_path)
    for t in range(3):
        back = engine.read(path, time_step=t)
        # a step that is one .pvtu reads as that file does: its own piece regions
        assert [r.name for r in back.regions] == ["piece_0", "piece_1", "piece_2"]
        assert _num_cells(back) == len(tris)
        assert float(np.ravel(back.field_data["meshio:time"])[0]) == 0.5 * t
        welded = meshioplusplus.clean(back, weld=True)
        assert len(welded.points) == len(base)
        np.testing.assert_allclose(welded.point_data["u"], t)
    assert meshioplusplus.read_metadata(path)["time_values"] == [0.0, 0.5, 1.0]


def test_the_piece_axis_selects_a_part_not_a_pvtu_piece(engine, tmp_path):
    """One .pvtu per step is one part; `piece=` is about the .pvd's own axis."""
    path, _base, tris = _partitioned_run(tmp_path)
    whole = engine.read(path, time_step=1, piece=0)
    assert _num_cells(whole) == len(tris)


def test_ghost_policy_reaches_the_children(engine, tmp_path):
    path, _base, tris = _partitioned_run(tmp_path, nsteps=2, ghost_layers=1)
    kept = engine.read(path)
    assert _num_cells(kept) > len(tris) and "vtkGhostType" in kept.cell_data
    dropped = engine.read(path, ghosts="drop")
    assert _num_cells(dropped) == len(tris)
    assert "vtkGhostType" not in dropped.cell_data


def test_ghost_policy_applies_to_a_plain_vtu_child(engine, tmp_path):
    piece = _step(0)
    piece.cell_data["vtkGhostType"] = [np.array([1], dtype=np.uint8)]
    meshioplusplus.vtu.write(tmp_path / "g.vtu", piece)
    _write_pvd(tmp_path / "g.pvd", [{"timestep": 0, "file": "g.vtu"}])
    assert _num_cells(engine.read(tmp_path / "g.pvd")) == 1
    assert _num_cells(engine.read(tmp_path / "g.pvd", ghosts="drop")) == 0


def test_a_vtm_child_is_read(engine, tmp_path):
    mesh = _step(0)
    mesh.cells.append(meshioplusplus.CellBlock("triangle", np.array([[0, 1, 2]])))
    meshioplusplus.vtm.write(tmp_path / "m.vtm", mesh)
    _write_pvd(tmp_path / "v.pvd", [{"timestep": 0, "file": "m.vtm"}])
    assert _num_cells(engine.read(tmp_path / "v.pvd")) == 2


# --- refusals ----------------------------------------------------------------


def test_a_missing_piece_names_the_attribute_and_the_file(engine, tmp_path):
    _write_pvd(tmp_path / "m.pvd", [{"timestep": 0, "file": "gone/x.vtu"}])
    with pytest.raises(ReadError, match=r"gone/x\.vtu.*\(file=\).*does not exist"):
        engine.read(tmp_path / "m.pvd")


def test_an_entry_without_a_file_is_refused(engine, tmp_path):
    _write_pvd(tmp_path / "f.pvd", [{"timestep": 0}])
    with pytest.raises(ReadError, match="missing its 'file' attribute"):
        engine.read(tmp_path / "f.pvd")


def test_legacy_vtk_and_foreign_entries_are_refused_by_name(engine, tmp_path):
    (tmp_path / "old.vtk").write_text("# vtk DataFile Version 3.0\n")
    _write_pvd(tmp_path / "l.pvd", [{"timestep": 0, "file": "old.vtk"}])
    with pytest.raises(
        ReadError, match=r"unsupported piece.*old\.vtk.*\.vtu/\.vtp/\.vtm/\.pvtu/\.pvtp"
    ):
        engine.read(tmp_path / "l.pvd")


def test_a_bad_timestep_is_refused(engine, tmp_path):
    _write_vtu(tmp_path / "a.vtu")
    _write_pvd(tmp_path / "b.pvd", [{"timestep": "soon", "file": "a.vtu"}])
    with pytest.raises(ReadError, match="bad timestep/part"):
        engine.read(tmp_path / "b.pvd")


def test_wrong_root_type_and_not_xml(engine, tmp_path):
    (tmp_path / "w.pvd").write_text(
        '<VTKFile type="vtkMultiBlockDataSet"><Collection/></VTKFile>'
    )
    with pytest.raises(ReadError, match="expected type Collection"):
        engine.read(tmp_path / "w.pvd")
    (tmp_path / "x.pvd").write_text("nope")
    with pytest.raises(ReadError):
        engine.read(tmp_path / "x.pvd")


def test_an_empty_collection_reads_as_an_empty_mesh(engine, tmp_path):
    _write_pvd(tmp_path / "e.pvd", [])
    back = engine.read(tmp_path / "e.pvd")
    assert len(back.points) == 0 and back.cells == []


def test_a_bad_ghost_policy_is_refused(engine, tmp_path):
    _write_vtu(tmp_path / "a.vtu")
    _write_pvd(tmp_path / "g.pvd", [{"timestep": 0, "file": "a.vtu"}])
    with pytest.raises(ValueError, match="ghosts must be 'keep' or 'drop'"):
        engine.read(tmp_path / "g.pvd", ghosts="x")


def test_buffers_are_refused(tmp_path):
    import io

    with pytest.raises(Exception, match="pvd"):
        meshioplusplus.write(io.StringIO(), _step(0), file_format="pvd")


# --- the engines agree -------------------------------------------------------


@pytest.mark.skipif(not hasattr(_core, "pvd_read"), reason="no pvd core")
@pytest.mark.parametrize(
    "binary, compression", [(False, None), (True, None), (True, "zlib")]
)
def test_cross_compat(binary, compression, tmp_path):
    set_strict_core(True)
    try:
        (tmp_path / "core").mkdir()
        (tmp_path / "py").mkdir()
        mesh = _step(3)
        mesh.field_data["meshio:time"] = np.array([1.5])
        meshioplusplus.pvd.write(
            tmp_path / "core" / "c.pvd", mesh, binary=binary, compression=compression
        )
        _pvd.write(
            tmp_path / "py" / "c.pvd", mesh, binary=binary, compression=compression
        )
        for p in (tmp_path / "core" / "c.pvd", tmp_path / "py" / "c.pvd"):
            a, b = meshioplusplus.pvd.read(p), _pvd.read(p)
            np.testing.assert_allclose(a.points, b.points)
            assert float(np.ravel(a.field_data["meshio:time"])[0]) == 1.5
            assert float(np.ravel(b.field_data["meshio:time"])[0]) == 1.5

        def body(p):
            return [
                ln for ln in p.read_text().splitlines() if not ln.startswith("<!--")
            ]

        assert body(tmp_path / "core" / "c.pvd") == body(tmp_path / "py" / "c.pvd")
    finally:
        set_strict_core(None)


# --- the file's own time: TimeValue (v15.0.0) ----------------------------------


def _write_timed(path, k, time_value=True, meshio_time=False):
    mesh = _step(k)
    if time_value:
        mesh.field_data["TimeValue"] = np.array([10.0 * k + 0.5])
    if meshio_time:
        mesh.field_data[TIME_KEY] = np.array([100.0 * k])
    meshioplusplus.vtu.write(path, mesh)


TIME_KEY = "meshio:time"


def test_an_entry_without_a_timestep_takes_its_time_from_the_files_time_value(
    engine, tmp_path
):
    """The per-file alternative to `timestep=` (VTK's "time in field data")."""
    _write_timed(tmp_path / "a.vtu", 3)  # TimeValue 30.5
    _write_timed(tmp_path / "b.vtu", 1, meshio_time=True)  # TimeValue wins: 10.5
    _write_timed(tmp_path / "c.vtu", 2, time_value=False, meshio_time=True)  # 200
    _write_timed(tmp_path / "e.vtu", 4, time_value=False)  # neither: 0
    _write_pvd(
        tmp_path / "t.pvd",
        [{"file": "a.vtu"}, {"file": "b.vtu"}, {"file": "c.vtu"}, {"file": "e.vtu"}],
    )
    times = [0.0, 10.5, 30.5, 200.0]
    assert meshioplusplus.read_metadata(tmp_path / "t.pvd")["time_values"] == times
    expected_x = [4.0, 1.0, 3.0, 2.0]  # ascending time: e, b, a, c
    for step in range(4):
        back = engine.read(tmp_path / "t.pvd", time_step=step)
        assert back.points[0][0] == expected_x[step]
        assert float(np.ravel(back.field_data[TIME_KEY])[0]) == times[step]


def test_an_explicit_timestep_outranks_the_files_own_time(engine, tmp_path):
    _write_timed(tmp_path / "a.vtu", 1, meshio_time=True)
    _write_pvd(tmp_path / "t.pvd", [{"timestep": 7, "file": "a.vtu"}])
    assert meshioplusplus.read_metadata(tmp_path / "t.pvd")["time_values"] == [7.0]


def test_a_missing_child_of_a_timeless_entry_names_the_file(engine, tmp_path):
    _write_pvd(tmp_path / "t.pvd", [{"file": "gone.vtu"}])
    with pytest.raises(ReadError, match=r"gone\.vtu"):
        engine.read(tmp_path / "t.pvd")


def test_a_pvd_of_a_pvtu_whose_pieces_carry_time_value(engine, tmp_path):
    """Every piece repeats it (union, not `0:TimeValue`), and the step is timed by it."""
    for step in range(2):
        pieces = []
        for _ in range(2):
            piece = _step(step)
            piece.field_data["TimeValue"] = np.array([0.5 * (step + 1)])
            pieces.append(piece)
        meshioplusplus.pvtu.write_pieces(tmp_path / f"s{step}.pvtu", pieces)
    _write_pvd(tmp_path / "run.pvd", [{"file": "s0.pvtu"}, {"file": "s1.pvtu"}])
    assert meshioplusplus.read_metadata(tmp_path / "run.pvd")["time_values"] == [
        0.5,
        1.0,
    ]
    assert (
        float(
            np.ravel(
                engine.read(tmp_path / "run.pvd", time_step=1).field_data[TIME_KEY]
            )[0]
        )
        == 1.0
    )


# --- fan-out writes the time into each file (v15.0.0) ---------------------------


def test_a_fan_out_to_vtu_records_each_steps_time_in_the_file(tmp_path):
    path = _series(tmp_path, times=(0.0, 0.5, 2.0))
    written = meshioplusplus.write_sequence(
        str(tmp_path / "out_{step}.vtu"), meshioplusplus.read_sequence(path)
    )
    times = [
        float(np.ravel(meshioplusplus.read(p).field_data[TIME_KEY])[0]) for p in written
    ]
    assert times == [0.0, 0.5, 2.0]
