"""Transient VTKHDF: the h5py writer/reader, the C++ writer, and the sequence engine.

Geometry x 1, fields x N: one static grid, then one step at a time. The C++ class is
reachable as ``_core.VtkhdfTimeSeriesWriter``; tests that need it skip on a build that
predates it.
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core, _sequence
from meshioplusplus.vtkhdf import TimeSeriesReader, TimeSeriesWriter, _vtkhdf

h5py = pytest.importorskip("h5py")

has_cpp_writer = hasattr(_core, "VtkhdfTimeSeriesWriter") and getattr(
    _core, "__has_hdf5__", False
)
needs_cpp_writer = pytest.mark.skipif(
    not has_cpp_writer, reason="_core predates VtkhdfTimeSeriesWriter / no HDF5"
)

POINTS = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float)
CELLS = [("tetra", np.array([[0, 1, 2, 3], [1, 2, 3, 4]]))]


def _step(k, field=False):
    point = {"u": np.full(5, float(k)), "vel": np.full((5, 3), float(k))}
    cell = {"p": [np.array([10.0 * k, 10.0 * k + 1])]}
    fd = {"gain": np.array([k, 2.0 * k, 3.0 * k])} if field else None
    return point, cell, fd


def _mesh_step(k, field=False):
    point, cell, fd = _step(k, field)
    return meshioplusplus.Mesh(
        POINTS.copy(),
        [(c.type, c.data) for c in _grid().cells],
        point_data=point,
        cell_data=cell,
        field_data=fd or {},
    )


def _grid():
    return meshioplusplus.Mesh(POINTS, CELLS)


def _write_py(path, n=4, field=False, **kw):
    with TimeSeriesWriter(path, **kw) as w:
        w.write_points_cells(POINTS, CELLS)
        for k in range(n):
            point, cell, fd = _step(k, field)
            w.write_data(0.5 * k, point_data=point, cell_data=cell, field_data=fd)
    return path


def _write_cpp(path, n=4, field=False, **kw):
    with _core.VtkhdfTimeSeriesWriter(str(path), **kw) as w:
        w.write_points_cells(_grid())
        for k in range(n):
            w.write_data(0.5 * k, _mesh_step(k, field))
    return path


WRITERS = [
    pytest.param(_write_py, id="python"),
    pytest.param(_write_cpp, id="cpp", marks=needs_cpp_writer),
]


@pytest.mark.parametrize("write", WRITERS)
def test_every_step_reads_back(write, tmp_path):
    p = write(tmp_path / "s.vtkhdf", 4, field=True)
    for k in range(4):
        m = meshioplusplus.read(p, time_step=k)
        assert m.point_data["u"].tolist() == [float(k)] * 5
        assert m.point_data["vel"].shape == (5, 3)
        assert m.cell_data["p"][0].tolist() == [10.0 * k, 10.0 * k + 1]
        assert m.field_data["gain"].tolist() == [k, 2.0 * k, 3.0 * k]
        assert m.field_data[_vtkhdf.TIME_KEY].tolist() == [0.5 * k]
        assert m.cells[0].data.tolist() == [[0, 1, 2, 3], [1, 2, 3, 4]]
    assert meshioplusplus.read(p, time_step=-1).point_data["u"][0] == 3.0
    with pytest.raises(meshioplusplus.ReadError, match="4 step"):
        meshioplusplus.read(p, time_step=4)


@pytest.mark.parametrize("write", WRITERS)
def test_on_disk_layout(write, tmp_path):
    p = write(tmp_path / "s.vtkhdf", 3, field=True)
    with h5py.File(p, "r") as f:
        g = f["VTKHDF"]
        assert g.attrs["Type"] == b"UnstructuredGrid"
        assert g.attrs["Version"].tolist() == [2, 0]
        assert g["NumberOfPoints"].shape == (1,) and g["Points"].shape == (
            5,
            3,
        )  # geometry x 1
        assert g["PointData/u"].shape == (15,) and g["PointData/vel"].shape == (
            15,
            3,
        )  # fields x N
        assert g["CellData/p"].shape == (6,)
        for name in ("PointData/u", "CellData/p", "Steps/Values", "FieldData/gain"):
            assert g[name].chunks is not None, name  # appended, so chunked
        s = g["Steps"]
        assert s.attrs["NSteps"] == 3
        assert s["Values"][:].tolist() == [0.0, 0.5, 1.0]
        for name in (
            "PartOffsets",
            "PointOffsets",
            "CellOffsets",
            "ConnectivityIdOffsets",
        ):
            assert s[name][:].tolist() == [0, 0, 0], name
        assert s["NumberOfParts"][:].tolist() == [1, 1, 1]
        assert s["PointDataOffsets/u"][:].tolist() == [0, 5, 10]
        assert s["CellDataOffsets/p"][:].tolist() == [0, 2, 4]
        assert s["FieldDataSizes/gain"].shape == (3, 2) and s["FieldDataSizes/gain"][
            0
        ].tolist() == [1, 3]
        assert s["FieldDataOffsets/gain"][:].tolist() == [0, 3, 6]


@needs_cpp_writer
def test_both_writers_produce_the_same_layout(tmp_path):
    a = _write_py(tmp_path / "py.vtkhdf", 3, field=True)
    b = _write_cpp(tmp_path / "cpp.vtkhdf", 3, field=True)

    def shape_of(path):
        out = {}
        with h5py.File(path, "r") as f:

            def visit(name, obj):
                if isinstance(obj, h5py.Dataset):
                    out[name] = (
                        obj.shape,
                        str(obj.dtype),
                        obj[()].tolist() if obj.size < 64 else "big",
                    )

            f["VTKHDF"].visititems(visit)
        return out

    sa, sb = shape_of(a), shape_of(b)
    assert sa.keys() == sb.keys()
    for name in sa:
        assert sa[name] == sb[name], name


def test_steps_are_durable_before_finalize(tmp_path):
    """Nothing is buffered: a run that is killed mid-series leaves a readable file."""
    p = tmp_path / "live.vtkhdf"
    with TimeSeriesWriter(p) as w:
        w.write_points_cells(POINTS, CELLS)
        for k in (1, 2):
            point, cell, _ = _step(k)
            w.write_data(0.5 * (k - 1), point_data=point, cell_data=cell)
            # still open, not finalized: yet every completed step is already readable
            live = meshioplusplus.read(p, time_step=-1)
            assert live.point_data["u"][0] == float(k)
            assert w.num_steps == k
    assert meshioplusplus.read(p, time_step=-1).point_data["u"][0] == 2.0


def test_gzip_option(tmp_path):
    p = _write_py(tmp_path / "gz.vtkhdf", 3, compression="gzip", compression_opts=4)
    with h5py.File(p, "r") as f:
        assert f["VTKHDF/PointData/u"].compression == "gzip"
    assert meshioplusplus.read(p, time_step=2).point_data["u"][0] == 2.0


def test_the_name_set_is_fixed_at_the_first_step(tmp_path):
    with TimeSeriesWriter(tmp_path / "s.vtkhdf") as w:
        w.write_points_cells(POINTS, CELLS)
        point, cell, _ = _step(0)
        w.write_data(0.0, point_data=point, cell_data=cell)
        with pytest.raises(meshioplusplus.WriteError, match="brand_new"):
            w.write_data(
                1.0, point_data={**point, "brand_new": np.zeros(5)}, cell_data=cell
            )
        with pytest.raises(meshioplusplus.WriteError, match="vel"):
            w.write_data(1.0, point_data={"u": point["u"]}, cell_data=cell)
        with pytest.raises(meshioplusplus.WriteError, match="component"):
            w.write_data(
                1.0,
                point_data={"u": point["u"], "vel": np.zeros((5, 2))},
                cell_data=cell,
            )
        assert w.num_steps == 1  # none of the refused steps landed


def test_call_order_and_shapes_are_enforced(tmp_path):
    w = TimeSeriesWriter(tmp_path / "s.vtkhdf")
    with pytest.raises(meshioplusplus.WriteError, match="not open"):
        w.write_points_cells(POINTS, CELLS)
    with w:
        with pytest.raises(
            meshioplusplus.WriteError, match="write_points_cells must be called"
        ):
            w.write_data(0.0, point_data={"u": np.zeros(5)})
        w.write_points_cells(POINTS, CELLS)
        with pytest.raises(meshioplusplus.WriteError, match="already called"):
            w.write_points_cells(POINTS, CELLS)
        with pytest.raises(meshioplusplus.WriteError, match="3 rows"):
            w.write_data(0.0, point_data={"u": np.zeros(3)})
    with pytest.raises(meshioplusplus.WriteError, match="not open"):
        w.write_data(0.0, point_data={"u": np.zeros(5)})


def test_polyhedral_grids_carry_the_extra_step_tables(tmp_path):
    cube = np.array(
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
    faces = [
        np.array(f)
        for f in (
            [0, 3, 2, 1],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
        )
    ]
    p = tmp_path / "poly.vtkhdf"
    with TimeSeriesWriter(p) as w:
        w.write_points_cells(cube, [("polyhedron8", [faces])])
        for k in range(2):
            w.write_data(k, point_data={"u": np.full(8, float(k))})
    with h5py.File(p, "r") as f:
        assert f["VTKHDF"].attrs["Version"].tolist() == [2, 5]
        for name in (
            "FaceConnectivityOffsets",
            "FaceOffsetsOffsets",
            "PolyhedronToFaceIdOffsets",
        ):
            assert f["VTKHDF/Steps"][name][:].tolist() == [0, 0], name
    m = meshioplusplus.read(p, time_step=1)
    assert m.cells[0].type == "polyhedron8" and m.point_data["u"].tolist() == [1.0] * 8


def test_append_continues_an_existing_series(tmp_path):
    p = tmp_path / "s.vtkhdf"
    _write_py(p, 2, field=True)
    with TimeSeriesWriter(p, mode="append") as w:
        assert w.num_steps == 2
        w.write_points_cells(POINTS, CELLS)  # only checks the grid
        point, cell, fd = _step(2, True)
        w.write_data(1.0, point_data=point, cell_data=cell, field_data=fd)
        with pytest.raises(meshioplusplus.WriteError, match="points"):
            TimeSeriesWriter(p, mode="append")  # constructing is fine...
            w._adopted = True
            w.write_points_cells(POINTS[:3], [("triangle", np.array([[0, 1, 2]]))])
    for k in range(3):
        m = meshioplusplus.read(p, time_step=k)
        assert m.point_data["u"].tolist() == [float(k)] * 5
        assert m.field_data["gain"].tolist() == [k, 2.0 * k, 3.0 * k]
    fresh = tmp_path / "fresh.vtkhdf"
    with TimeSeriesWriter(
        fresh, mode="append"
    ) as w:  # a missing path is just a fresh series
        assert w.num_steps == 0
        w.write_points_cells(POINTS, CELLS)
        w.write_data(0.0, point_data={"u": np.zeros(5)})
    assert meshioplusplus.read(fresh).point_data["u"].tolist() == [0.0] * 5


def test_append_refuses_a_file_it_cannot_continue(tmp_path):
    p = tmp_path / "static.vtkhdf"
    meshioplusplus.vtkhdf.write(p, meshioplusplus.Mesh(POINTS, CELLS))
    with pytest.raises(
        meshioplusplus.WriteError, match="not a transient UnstructuredGrid"
    ):
        with TimeSeriesWriter(p, mode="append"):
            pass


def test_time_series_reader(tmp_path):
    p = _write_py(tmp_path / "s.vtkhdf", 3)
    with TimeSeriesReader(p) as r:
        assert r.num_steps == 3 and r.times == [0.0, 0.5, 1.0]
        pts, cells = r.read_points_cells()
        assert pts.shape == (5, 3) and cells[0].type == "tetra"
        t, point, cell = r.read_data(2)
        assert (
            t == 1.0
            and point["u"].tolist() == [2.0] * 5
            and cell["p"][0].tolist() == [20.0, 21.0]
        )
        assert r.read_data(-1)[0] == 1.0
    with pytest.raises(meshioplusplus.ReadError):
        TimeSeriesReader(tmp_path / "missing.vtkhdf")


def test_time_values_reach_read_metadata(tmp_path):
    p = _write_py(tmp_path / "s.vtkhdf", 3)
    meta = meshioplusplus.read_metadata(p)
    assert meta["time_values"] == [0.0, 0.5, 1.0]
    assert meta["num_points"] == 5 and meta["cell_blocks"][0]["type"] == "tetra"


# ---- the sequence engine: the on-disk form of fan-in / fan-out ------------------ #
def test_fan_in_writes_one_transient_file(tmp_path):
    steps = ((0.5 * k, _mesh_step(k, field=True)) for k in range(4))
    out = tmp_path / "series.vtkhdf"
    written = meshioplusplus.write_sequence(out, steps)
    assert written == [str(out)]
    assert _sequence.num_steps(out) == 4
    got = list(meshioplusplus.read_sequence(out))
    assert [t for t, _ in got] == [0.0, 0.5, 1.0, 1.5]
    for k, (_, m) in enumerate(got):
        assert m.point_data["u"].tolist() == [float(k)] * 5
        assert m.field_data["gain"].tolist() == [k, 2.0 * k, 3.0 * k]


def test_fan_out_writes_one_file_per_step(tmp_path):
    src = _write_py(tmp_path / "s.vtkhdf", 3)
    out = meshioplusplus.write_sequence(
        str(tmp_path / "step_{step}.vtu"),
        meshioplusplus.read_sequence(src),
        file_format=None,
    )
    assert len(out) == 3
    for k, path in enumerate(out):
        assert meshioplusplus.read(path).point_data["u"].tolist() == [float(k)] * 5


def test_fan_in_from_per_step_files(tmp_path):
    for k in range(3):
        meshioplusplus.write(tmp_path / f"in_{k}.vtu", _mesh_step(k))
    out = tmp_path / "combined.vtkhdf"
    meshioplusplus.write_sequence(
        out, meshioplusplus.read_sequence(str(tmp_path / "in_*.vtu"))
    )
    assert _sequence.num_steps(out) == 3
    assert meshioplusplus.read(out, time_step=2).point_data["u"].tolist() == [2.0] * 5


def test_a_multi_step_file_is_refused_by_a_single_file_target(tmp_path):
    src = _write_py(tmp_path / "s.vtkhdf", 3)
    with pytest.raises(Exception, match="multi-step"):
        meshioplusplus.write_sequence(
            tmp_path / "x.stl", meshioplusplus.read_sequence(src)
        )
