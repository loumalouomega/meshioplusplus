"""Native/reference sets and stateful Exodus output, including lifecycle errors."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core
from meshioplusplus.exodus import _exodus as reference

from . import helpers

netCDF4 = pytest.importorskip("netCDF4")


def _mesh():
    mesh = helpers.tet_mesh.copy()
    mesh.point_data["temperature"] = np.arange(len(mesh.points), dtype=float)
    mesh.cell_data["stress"] = [np.ones(len(c)) for c in mesh.cells]
    mesh.cell_data["exodus:attr:RADIUS"] = [np.full(len(c), 0.5) for c in mesh.cells]
    mesh.regions = [
        mp.Region("nodes with a long name " + "x" * 40, "point", [0, 2], tag=41),
        mp.Region("empty nodes", "point", [], tag=42),
        mp.Region("wall", "side", [[0, 1], [0, 3]], tag=91),
        mp.Region("empty sides", "side", np.empty((0, 2), dtype=np.int64), tag=92),
    ]
    return mesh


def _region_key(mesh):
    return sorted((r.kind, r.name, r.tag, r.entries.tolist()) for r in mesh.regions if r.kind != "cell")


@pytest.fixture(params=["reference", "native"])
def writer(request):
    if request.param == "reference":
        return reference.write
    if not getattr(_core, "__has_netcdf__", False):
        pytest.skip("requires netCDF")
    return lambda path, mesh: _core.exodus_write(str(path), mesh)


def test_sets_preserve_ids_names_empty_groups_and_side_numbers(tmp_path, writer):
    mesh = _mesh()
    path = tmp_path / "sets.e"
    writer(path, mesh)
    for reader in [reference.read, mp.exodus.read]:
        assert _region_key(reader(path)) == _region_key(mesh)
    with netCDF4.Dataset(path) as file:
        assert sorted(file.variables["ns_prop1"][:]) == [41, 42]
        assert sorted(file.variables["ss_prop1"][:]) == [91, 92]
        assert sorted(file.variables["side_ss2"][:]) == [2, 4]


@pytest.mark.parametrize("kind,entries,error", [
    ("point", [99], "invalid point"),
    ("side", [[99, 0]], "invalid cell"),
    ("side", [[0, 99]], "unsupported facet"),
])
def test_invalid_sets_do_not_truncate_target(tmp_path, writer, kind, entries, error):
    mesh = _mesh()
    mesh.regions = [mp.Region("bad", kind, entries)]
    path = tmp_path / "existing.e"
    path.write_bytes(b"keep")
    with pytest.raises(Exception, match=error):
        writer(path, mesh)
    assert path.read_bytes() == b"keep"


@pytest.fixture(params=["reference", "native"])
def series_cls(request):
    if request.param == "reference":
        return reference.TimeSeriesWriter
    if not hasattr(_core, "ExodusTimeSeriesWriter"):
        pytest.skip("requires rebuilt Exodus series writer")
    return _core.ExodusTimeSeriesWriter


def test_series_steps_static_attributes_and_lifecycle(tmp_path, series_cls):
    path = tmp_path / "series.e"
    mesh = _mesh()
    with series_cls(str(path)) as series:
        with pytest.raises(Exception, match="write_points_cells"):
            series.write_data(0, mesh)
        series.write_points_cells(mesh)
        with pytest.raises(Exception, match="new, open"):
            series.write_points_cells(mesh)
        for k, time in enumerate([0.123456789012345, 2.5, 9.0]):
            step = mesh.copy()
            step.point_data["temperature"] += k
            step.cell_data["stress"][0] += k
            step.cell_data["exodus:attr:RADIUS"][0][:] = 88
            series.write_data(time, step)
        assert series.num_steps == 3
        bad = mesh.copy()
        bad.point_data["new"] = np.ones(len(mesh.points))
        with pytest.raises(Exception, match="schema changed"):
            series.write_data(10, bad)
        with pytest.raises(Exception, match="finite"):
            series.write_data(float("nan"), mesh)
        assert series.num_steps == 3
        series.flush()
    assert series.finalized
    series.finalize()
    series.flush()
    with pytest.raises(Exception, match="open series"):
        series.write_data(10, mesh)
    with netCDF4.Dataset(path) as file:
        np.testing.assert_array_equal(file.variables["time_whole"][:], [0.123456789012345, 2.5, 9.0])
    for reader in [reference.read, mp.exodus.read]:
        back = reader(path, time_step=-1)
        np.testing.assert_array_equal(back.point_data["temperature"], mesh.point_data["temperature"] + 2)
        np.testing.assert_array_equal(back.cell_data["exodus:attr:RADIUS"][0], mesh.cell_data["exodus:attr:RADIUS"][0])
        assert _region_key(back) == _region_key(mesh)
