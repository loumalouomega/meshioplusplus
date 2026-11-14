"""UGRID NetCDF (roadmap 2.4.3.2): files built here with netCDF4, read in both
engines, and cross-checked with xugrid."""

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _fallback
from meshioplusplus.ugrid_nc import _ugrid_nc

netCDF4 = pytest.importorskip("netCDF4")


@pytest.fixture(params=["core", "python"])
def reader(request, monkeypatch):
    if request.param == "core":
        if not getattr(mio._core, "__has_netcdf__", False):
            pytest.skip("this build has no netCDF")
        monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
        monkeypatch.setattr(_fallback, "_strict", True)
        return mio.ugrid_nc.read
    return _ugrid_nc.read


def _mixed_file(path, transposed=False, start=1):
    """Two triangles and a quad (fill-padded), two edges, two time steps."""
    x = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0])
    y = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    faces = np.array([[0, 1, 4, 3], [1, 2, 5, -9], [1, 5, 4, -9]])
    with netCDF4.Dataset(path, "w") as nc:
        nc.createDimension("nNodes", 6)
        nc.createDimension("nFaces", 3)
        nc.createDimension("nMax", 4)
        nc.createDimension("nEdges", 2)
        nc.createDimension("two", 2)
        nc.createDimension("time", 2)
        topo = nc.createVariable("Mesh2", "i4")
        topo.cf_role = "mesh_topology"
        topo.topology_dimension = 2
        topo.node_coordinates = "Mesh2_node_x Mesh2_node_y"
        topo.face_node_connectivity = "Mesh2_face_nodes"
        topo.face_dimension = "nFaces"
        topo.edge_node_connectivity = "Mesh2_edge_nodes"
        nc.createVariable("Mesh2_node_x", "f8", ("nNodes",))[:] = x
        nc.createVariable("Mesh2_node_y", "f8", ("nNodes",))[:] = y
        dims = ("nMax", "nFaces") if transposed else ("nFaces", "nMax")
        fv = nc.createVariable("Mesh2_face_nodes", "i4", dims, fill_value=-9)
        fv.cf_role = "face_node_connectivity"
        fv.start_index = start
        shifted = np.where(faces >= 0, faces + start, -9)
        fv[:] = shifted.T if transposed else shifted
        ev = nc.createVariable("Mesh2_edge_nodes", "i4", ("nEdges", "two"))
        ev.cf_role = "edge_node_connectivity"
        ev.start_index = start
        ev[:] = np.array([[0, 3], [2, 5]]) + start
        nc.createVariable("time", "f8", ("time",))[:] = [10.0, 20.0]
        h = nc.createVariable("waterlevel", "f8", ("time", "nNodes"))
        h.mesh = "Mesh2"
        h.location = "node"
        h[:] = np.vstack([np.arange(6.0), np.arange(6.0) + 100])
        u = nc.createVariable("velocity", "f8", ("time", "nFaces"))
        u.mesh = "Mesh2"
        u.location = "face"
        u[:] = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        q = nc.createVariable("discharge", "f8", ("nEdges",))
        q.mesh = "Mesh2"
        q.location = "edge"
        q[:] = [7.0, 8.0]


@pytest.mark.parametrize("transposed", [False, True])
@pytest.mark.parametrize("start", [0, 1])
def test_mixed_faces_edges_and_steps(tmp_path, reader, transposed, start):
    path = tmp_path / "m.nc"
    _mixed_file(path, transposed, start)
    mesh = reader(path)
    assert [c.type for c in mesh.cells] == ["triangle", "quad", "line"]
    assert np.asarray(mesh.cells[0].data).tolist() == [[1, 2, 5], [1, 5, 4]]
    assert np.asarray(mesh.cells[1].data).tolist() == [[0, 1, 4, 3]]
    assert np.asarray(mesh.cells[2].data).tolist() == [[0, 3], [2, 5]]
    assert np.allclose(mesh.point_data["waterlevel"], np.arange(6.0))
    vel = mesh.cell_data["velocity"]
    assert np.allclose(vel[0], [2.0, 3.0]) and np.allclose(vel[1], [1.0])
    assert np.isnan(vel[2]).all()
    dis = mesh.cell_data["discharge"]
    assert np.isnan(dis[0]).all() and np.allclose(dis[2], [7.0, 8.0])
    assert np.allclose(mesh.field_data["meshio:time"], [10.0])
    last = reader(path, time_step=-1)
    assert np.allclose(last.point_data["waterlevel"], np.arange(6.0) + 100)
    assert np.allclose(last.field_data["meshio:time"], [20.0])
    with pytest.raises(mio.ReadError, match="out of range"):
        reader(path, time_step=5)


def test_time_values(tmp_path):
    path = tmp_path / "m.nc"
    _mixed_file(path)
    assert mio.ugrid_nc.time_values(path) == [10.0, 20.0]
    assert _ugrid_nc.time_values(path) == [10.0, 20.0]


def test_an_exodus_file_is_named(tmp_path, reader):
    path = tmp_path / "e.nc"
    with netCDF4.Dataset(path, "w") as nc:
        nc.createDimension("num_nodes", 3)
        nc.createDimension("num_dim", 3)
    with pytest.raises(mio.ReadError, match="exodus"):
        reader(path)


@pytest.mark.parametrize("writer", ["core", "python"])
def test_round_trip(tmp_path, writer):
    mesh = mio.Mesh(
        np.array(
            [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0.5, 0]], dtype=float
        ),
        [
            ("triangle", np.array([[1, 4, 2]])),
            ("quad", np.array([[0, 1, 2, 3]])),
            ("line", np.array([[0, 3]])),
        ],
        point_data={"h": np.arange(5.0), "v": np.ones((5, 2))},
        cell_data={"c": [np.array([1.0]), np.array([2.0]), np.array([3.0])]},
    )
    path = tmp_path / "w.nc"
    (mio.ugrid_nc.write if writer == "core" else _ugrid_nc.write)(path, mesh)
    for read in (mio.ugrid_nc.read, _ugrid_nc.read):
        back = read(path)
        assert [c.type for c in back.cells] == ["triangle", "quad", "line"]
        assert np.allclose(back.points, mesh.points)
        assert np.allclose(back.point_data["h"], mesh.point_data["h"])
        assert np.allclose(back.point_data["v"], mesh.point_data["v"])
        assert [np.asarray(a).tolist() for a in back.cell_data["c"]] == [
            [1.0],
            [2.0],
            [3.0],
        ]


def test_volume_cells_are_refused(tmp_path):
    mesh = mio.Mesh(np.eye(4, 3), [("tetra", np.array([[0, 1, 2, 3]]))])
    for write in (mio.ugrid_nc.write, _ugrid_nc.write):
        with pytest.raises(mio.WriteError, match="tetra"):
            write(tmp_path / "t.nc", mesh)


def test_xugrid_reads_what_we_write_and_the_reverse(tmp_path):
    xu = pytest.importorskip("xugrid")
    mesh = mio.Mesh(
        np.array([[0, 0], [1, 0], [1, 1], [0, 1], [2, 0.5]], dtype=float),
        [("triangle", np.array([[1, 4, 2]])), ("quad", np.array([[0, 1, 2, 3]]))],
        cell_data={"c": [np.array([1.0]), np.array([2.0])]},
    )
    ours = tmp_path / "ours.nc"
    mio.write(ours, mesh)
    ds = xu.open_dataset(ours)
    grid = ds.ugrid.grid
    assert grid.n_node == 5 and grid.n_face == 2
    assert np.allclose(np.sort(ds["c"].values), [1.0, 2.0])
    theirs = tmp_path / "theirs.nc"
    ds.ugrid.to_netcdf(theirs)
    for read in (mio.ugrid_nc.read, _ugrid_nc.read):
        back = read(theirs)
        assert len(back.points) == 5
        # xugrid also writes the faces' edges, which read back as a line block.
        assert [c.type for c in back.cells][:2] == ["triangle", "quad"]


def test_a_coordinate_longer_than_the_node_dimension_is_refused(tmp_path, reader):
    """fuzz finding: the variable was read whole into a buffer sized by the node
    dimension, past its end (tests/fuzz/regressions/ugrid_nc)."""
    path = str(tmp_path / "long.nc")
    with netCDF4.Dataset(path, "w") as nc:
        nc.createDimension("nNodes", 3)
        nc.createDimension("nLong", 8)
        topo = nc.createVariable("Mesh2", "i4")
        topo.cf_role = "mesh_topology"
        topo.topology_dimension = 1
        topo.node_coordinates = "x y"
        nc.createVariable("x", "f8", ("nNodes",))[:] = [0.0, 1.0, 2.0]
        nc.createVariable("y", "f8", ("nLong",))[:] = np.arange(8.0)
    with pytest.raises(mio.ReadError, match="holds 8 values"):
        reader(path)


def test_a_time_variable_longer_than_its_dimension_is_refused(tmp_path, monkeypatch):
    """The core reads the variable named after the time dimension whole: a longer one
    must be refused, not read past the step buffer."""
    if not getattr(mio._core, "__has_netcdf__", False):
        pytest.skip("this build has no netCDF")
    monkeypatch.setattr(_fallback, "_strict", True)
    path = str(tmp_path / "time.nc")
    with netCDF4.Dataset(path, "w") as nc:
        nc.createDimension("nNodes", 2)
        nc.createDimension("time", 2)
        nc.createDimension("longer", 5)
        topo = nc.createVariable("Mesh2", "i4")
        topo.cf_role = "mesh_topology"
        topo.topology_dimension = 1
        topo.node_coordinates = "x y"
        nc.createVariable("x", "f8", ("nNodes",))[:] = [0.0, 1.0]
        nc.createVariable("y", "f8", ("nNodes",))[:] = [0.0, 1.0]
        nc.createVariable("time", "f8", ("longer",))[:] = np.arange(5.0)
        h = nc.createVariable("h", "f8", ("time", "nNodes"))
        h.mesh = "Mesh2"
        h.location = "node"
        h[:] = np.zeros((2, 2))
    with pytest.raises(mio.ReadError, match="holds 5 values"):
        mio.ugrid_nc.read(path)
