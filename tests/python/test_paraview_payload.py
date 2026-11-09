"""The ParaView plugin's conversions (``meshioplusplus._paraview``).

The plugin itself only declares proxies; these tests cover what it does to a
mesh, with the PyPI ``vtk`` wheel standing in for ParaView (the tests skip
without it). ``tests/python/test_paraview_plugin.py`` runs the plugin under a
real ``pvpython``.
"""

import pathlib
import warnings

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _paraview

vtk = pytest.importorskip("vtk")

MESHES = pathlib.Path(__file__).resolve().parent / "meshes"


def _wedge_mesh():
    """One unit wedge, positive in meshio++'s (Gmsh) node order, which is
    VTK's (tests/python/test_vtk_wedge_orientation.py)."""
    points = np.array(
        [[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]]
    )
    return meshioplusplus.Mesh(points, [("wedge", [[0, 1, 2, 3, 4, 5]])])


def _vtk_volumes(grid):
    """Signed cell volumes as VTK computes them."""
    from vtkmodules.util.numpy_support import vtk_to_numpy
    from vtkmodules.vtkFiltersVerdict import vtkCellSizeFilter

    sizes = vtkCellSizeFilter()
    sizes.SetInputData(grid)
    sizes.Update()
    return vtk_to_numpy(sizes.GetOutput().GetCellData().GetArray("Volume"))


def _num_inverted(mesh):
    return meshioplusplus.compute_quality(mesh)["num_inverted"]


@pytest.mark.parametrize("source", ["unit", "toroid"])
def test_reader_wedges_are_not_inside_out(source):
    """roadmap 1.1: a wedge positive to meshio++ is positive to VTK."""
    if source == "unit":
        mesh = _wedge_mesh()
    else:
        mesh = meshioplusplus.read(MESHES / "mfem" / "toroid-wedge.mesh")
    assert _num_inverted(mesh) == 0
    volumes = _vtk_volumes(_paraview.grid_from_mesh(mesh))
    assert len(volumes) == sum(len(b) for b in mesh.cells)
    solids = volumes[volumes != 0]  # the toroid's quadrilateral boundary has none
    assert len(solids) and (solids > 0).all()


def test_writer_wedges_are_not_inverted():
    """roadmap 1.1: a positive VTK wedge is written as a positive meshio++ one."""
    from vtkmodules.vtkCommonDataModel import VTK_WEDGE, vtkUnstructuredGrid

    mesh = _wedge_mesh()
    grid = vtkUnstructuredGrid()
    grid.SetPoints(_paraview.grid_from_mesh(mesh).GetPoints())
    grid.InsertNextCell(VTK_WEDGE, 6, list(range(6)))
    assert _vtk_volumes(grid)[0] > 0
    back = _paraview.mesh_from_grid(grid)
    assert back.cells[0].type == "wedge"
    np.testing.assert_array_equal(back.cells[0].data, [[0, 1, 2, 3, 4, 5]])
    assert _num_inverted(back) == 0


def test_polygon_n_blocks_are_kept():
    """roadmap 1.2: a uniform n-gon block is named polygon5, and is a polygon."""
    points = np.array(
        [[0.0, 0, 0], [1, 0, 0], [1.5, 1, 0], [0.5, 1.6, 0], [-0.5, 1, 0], [3, 0, 0]]
    )
    mesh = meshioplusplus.Mesh(
        points,
        [("polygon5", [[0, 1, 2, 3, 4]]), ("triangle", [[1, 5, 2]])],
        cell_data={"c": [np.array([7.0]), np.array([8.0])]},
    )
    grid = _paraview.grid_from_mesh(mesh)
    assert grid.GetNumberOfCells() == 2
    assert grid.GetCellType(0) == vtk.VTK_POLYGON
    assert grid.GetCell(0).GetNumberOfPoints() == 5
    from vtkmodules.util.numpy_support import vtk_to_numpy

    np.testing.assert_array_equal(
        vtk_to_numpy(grid.GetCellData().GetArray("c")), [7.0, 8.0]
    )


def test_one_dimensional_points_are_padded():
    """roadmap 1.2: a 1-D mesh reaches VTK with three point columns."""
    mesh = meshioplusplus.Mesh(
        np.array([[0.0], [1.0], [3.0]]), [("line", [[0, 1], [1, 2]])]
    )
    grid = _paraview.grid_from_mesh(mesh)
    assert grid.GetNumberOfPoints() == 3
    assert grid.GetNumberOfCells() == 2
    assert grid.GetPoint(2) == (3.0, 0.0, 0.0)


def test_round_trip_keeps_cell_order_and_data():
    """The writer keeps the grid's cell order: runs become blocks."""
    points = np.array(
        [
            [0.0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
        ]
    )
    mesh = meshioplusplus.Mesh(
        points,
        [
            ("tetra", [[0, 1, 2, 3]]),
            ("wedge", [[0, 1, 2, 3, 4, 5]]),
            ("tetra", [[1, 2, 6, 4], [3, 4, 5, 6]]),
        ],
        point_data={"u": np.arange(7.0), "v": np.ones((7, 3), dtype=np.float32)},
        cell_data={"c": [np.array([1]), np.array([2]), np.array([3, 4])]},
        field_data={"f": np.array([1.5, 2.5])},
    )
    back = _paraview.mesh_from_grid(_paraview.grid_from_mesh(mesh))
    assert [(b.type, len(b)) for b in back.cells] == [
        ("tetra", 1),
        ("wedge", 1),
        ("tetra", 2),
    ]
    for got, want in zip(back.cells, mesh.cells):
        np.testing.assert_array_equal(got.data, want.data)
    np.testing.assert_array_equal(back.points, points)
    np.testing.assert_array_equal(back.point_data["u"], np.arange(7.0))
    assert back.point_data["v"].dtype == np.float32
    assert back.point_data["v"].shape == (7, 3)
    assert [list(c) for c in back.cell_data["c"]] == [[1], [2], [3, 4]]
    np.testing.assert_array_equal(back.field_data["f"], [1.5, 2.5])


def test_string_arrays_are_skipped_with_a_warning():
    """A string point/cell array has no meshio++ counterpart; field strings stay."""
    from vtkmodules.vtkCommonCore import vtkStringArray

    grid = _paraview.grid_from_mesh(_wedge_mesh())
    for attributes, name in (
        (grid.GetCellData(), "labels"),
        (grid.GetFieldData(), "QA Records"),
    ):
        strings = vtkStringArray()
        strings.SetName(name)
        strings.InsertNextValue("abc")
        attributes.AddArray(strings)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        back = _paraview.mesh_from_grid(grid)
    assert "labels" not in back.cell_data
    assert list(back.field_data["QA Records"]) == ["abc"]


def test_field_strings_reach_the_grid():
    mesh = _wedge_mesh()
    mesh.field_data["title"] = np.array(["a mesh"])
    grid = _paraview.grid_from_mesh(mesh)
    array = grid.GetFieldData().GetAbstractArray("title")
    assert array.GetValue(0) == "a mesh"


def test_empty_mesh():
    mesh = meshioplusplus.Mesh(np.zeros((0, 3)), [])
    grid = _paraview.grid_from_mesh(mesh)
    assert grid.GetNumberOfCells() == 0
    back = _paraview.mesh_from_grid(grid)
    assert len(back.cells) == 0


# --------------------------------------------------------------------------- #
# polyhedra (roadmap 6.2)                                                     #
# --------------------------------------------------------------------------- #
#: The tetrahedron (0, 1, 3, 4) of the cube below, faces outward.
_TET_FACES = [[0, 3, 1], [0, 1, 4], [1, 3, 4], [0, 4, 3]]
#: A unit cube's six faces, outward.
_CUBE_FACES = [
    [0, 3, 2, 1],
    [4, 5, 6, 7],
    [0, 1, 5, 4],
    [1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 0, 4, 7],
]


def _polyhedra_mesh():
    """Polyhedra of 4 and 8 points around a tetrahedron, with cell data."""
    points = np.array(
        [
            [0.0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ]
    )
    tet = [np.array(f) for f in _TET_FACES]
    cube = [np.array(f) for f in _CUBE_FACES]
    return meshioplusplus.Mesh(
        points,
        [
            ("polyhedron4", [tet, tet]),
            ("tetra", [[1, 2, 3, 6]]),
            ("polyhedron8", [cube]),
        ],
        cell_data={"c": [np.array([1.0, 2.0]), np.array([3.0]), np.array([4.0])]},
    )


def _assert_same_cells(back, mesh):
    assert [(b.type, len(b)) for b in back.cells] == [
        (b.type, len(b)) for b in mesh.cells
    ]
    for got, want in zip(back.cells, mesh.cells):
        if want.type.startswith("polyhedron"):
            assert [[f.tolist() for f in c] for c in got.data] == [
                [np.asarray(f).tolist() for f in c] for c in want.data
            ]
        else:
            np.testing.assert_array_equal(got.data, want.data)


@pytest.mark.parametrize("legacy", [False, True], ids=["vtk94", "legacy"])
def test_polyhedra_round_trip(legacy):
    """Both of VTK's polyhedron APIs: SetPolyhedralCells (9.4 and later) and
    the face stream SetCells took before."""
    from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

    mesh = _polyhedra_mesh()
    if legacy:
        from meshioplusplus._interop import _to_vtk_payload

        grid = vtkUnstructuredGrid()
        grid.SetPoints(_paraview.grid_from_mesh(mesh).GetPoints())
        payload = _to_vtk_payload(mesh)
        try:
            _paraview.set_cells(grid, payload, legacy=True)
        except TypeError:
            pytest.skip("this VTK has no legacy face-stream SetCells")
        grid.GetCellData().AddArray(_paraview._vtk_array(payload.cell_data["c"], "c"))
    else:
        grid = _paraview.grid_from_mesh(mesh)
    assert grid.GetNumberOfCells() == 4
    volumes = _vtk_volumes(grid)
    np.testing.assert_allclose(volumes, [1 / 6, 1 / 6, 1 / 6, 1.0])
    back = _paraview.mesh_from_grid(grid)
    _assert_same_cells(back, mesh)
    assert [c.tolist() for c in back.cell_data["c"]] == [[1.0, 2.0], [3.0], [4.0]]


def test_legacy_face_stream_round_trip():
    from meshioplusplus._interop import _to_vtk_payload

    payload = _to_vtk_payload(_polyhedra_mesh())
    stream, where = _paraview.legacy_face_stream(
        payload.face_connectivity, payload.face_offsets, payload.face_location_offsets
    )
    tet = [4] + sum(([3] + f for f in _TET_FACES), [])
    assert stream[:17].tolist() == tet
    assert stream[17:34].tolist() == tet
    assert where.tolist() == [0, 17, -1, 34]
    assert stream[34] == 6 and stream[35] == 4
    flat, offsets, locations = _paraview.faces_from_legacy_stream(stream, where)
    np.testing.assert_array_equal(flat, payload.face_connectivity)
    np.testing.assert_array_equal(offsets, payload.face_offsets)
    np.testing.assert_array_equal(locations, payload.face_location_offsets)


def test_polyhedron_fixture_keeps_its_cells():
    """The roadmap's probe: a two-piece VTU of polyhedra has the cell count
    meshioplusplus.read gives once through the plugin's conversions."""
    path = MESHES / "vtu" / "two_pieces_polyhedron.vtu"
    if path.read_bytes().startswith(b"version https://git-lfs"):
        pytest.skip("LFS content not fetched")
    mesh = meshioplusplus.read(path)
    grid = _paraview.grid_from_mesh(mesh)
    assert grid.GetNumberOfCells() == sum(len(b) for b in mesh.cells)
    back = _paraview.mesh_from_grid(grid)
    assert sum(len(b) for b in back.cells) == sum(len(b) for b in mesh.cells)


# --------------------------------------------------------------------------- #
# higher order (roadmap 6.3.1)                                                #
# --------------------------------------------------------------------------- #
def test_degrees_and_weights_are_the_active_attributes():
    """HigherOrderDegrees and RationalWeights reach VTK as the attributes its
    Lagrange and Bézier cells read, and come back as arrays."""
    # An order-2 Bézier triangle: corners, then the mid-edge control points.
    points = np.array(
        [
            [0.0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0.5, -0.1, 0],
            [0.6, 0.6, 0],
            [-0.1, 0.5, 0],
        ]
    )
    mesh = meshioplusplus.Mesh(
        points,
        [("VTK_BEZIER_TRIANGLE", [[0, 1, 2, 3, 4, 5]])],
        point_data={"RationalWeights": np.array([1.0, 1, 1, 0.8, 0.8, 0.8])},
        cell_data={"HigherOrderDegrees": [np.array([[2.0, 2, 0]])]},
    )
    grid = _paraview.grid_from_mesh(mesh)
    assert grid.GetCellType(0) == vtk.VTK_BEZIER_TRIANGLE
    degrees = grid.GetCellData().GetHigherOrderDegrees()
    weights = grid.GetPointData().GetRationalWeights()
    assert degrees is not None and degrees.GetName() == "HigherOrderDegrees"
    assert weights is not None and weights.GetName() == "RationalWeights"
    back = _paraview.mesh_from_grid(grid)
    np.testing.assert_array_equal(
        back.point_data["RationalWeights"], [1, 1, 1, 0.8, 0.8, 0.8]
    )
    np.testing.assert_array_equal(back.cell_data["HigherOrderDegrees"][0], [[2, 2, 0]])


# --------------------------------------------------------------------------- #
# time (roadmap 6.4)                                                          #
# --------------------------------------------------------------------------- #
_TIMES = [0.0, 0.5, 2.0]


def _step_meshes():
    points = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    return [
        (
            t,
            meshioplusplus.Mesh(
                points,
                [("tetra", [[0, 1, 2, 3]])],
                point_data={"u": np.full(4, 10.0 * t)},
            ),
        )
        for t in _TIMES
    ]


def test_time_plan_of_a_multi_step_file(tmp_path):
    """The roadmap's probe, without ParaView: a three-step frd has three
    time steps, and the time 2.0 reads the last step."""
    path = tmp_path / "steps.frd"
    meshioplusplus.write_sequence(path, _step_meshes())
    plan = _paraview.time_plan([path])
    times = _paraview.time_values(plan)
    np.testing.assert_allclose(times, _TIMES)
    index = _paraview.step_index(times, 2.0)
    assert index == 2
    mesh = _paraview.read_step(plan, index)
    np.testing.assert_allclose(mesh.point_data["u"], 20.0)
    assert _paraview.time_plan([path]) is plan  # cached


def test_time_plan_of_a_file_series(tmp_path):
    """A series of files with no number in their names takes each file's own
    time (the meshio:time field data write_sequence stores)."""
    paths = []
    for name, (t, mesh) in zip(["first", "second", "third"], _step_meshes()):
        mesh.field_data["meshio:time"] = np.array([t])
        paths.append(tmp_path / f"{name}.vtu")
        meshioplusplus.write(paths[-1], mesh)
    plan = _paraview.time_plan(paths)
    np.testing.assert_allclose(_paraview.time_values(plan), _TIMES)
    np.testing.assert_allclose(_paraview.read_step(plan, 1).point_data["u"], 5.0)


def test_a_rewritten_file_is_planned_again(tmp_path):
    path = tmp_path / "steps.frd"
    meshioplusplus.write_sequence(path, _step_meshes())
    first = _paraview.time_plan([path])
    meshioplusplus.write_sequence(path, _step_meshes()[:2])
    import os

    os.utime(path, ns=(1, 1))  # a different identity even on a coarse clock
    assert len(_paraview.time_plan([path])) == 2
    assert len(first) == 3


def test_a_single_step_has_no_time():
    assert _paraview.time_values([{"path": "x", "step": 0, "time": 0.0}]) == []


def test_step_index_takes_the_last_step_at_or_before():
    assert _paraview.step_index(_TIMES, -1.0) == 0
    assert _paraview.step_index(_TIMES, 0.5) == 1
    assert _paraview.step_index(_TIMES, 1.9) == 1
    assert _paraview.step_index(_TIMES, 5.0) == 2
    assert _paraview.step_index([], 5.0) == 0


def test_series_path_and_one_file_per_step(tmp_path):
    assert _paraview.series_path(tmp_path / "a.frd") == str(tmp_path / "a.frd")
    template = _paraview.series_path(tmp_path / "a.vtu")
    assert template == str(tmp_path / "a_{step}.vtu")
    for i, (t, mesh) in enumerate(_step_meshes()):
        _paraview.write_one_step(template, i, 3, t, mesh)
    written = sorted(tmp_path.glob("a_*.vtu"))
    assert [p.name for p in written] == ["a_0000.vtu", "a_0001.vtu", "a_0002.vtu"]
    # Each file keeps its time, though a number in a file name takes
    # precedence when the series is read back (the sequence driver's rule).
    for path, t in zip(written, _TIMES):
        assert meshioplusplus.read(path).field_data["meshio:time"][0] == t
    np.testing.assert_allclose(
        _paraview.time_values(_paraview.time_plan(written)), [0, 1, 2]
    )


# --------------------------------------------------------------------------- #
# read options and status (roadmap 6.5)                                       #
# --------------------------------------------------------------------------- #
def test_read_kwargs_pass_only_what_differs():
    assert _paraview.read_kwargs(["u", "v"], ["c"]) == {}
    assert _paraview.read_kwargs(["u", "v"], ["c"], ["u", "v"], ["c"]) == {}
    assert _paraview.read_kwargs(["u", "v"], ["c"], ["u"], ["c"]) == {
        "arrays": ["c", "u"]
    }
    assert _paraview.read_kwargs(["u"], ["c"], [], []) == {"arrays": []}
    assert _paraview.read_kwargs(points_only=True, piece=2, ghosts="drop") == {
        "points_only": True,
        "piece": 2,
        "ghosts": "drop",
    }
    assert _paraview.read_kwargs(piece=-1) == {}
    with pytest.raises(ValueError, match="ghosts"):
        _paraview.read_kwargs(ghosts="only")


def test_an_unticked_array_is_not_read(tmp_path):
    mesh = _step_meshes()[0][1]
    mesh.point_data["v"] = np.arange(4.0)
    path = tmp_path / "m.vtu"
    meshioplusplus.write(path, mesh)
    meta = _paraview.metadata(path)
    assert sorted(meta["point_data_names"]) == ["u", "v"]
    kwargs = _paraview.read_kwargs(meta["point_data_names"], [], ["v"], [])
    plan = _paraview.time_plan([path])
    back = _paraview.read_step(plan, 0, **kwargs)
    assert list(back.point_data) == ["v"]


def test_cad_formats_are_recognised():
    assert _paraview.cad_format("part.STEP") == "step"
    assert _paraview.cad_format("part.igs") == "iges"
    assert _paraview.cad_format("part.vtu") is None
    assert _paraview.cad_format("part.dat", "step") == "step"


def test_addon_status_names_the_three_addons():
    status = _paraview.addon_status()
    assert set(status) == {"occt", "netgen", "mmg"}
    assert all(isinstance(v, bool) for v in status.values())


# --------------------------------------------------------------------------- #
# named regions (roadmap 6.6)                                                 #
# --------------------------------------------------------------------------- #
def _regions_mesh():
    from meshioplusplus._regions import Region

    return meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]),
        [("triangle", [[0, 1, 2]]), ("tetra", [[0, 1, 2, 3], [1, 2, 3, 4]])],
        regions=[
            Region("inlet", "point", [0, 1], 0, 5),
            Region("steel", "cell", [1, 2], 3, 7),
            Region("skin", "cell", [0], 2, 8),
            Region("wall", "side", [[1, 0], [2, 3]], 2, 9),
        ],
    )


def _regions_of(mesh):
    return sorted(
        (r.name, r.kind, np.asarray(r.entries).tolist(), r.dim, r.tag)
        for r in mesh.regions
    )


def test_regions_round_trip_as_arrays():
    """roadmap 6.6.1: point and cell regions as region: masks with the
    sidecar, side regions as the VTU convention's field arrays; all come
    back with their dim and tag, and no carrier array is left over."""
    mesh = _regions_mesh()
    grid = _paraview.grid_from_mesh(mesh)
    cell_names = {
        grid.GetCellData().GetArrayName(i)
        for i in range(grid.GetCellData().GetNumberOfArrays())
    }
    assert cell_names == {"region:steel", "region:skin"}
    assert grid.GetFieldData().GetAbstractArray("region:side:wall") is not None
    back = _paraview.mesh_from_grid(grid)
    assert _regions_of(back) == _regions_of(mesh)
    assert dict(back.field_data) == {} and back.cell_data == {}
    assert back.point_data == {}


def test_regions_from_a_vtu_paraview_opened_itself(tmp_path):
    """A meshio++ VTU read by ParaView's own reader carries its regions as
    region:<kind>:<name> field arrays; the writer turns them back into
    regions."""
    path = tmp_path / "r.vtu"
    meshioplusplus.write(path, _regions_mesh())
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    back = _paraview.mesh_from_grid(reader.GetOutput())
    assert _regions_of(back) == _regions_of(_regions_mesh())


def test_regions_as_blocks(tmp_path):
    """roadmap 6.6.2: the roadmap's probe -- a Gmsh file with two physical
    groups shows both names in the assembly, each with its cell count."""
    from meshioplusplus._regions import Region

    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]),
        [("tetra", [[0, 1, 2, 3], [1, 2, 3, 4]]), ("triangle", [[0, 1, 2]])],
        regions=[
            Region("solid", "cell", [0, 1], 3, 1),
            Region("bottom", "cell", [2], 2, 2),
        ],
    )
    path = tmp_path / "groups.msh"
    meshioplusplus.write(path, mesh)
    read = meshioplusplus.read(path)
    collection = _paraview.fill_collection(vtk.vtkPartitionedDataSetCollection(), read)
    names = {
        collection.GetMetaData(i)
        .Get(vtk.vtkCompositeDataSet.NAME()): collection.GetPartition(i, 0)
        .GetNumberOfCells()
        for i in range(collection.GetNumberOfPartitionedDataSets())
    }
    expected = {r.name: len(r.entries) for r in read.regions if r.kind == "cell"}
    assert names == expected
    assert names["solid"] == 2 and names["bottom"] == 1
    assembly = collection.GetDataAssembly()
    group = assembly.FindFirstNodeWithName("cell_regions")
    children = {assembly.GetNodeName(n) for n in assembly.GetChildNodes(group)}
    assert children == set(expected)


def test_region_pieces_of_a_mesh_without_regions():
    mesh = _step_meshes()[0][1]
    assert [(g, n) for g, n, _ in _paraview.region_pieces(mesh)] == [("mesh", "mesh")]


# --------------------------------------------------------------------------- #
# the writer's inputs and options (roadmap 6.7)                               #
# --------------------------------------------------------------------------- #
def _lz4_available():
    from meshioplusplus import _core

    if getattr(_core, "__has_lz4__", False):
        return True
    try:
        import lz4.block  # noqa: F401
    except ImportError:
        return False
    return True


def test_write_kwargs_go_through_the_pipeline_mapping(tmp_path):
    vtu = tmp_path / "x.vtu"
    assert _paraview.write_kwargs(vtu) == {}
    assert _paraview.write_kwargs(vtu, codec="lz4") == {
        "binary": True,
        "compression": "lz4",
    }
    assert _paraview.write_kwargs(vtu, encoding="ascii", float_format=".6e") == {
        "binary": False,
        "float_fmt": ".6e",
    }
    assert _paraview.write_kwargs(tmp_path / "x.dat", "vtu", "binary") == {
        "binary": True,
        "file_format": "vtu",
    }
    with pytest.raises(ValueError, match="ascii"):
        _paraview.write_kwargs(vtu, encoding="ascii", codec="zlib")


def test_lz4_reaches_the_vtu_header(tmp_path):
    """The roadmap's probe: Codec=lz4 puts vtkLZ4DataCompressor in the header."""
    if not _lz4_available():
        pytest.skip("no lz4 in this build")
    path = tmp_path / "x.vtu"
    meshioplusplus.write(
        path, _wedge_mesh(), **_paraview.write_kwargs(path, codec="lz4")
    )
    assert b"vtkLZ4DataCompressor" in path.read_bytes()[:400]


def test_provenance_names_paraview(tmp_path):
    path = tmp_path / "x.vtu"
    with _paraview.provenance_scope("ParaView 6.1.0", path):
        meshioplusplus.write(path, _wedge_mesh())
    lines = meshioplusplus.read_metadata(path)["provenance"]
    assert any("ParaView 6.1.0" in line for line in lines)
    with _paraview.provenance_scope(None, path):
        meshioplusplus.write(path, _wedge_mesh())
    assert not any(
        "ParaView" in line for line in meshioplusplus.read_metadata(path)["provenance"]
    )


def test_poly_data_input():
    """The roadmap's probe without ParaView: a sphere source (poly data)
    becomes triangles."""
    sphere = vtk.vtkSphereSource()
    sphere.Update()
    mesh = _paraview.mesh_from_data_object(sphere.GetOutput())
    assert [(b.type, len(b)) for b in mesh.cells] == [
        ("triangle", sphere.GetOutput().GetNumberOfCells())
    ]
    assert mesh.regions == []


def _block(offset):
    return _paraview.grid_from_mesh(
        meshioplusplus.Mesh(
            np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]) + offset,
            [("tetra", [[0, 1, 2, 3]])],
            cell_data={"c": [np.array([offset])]},
        )
    )


def test_a_multiblock_becomes_one_mesh_with_a_region_per_block(tmp_path):
    """The roadmap's probe: a two-block multiblock writes .msh with two
    physical groups."""
    blocks = vtk.vtkMultiBlockDataSet()
    for i, name in enumerate(["left", "right"]):
        blocks.SetBlock(i, _block(2.0 * i))
        blocks.GetMetaData(i).Set(vtk.vtkCompositeDataSet.NAME(), name)
    mesh = _paraview.mesh_from_data_object(blocks)
    assert [(b.type, len(b)) for b in mesh.cells] == [("tetra", 1), ("tetra", 1)]
    assert [(r.name, r.kind, list(r.entries)) for r in mesh.regions] == [
        ("left", "cell", [0]),
        ("right", "cell", [1]),
    ]
    assert [list(c) for c in mesh.cell_data["c"]] == [[0.0], [2.0]]
    path = tmp_path / "blocks.msh"
    meshioplusplus.write(path, mesh)
    groups = {
        r.name: len(r.entries)
        for r in meshioplusplus.read(path).regions
        if r.kind == "cell" and r.dim == 3
    }
    assert groups == {"left": 1, "right": 1}
