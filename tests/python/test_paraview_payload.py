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
