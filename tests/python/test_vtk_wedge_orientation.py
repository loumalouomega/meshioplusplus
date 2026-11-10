"""Linear wedges keep their orientation through VTK's formats (roadmap 1.1).

meshio++ stores a ``wedge`` in Gmsh's prism order. VTK's own ``vtkWedge``
uses the same orientation: its reference element is (0,0,0), (1,0,0), (0,1,0)
under (0,0,1), (1,0,1), (0,1,1), its face table lists the base as ``[0, 2, 1]``
(outward), and its Jacobian and the verdict volume ParaView's Mesh Quality
reports are positive for it. meshio++ used to permute the wedge with
``[0, 2, 1, 3, 5, 4]`` on the way to VTK, which turned every wedge inside out.
These tests hold the sign to VTK itself, so they skip without the ``vtk``
wheel.
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus.vtk import _main as _vtk_py
from meshioplusplus.vtu import _vtu

vtk = pytest.importorskip("vtk")
from vtkmodules.util.numpy_support import vtk_to_numpy  # noqa: E402

_POINTS = np.array(
    [[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 0]]
)


def _wedge_mesh():
    """A wedge positive to meshio++, beside a tetrahedron as a control."""
    return meshioplusplus.Mesh(
        _POINTS,
        [("wedge", [[0, 1, 2, 3, 4, 5]]), ("tetra", [[1, 6, 2, 4]])],
    )


def _vtk_sizes(grid):
    """Signed cell volumes as VTK computes them."""
    sizes = vtk.vtkCellSizeFilter()
    sizes.SetInputData(grid)
    sizes.Update()
    return vtk_to_numpy(sizes.GetOutput().GetCellData().GetArray("Volume"))


def _vtk_read(path):
    if path.suffix == ".vtu":
        reader = vtk.vtkXMLUnstructuredGridReader()
    else:
        reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


def test_meshio_wedge_is_positive_to_meshio():
    assert meshioplusplus.compute_quality(_wedge_mesh())["num_inverted"] == 0


def test_vtk_reference_wedge_is_positive():
    """VTK's own wedge, nodes in VTK's order, has a positive volume: the
    convention the tests below rely on."""
    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    for p in _POINTS[:6]:
        points.InsertNextPoint(*p)
    grid.SetPoints(points)
    grid.InsertNextCell(vtk.VTK_WEDGE, 6, list(range(6)))
    assert _vtk_sizes(grid)[0] > 0


@pytest.mark.parametrize(
    "suffix, writer",
    [
        (".vtu", meshioplusplus.write),
        (".vtu", _vtu.write),
        (".vtk", meshioplusplus.write),
        (".vtk", _vtk_py.write),
    ],
    ids=["vtu-core", "vtu-python", "vtk-core", "vtk-python"],
)
def test_written_wedges_are_positive_in_vtk(tmp_path, suffix, writer):
    path = tmp_path / f"wedge{suffix}"
    writer(path, _wedge_mesh())
    sizes = _vtk_sizes(_vtk_read(path))
    assert len(sizes) == 2
    assert (sizes > 0).all(), sizes


@pytest.mark.parametrize("suffix", [".vtu", ".vtk"])
def test_vtk_made_wedges_read_back_positive(tmp_path, suffix):
    """A wedge VTK wrote is not inverted once meshio++ has read it."""
    grid = _vtk_read_back(tmp_path, suffix)
    assert _vtk_sizes(grid)[0] > 0
    for reader in (meshioplusplus.read, _vtu.read if suffix == ".vtu" else None):
        if reader is None:
            continue
        mesh = reader(tmp_path / f"vtk{suffix}")
        assert [b.type for b in mesh.cells] == ["wedge"]
        np.testing.assert_array_equal(mesh.cells[0].data, [[0, 1, 2, 3, 4, 5]])
        assert meshioplusplus.compute_quality(mesh)["num_inverted"] == 0


def _vtk_read_back(tmp_path, suffix):
    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    for p in _POINTS[:6]:
        points.InsertNextPoint(*p)
    grid.SetPoints(points)
    grid.InsertNextCell(vtk.VTK_WEDGE, 6, list(range(6)))
    if suffix == ".vtu":
        writer = vtk.vtkXMLUnstructuredGridWriter()
    else:
        writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(str(tmp_path / f"vtk{suffix}"))
    writer.SetInputData(grid)
    writer.Write()
    return grid
