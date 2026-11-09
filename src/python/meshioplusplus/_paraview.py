"""ParaView bridge: the conversions behind the ParaView plugin.

The plugin (``tools/paraview-meshioplusplus-plugin.py``, installed with the
wheel under ``share/meshioplusplus/paraview/``) declares ParaView proxies and
nothing else; everything it does to a mesh lives here, so it can be tested
without ParaView. See ``doc/paraview_plugin.md``.

Architecture
------------
Split the way :mod:`meshioplusplus._blender` is, in two layers:

* a **payload** layer that imports neither VTK nor ParaView. The mesh-to-VTK
  direction is :func:`meshioplusplus._interop._to_vtk_payload`, the one
  ``to_pyvista`` uses (vectorised, every ``polygon*`` block mapped to
  ``VTK_POLYGON``, points of any width padded to three columns); the
  VTK-to-mesh direction is
  :func:`meshioplusplus._vtk_common.vtk_cells_from_data`, the VTU reader's
  (cells split into runs in their original order). Node orders need no
  permutation: meshio++'s are VTK's. Neither is copied here.
* a **VTK adapter** layer, :func:`grid_from_mesh` and :func:`mesh_from_grid`,
  which imports ``vtkmodules`` inside each function and moves whole arrays
  through ``numpy_support``: no per-cell Python loop in either direction. It
  needs VTK but not ParaView, so the PyPI ``vtk`` wheel tests it.

Notes about what could not be carried (a cell block VTK has no type for, a
string array) are surfaced as warnings, the interop layer's convention.
"""

from __future__ import annotations

import numpy as np

#: Operation names the notes are reported under.
READER_OP = "ParaView reader"
WRITER_OP = "ParaView writer"


# --------------------------------------------------------------------------- #
# VTK adapter: mesh -> vtkUnstructuredGrid                                    #
# --------------------------------------------------------------------------- #
def _vtk_array(values, name=None):
    """A VTK data array sharing ``values`` (kept alive by the array)."""
    from vtkmodules.util.numpy_support import numpy_to_vtk

    values = np.ascontiguousarray(values)
    array = numpy_to_vtk(values, deep=False)
    # numpy_to_vtk attaches the numpy array to the VTK object, which holds the
    # buffer for as long as VTK holds the array.
    if name is not None:
        array.SetName(name)
    return array


def _id_array(values):
    """A ``vtkIdTypeArray`` over int64 ``values``, shared."""
    from vtkmodules.util.numpy_support import numpy_to_vtkIdTypeArray

    return numpy_to_vtkIdTypeArray(
        np.ascontiguousarray(values, dtype=np.int64), deep=False
    )


def _field_array(name, raw, notes):
    """One field-data array: numbers as a data array, strings as a string array."""
    a = np.asarray(raw)
    if a.dtype.kind in "SU" or a.dtype == object:
        from vtkmodules.vtkCommonCore import vtkStringArray

        strings = vtkStringArray()
        strings.SetName(name)
        for value in a.reshape(-1):
            strings.InsertNextValue(str(value))
        return strings
    if np.issubdtype(a.dtype, np.complexfloating):
        notes.append(f"field_data:{name} has dtype {a.dtype} and was dropped")
        return None
    if a.dtype == bool:
        a = a.astype(np.uint8)
    if a.ndim == 0:
        a = a.reshape(1)
    return _vtk_array(a, name)


def fill_grid(grid, mesh, op=READER_OP):
    """Fill the ``vtkUnstructuredGrid`` ``grid`` from ``mesh``, in place.

    The plugin's reader passes the grid ParaView hands it; :func:`grid_from_mesh`
    passes a new one. Returns the payload's notes, already emitted.
    """
    from vtkmodules.util.numpy_support import numpy_to_vtk
    from vtkmodules.vtkCommonCore import VTK_UNSIGNED_CHAR, vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray

    from ._interop import _emit, _to_vtk_payload

    payload = _to_vtk_payload(mesh)
    notes = list(payload.notes)

    points = vtkPoints()
    points.SetData(_vtk_array(payload.points))
    grid.SetPoints(points)

    cells = vtkCellArray()
    cells.SetData(_id_array(payload.offsets), _id_array(payload.connectivity))
    types = numpy_to_vtk(
        np.ascontiguousarray(payload.celltypes, dtype=np.uint8),
        deep=False,
        array_type=VTK_UNSIGNED_CHAR,
    )
    grid.SetCells(types, cells)

    for name, values in payload.point_data.items():
        grid.GetPointData().AddArray(_vtk_array(values, name))
    for name, values in payload.cell_data.items():
        grid.GetCellData().AddArray(_vtk_array(values, name))
    for name in sorted(mesh.field_data):
        array = _field_array(name, mesh.field_data[name], notes)
        if array is not None:
            grid.GetFieldData().AddArray(array)

    _emit(op, notes)
    return notes


def grid_from_mesh(mesh, op=READER_OP):
    """A new ``vtkUnstructuredGrid`` holding ``mesh``."""
    from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

    grid = vtkUnstructuredGrid()
    fill_grid(grid, mesh, op)
    return grid


# --------------------------------------------------------------------------- #
# VTK adapter: vtkUnstructuredGrid -> mesh                                    #
# --------------------------------------------------------------------------- #
def _numpy(array):
    """A VTK data array as numpy, through a copy into AOS storage if needed.

    ``vtk_to_numpy`` needs a contiguous array-of-structures buffer; VTK can hand
    back other storage (an implicit or affine array, say) for offsets or
    connectivity, which a deep copy into a plain array of the same type
    converts.
    """
    from vtkmodules.util.numpy_support import vtk_to_numpy
    from vtkmodules.vtkCommonCore import vtkDataArray

    try:
        return vtk_to_numpy(array)
    except Exception:
        plain = vtkDataArray.CreateDataArray(array.GetDataType())
        plain.DeepCopy(array)
        return vtk_to_numpy(plain)


def _arrays(attributes, what, notes):
    """``{name: numpy}`` for the numeric arrays of a ``vtkFieldData``.

    String and other non-numeric arrays (the QA records and information lines
    ParaView's Exodus reader adds as field data) have no meshio++ counterpart
    among point and cell data and are skipped with a note; for field data
    a string array is kept as a numpy string array.
    """
    from vtkmodules.vtkCommonCore import vtkDataArray, vtkStringArray

    out = {}
    for i in range(attributes.GetNumberOfArrays()):
        array = attributes.GetAbstractArray(i)
        if array is None:
            continue
        name = array.GetName() or f"{what}_{i}"
        if isinstance(array, vtkDataArray):
            out[name] = _numpy(array)
        elif what == "field_data" and isinstance(array, vtkStringArray):
            out[name] = np.array(
                [array.GetValue(k) for k in range(array.GetNumberOfValues())]
            )
        else:
            notes.append(
                f"{what}:{name} is a {array.GetClassName()} and was skipped "
                "(meshio++ holds numeric arrays only)"
            )
    return out


def mesh_from_grid(grid, op=WRITER_OP):
    """A mesh from a ``vtkUnstructuredGrid``, cells in the grid's order.

    Consecutive cells of one type (and, for polygons, one size) form a cell
    block, so the grid's cell order is kept: a mixed grid becomes several
    blocks rather than one per type.
    """
    from ._interop import _emit
    from ._mesh import Mesh
    from ._vtk_common import vtk_cells_from_data

    notes: list = []
    if grid.GetPoints() is not None and grid.GetNumberOfPoints():
        points = _numpy(grid.GetPoints().GetData())
    else:
        points = np.zeros((0, 3))

    point_data = _arrays(grid.GetPointData(), "point_data", notes)
    cell_flat = _arrays(grid.GetCellData(), "cell_data", notes)
    field_data = _arrays(grid.GetFieldData(), "field_data", notes)

    if grid.GetNumberOfCells():
        types = _numpy(grid.GetCellTypesArray()).astype(np.uint8)
        cell_array = grid.GetCells()
        offsets = _numpy(cell_array.GetOffsetsArray()).astype(np.int64)
        connectivity = _numpy(cell_array.GetConnectivityArray()).astype(np.int64)
        # vtk_cells_from_data takes END offsets (the VTU convention); a VTK 9
        # cell array's offsets start at 0 and have num_cells + 1 entries.
        cells, cell_data = vtk_cells_from_data(
            connectivity, offsets[1:], types, cell_flat
        )
    else:
        cells, cell_data = [], {}

    _emit(op, notes)
    return Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
    )
