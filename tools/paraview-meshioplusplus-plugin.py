import sys

import numpy as np
from paraview.util.vtkAlgorithm import (
    VTKPythonAlgorithmBase,
    smdomain,
    smhint,
    smproperty,
    smproxy,
)
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

import meshioplusplus

paraview_plugin_version = meshioplusplus.__version__
vtk_to_meshio_type = meshioplusplus._vtk_common.vtk_to_meshio_type
meshio_to_vtk_type = meshioplusplus._vtk_common.meshio_to_vtk_type
meshio_input_filetypes = list(meshioplusplus._helpers.reader_map.keys())
meshio_extensions = [ext[1:] for ext in meshioplusplus.extension_to_filetypes.keys()]
meshio_input_filetypes = ["automatic"] + meshio_input_filetypes


@smproxy.reader(
    name="meshio++ reader",
    extensions=meshio_extensions,
    file_description="meshio++-supported files",
    support_reload=False,
)
class MeshioReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkUnstructuredGrid"
        )
        self._filename = None
        self._file_format = None

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(
        extensions=meshio_extensions, file_description="meshio++-supported files"
    )
    def SetFileName(self, filename):
        if self._filename != filename:
            self._filename = filename
            self.Modified()

    @smproperty.stringvector(name="StringInfo", information_only="1")
    def GetStrings(self):
        return meshio_input_filetypes

    @smproperty.stringvector(name="FileFormat", number_of_elements="1")
    @smdomain.xml(
        """
        <StringListDomain name="list">
            <RequiredProperties>
                <Property name="StringInfo" function="StringInfo"/>
            </RequiredProperties>
        </StringListDomain>
        """
    )
    def SetFileFormat(self, file_format):
        # Automatically deduce input format
        if file_format == "automatic":
            file_format = None

        if self._file_format != file_format:
            self._file_format = file_format
            self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        output = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(outInfoVec))

        # Use meshio++ to read the mesh
        mesh = meshioplusplus.read(self._filename, self._file_format)
        points, cells = mesh.points, mesh.cells

        # Points
        if points.shape[1] == 2:
            points = np.hstack([points, np.zeros((len(points), 1))])
        output.SetPoints(points)

        # Cells, in VTK's legacy (count, ids...) layout. Ragged polygon rows
        # are written one by one; a block VTK has no id for (polyhedra, which
        # need a face stream) is skipped with a message rather than failing
        # the whole read.
        cell_types, cell_offsets, cell_conn, kept = [], [], [], []
        size = 0
        for block in cells:
            vtk_type = meshio_to_vtk_type.get(block.type)
            if vtk_type is None or block.type.startswith("polyhedron"):
                print(
                    f"meshio++ reader: skipping '{block.type}' cells", file=sys.stderr
                )
                kept.append(False)
                continue
            kept.append(True)
            for row in block.data:
                row = np.asarray(row, dtype=np.int64)
                cell_types.append(vtk_type)
                cell_offsets.append(size)
                cell_conn.append(len(row))
                cell_conn.extend(row.tolist())
                size += 1 + len(row)
        output.SetCells(
            np.asarray(cell_types, dtype=np.ubyte),
            np.asarray(cell_offsets, dtype=np.int64),
            np.asarray(cell_conn, dtype=np.int64),
        )

        # Point data
        for name, array in mesh.point_data.items():
            output.PointData.append(array, name)

        # Cell data, for the blocks that became cells
        for name, data in mesh.cell_data.items():
            parts = [np.asarray(d) for d, k in zip(data, kept) if k]
            if parts:
                output.CellData.append(np.concatenate(parts), name)

        # Field data
        for name, array in mesh.field_data.items():
            output.FieldData.append(array, name)

        return 1


@smproxy.writer(
    name="meshio++ Writer",
    extensions=meshio_extensions,
    file_description="meshio++-supported files",
    support_reload=False,
)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
class MeshioWriter(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=1, nOutputPorts=0, inputType="vtkUnstructuredGrid"
        )
        self._filename = None

    @smproperty.stringvector(name="FileName", panel_visibility="never")
    @smdomain.filelist()
    def SetFileName(self, filename):
        if self._filename != filename:
            self._filename = filename
            self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        mesh = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(inInfoVec[0]))

        # Read points
        points = np.asarray(mesh.GetPoints())

        # Cells, grouped by (VTK type, node count) in first-seen order: a
        # meshio++ cell block is one type with one row length.
        grid = mesh.VTKObject
        types = (
            vtk_to_numpy(grid.GetCellTypesArray()) if grid.GetNumberOfCells() else []
        )
        array = grid.GetCells()
        offsets = vtk_to_numpy(array.GetOffsetsArray())
        conn = vtk_to_numpy(array.GetConnectivityArray())
        groups = {}
        for c, vtk_type in enumerate(types):
            row = conn[offsets[c] : offsets[c + 1]]
            groups.setdefault((int(vtk_type), len(row)), []).append((c, row))
        cells, members = [], []
        for (vtk_type, _), rows in groups.items():
            cells.append(
                meshioplusplus.CellBlock(
                    vtk_to_meshio_type[vtk_type], np.array([r for _, r in rows])
                )
            )
            members.append(np.array([c for c, _ in rows]))

        # Read point and field data
        # Adapted from test/legacy_reader.py
        def _read_data(data):
            out = {}
            for i in range(data.VTKObject.GetNumberOfArrays()):
                name = data.VTKObject.GetArrayName(i)
                array = np.asarray(data.GetArray(i))
                out[name] = array
            return out

        point_data = _read_data(mesh.GetPointData())
        field_data = _read_data(mesh.GetFieldData())

        # Cell data, split the same way as the cells
        cell_data = {
            name: [array[m] for m in members]
            for name, array in _read_data(mesh.GetCellData()).items()
        }

        meshioplusplus.write(
            self._filename,
            meshioplusplus.Mesh(
                points,
                cells,
                point_data=point_data,
                cell_data=cell_data,
                field_data=field_data,
            ),
        )
        return 1

    def Write(self):
        self.Modified()
        self.Update()
