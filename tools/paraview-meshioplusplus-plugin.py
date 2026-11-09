"""meshio++ in ParaView: a reader and a writer of unstructured grids.

Only the proxy declarations live here; every conversion is in
``meshioplusplus._paraview``, which is tested without ParaView. See
doc/paraview_plugin.md.
"""

from paraview.util.vtkAlgorithm import (
    VTKPythonAlgorithmBase,
    smdomain,
    smhint,
    smproperty,
    smproxy,
)
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

import meshioplusplus
from meshioplusplus import _paraview

paraview_plugin_version = meshioplusplus.__version__
_formats = meshioplusplus.formats()
meshio_extensions = [ext[1:] for ext in _formats["extensions"]]
meshio_input_filetypes = ["automatic"] + _formats["readable"]


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
        mesh = meshioplusplus.read(self._filename, self._file_format)
        _paraview.fill_grid(vtkUnstructuredGrid.GetData(outInfoVec), mesh)
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
        grid = vtkUnstructuredGrid.GetData(inInfoVec[0])
        meshioplusplus.write(self._filename, _paraview.mesh_from_grid(grid))
        return 1

    def Write(self):
        self.Modified()
        self.Update()
