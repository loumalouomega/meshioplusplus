"""meshio++ in ParaView: a reader and a writer of unstructured grids.

Only the proxy declarations live here; every conversion is in
``meshioplusplus._paraview``, which is tested without ParaView. See
doc/paraview_plugin.md.

Proxy names are stable identifiers (``MeshioPlusPlusReader``,
``MeshioPlusPlusWriter``); ``paraview.simple`` names its functions after the
labels, which keep the spellings scripts already use (``meshioreader``,
``meshioWriter``). The names the plugin used before are registered as aliases,
so state files saved with them still load.
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


def _alias(proxy, name, decorator):
    """Register ``proxy``'s class a second time under its old proxy ``name``.

    ParaView collects a proxy's properties from its class's own ``__dict__``,
    so the alias copies the decorated methods rather than only inheriting
    them; it drops the reader/writer factory hints, so the file dialogs offer
    one meshio++ reader and one writer.
    """
    base = proxy._pv_original_func
    namespace = {k: v for k, v in vars(base).items() if not k.startswith("__")}
    namespace["_pvsm_hints_xmls"] = []
    alias = type(base.__name__ + "OldName", (base,), namespace)
    globals()[alias.__name__] = decorator(
        name=name, label=name + " (old name)", support_reload=False
    )(alias)


@smproxy.reader(
    name="MeshioPlusPlusReader",
    label="meshio++ reader",
    extensions=meshio_extensions,
    file_description="meshio++-supported files",
    support_reload=False,
)
class MeshioPlusPlusReader(VTKPythonAlgorithmBase):
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
    name="MeshioPlusPlusWriter",
    label="meshio++ Writer",
    extensions=meshio_extensions,
    file_description="meshio++-supported files",
    support_reload=False,
)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
class MeshioPlusPlusWriter(VTKPythonAlgorithmBase):
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


_alias(MeshioPlusPlusReader, "meshio++ reader", smproxy.source)
_alias(
    MeshioPlusPlusWriter,
    "meshio++ Writer",
    lambda **kw: smproxy.source(group="writers", type="WriterProxy", **kw),
)
