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
from vtkmodules.vtkCommonDataModel import vtkDataObject, vtkUnstructuredGrid

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
        self._filenames = []
        self._file_format = None

    # Repeatable, so ParaView hands over a whole file series (out_0.vtu,
    # out_1.vtu, ...) as one reader; a script may still pass one string.
    @smproperty.stringvector(
        name="FileName",
        command="AddFileName",
        clean_command="RemoveAllFileNames",
        repeat_command="1",
        number_of_elements_per_command="1",
    )
    @smdomain.filelist()
    @smhint.filechooser(
        extensions=meshio_extensions, file_description="meshio++-supported files"
    )
    def AddFileName(self, filename):
        if filename and filename not in self._filenames:
            self._filenames.append(filename)
            self.Modified()

    def RemoveAllFileNames(self):
        if self._filenames:
            self._filenames = []
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

    @smproperty.xml(
        """<Property name="Refresh" command="Refresh" panel_widget="command_button">
             <Documentation>Plan the file's time steps again.</Documentation>
           </Property>"""
    )
    def Refresh(self):
        _paraview.forget_plans()
        self.Modified()

    def _plan(self):
        if not self._filenames:
            return []
        return _paraview.time_plan(self._filenames, self._file_format)

    @smproperty.doublevector(
        name="TimestepValues",
        information_only="1",
        si_class="vtkSITimeStepsProperty",
    )
    def GetTimestepValues(self):
        return _paraview.time_values(self._plan())

    def RequestInformation(self, request, inInfoVec, outInfoVec):
        executive = self.GetExecutive()
        info = outInfoVec.GetInformationObject(0)
        info.Remove(executive.TIME_STEPS())
        info.Remove(executive.TIME_RANGE())
        times = _paraview.time_values(self._plan())
        for t in times:
            info.Append(executive.TIME_STEPS(), t)
        if times:
            info.Append(executive.TIME_RANGE(), times[0])
            info.Append(executive.TIME_RANGE(), times[-1])
        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        executive = self.GetExecutive()
        info = outInfoVec.GetInformationObject(0)
        plan = self._plan()
        times = _paraview.time_values(plan)
        t = None
        if times and info.Has(executive.UPDATE_TIME_STEP()):
            t = info.Get(executive.UPDATE_TIME_STEP())
        index = _paraview.step_index(times, t)
        mesh = _paraview.read_step(plan, index, self._file_format)
        output = vtkUnstructuredGrid.GetData(outInfoVec)
        _paraview.fill_grid(output, mesh)
        if times:
            output.GetInformation().Set(vtkDataObject.DATA_TIME_STEP(), times[index])
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
        self._all_steps = False
        self._times = []
        self._index = 0
        self._steps = []

    @smproperty.stringvector(name="FileName", panel_visibility="never")
    @smdomain.filelist()
    def SetFileName(self, filename):
        if self._filename != filename:
            self._filename = filename
            self.Modified()

    @smproperty.intvector(name="WriteTimeSteps", default_values="0")
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetWriteTimeSteps(self, value):
        if self._all_steps != bool(value):
            self._all_steps = bool(value)
            self.Modified()

    def RequestInformation(self, request, inInfoVec, outInfoVec):
        executive = self.GetExecutive()
        info = inInfoVec[0].GetInformationObject(0)
        self._times = []
        if self._all_steps and info.Has(executive.TIME_STEPS()):
            self._times = list(info.Get(executive.TIME_STEPS()))
        return 1

    def RequestUpdateExtent(self, request, inInfoVec, outInfoVec):
        if self._times:
            executive = self.GetExecutive()
            info = inInfoVec[0].GetInformationObject(0)
            info.Set(executive.UPDATE_TIME_STEP(), self._times[self._index])
        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        mesh = _paraview.mesh_from_grid(vtkUnstructuredGrid.GetData(inInfoVec[0]))
        if len(self._times) < 2:
            meshioplusplus.write(self._filename, mesh)
            return 1
        # Every step, one pass of the pipeline each (CONTINUE_EXECUTING): a
        # one-file-per-step format is written as each step arrives; a format
        # holding a series in one file is written once the last one has.
        target = _paraview.series_path(self._filename)
        t, n = self._times[self._index], len(self._times)
        if target == str(self._filename):
            self._steps.append((t, mesh))
        else:
            _paraview.write_one_step(target, self._index, n, t, mesh)
        executive = self.GetExecutive()
        if self._index + 1 < n:
            self._index += 1
            request.Set(executive.CONTINUE_EXECUTING(), 1)
            return 1
        request.Remove(executive.CONTINUE_EXECUTING())
        if self._steps:
            meshioplusplus.write_sequence(self._filename, self._steps)
        self._index, self._steps = 0, []
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
