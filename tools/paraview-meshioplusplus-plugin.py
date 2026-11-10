"""meshio++ in ParaView: a reader, a writer and filters.

Only the proxy declarations live here; every conversion is in
``meshioplusplus._paraview``, which is tested without ParaView. See
doc/paraview_plugin.md.

Proxy names are stable identifiers (``MeshioPlusPlusReader``,
``MeshioPlusPlusWriter``); ``paraview.simple`` names its functions after the
labels, which keep the spellings scripts already use (``meshioreader``,
``meshioWriter``). The names the plugin used before are registered as aliases,
so state files saved with them still load.

The self-contained bundle (src/paraview/build_bundle.py) puts this file beside
a site-packages folder holding meshio++; it is put on the path first.
"""

import html
import os
import sys

from paraview.util.vtkAlgorithm import (
    VTKPythonAlgorithmBase,
    smdomain,
    smhint,
    smproperty,
    smproxy,
)
from vtkmodules.vtkCommonCore import vtkDataArraySelection
from vtkmodules.vtkCommonDataModel import (
    vtkDataObject,
    vtkPartitionedDataSet,
    vtkPartitionedDataSetCollection,
    vtkUnstructuredGrid,
)

_BUNDLED = os.path.join(os.path.dirname(os.path.abspath(__file__)), "site-packages")
if os.path.isdir(_BUNDLED) and _BUNDLED not in sys.path:
    sys.path.insert(0, _BUNDLED)

import meshioplusplus  # noqa: E402  isort:skip
from meshioplusplus import _paraview  # noqa: E402  isort:skip

paraview_plugin_version = meshioplusplus.__version__
_formats = meshioplusplus.formats()
meshio_extensions = [ext[1:] for ext in _formats["extensions"]]
meshio_input_filetypes = ["automatic"] + _formats["readable"]
meshio_output_filetypes = ["automatic"] + _formats["writable"]


def _paraview_version():
    """This ParaView's version, for the written file's provenance."""
    try:
        from paraview.modules.vtkRemotingCore import vtkPVVersion
    except ImportError:
        return "(unknown version)"
    return vtkPVVersion.GetPVVersionFull()


def _modified_callback(algorithm):
    """An observer that marks ``algorithm`` modified (an array selection
    changed), holding it weakly so the observer keeps no cycle alive."""
    import weakref

    ref = weakref.ref(algorithm)

    def mark(*args):
        target = ref()
        if target is not None:
            target.Modified()

    return mark


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
        # vtkDataObject: RequestDataObject picks an unstructured grid, or a
        # partitioned dataset collection when RegionsAsBlocks is on.
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkDataObject"
        )
        self._regions_as_blocks = False
        self._filenames = []
        self._file_format = None
        self._points_only = False
        self._piece = -1
        self._ghosts = "keep"
        self._deflection = 0.0
        self._angular_deflection = 0.5
        self._point_arrays = vtkDataArraySelection()
        self._cell_arrays = vtkDataArraySelection()
        for selection in (self._point_arrays, self._cell_arrays):
            selection.AddObserver("ModifiedEvent", _modified_callback(self))

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

    @smproperty.dataarrayselection(name="PointArrays")
    def GetPointArraySelection(self):
        return self._point_arrays

    @smproperty.dataarrayselection(name="CellArrays")
    def GetCellArraySelection(self):
        return self._cell_arrays

    @smproperty.intvector(name="PointsOnly", default_values="0")
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetPointsOnly(self, value):
        if self._points_only != bool(value):
            self._points_only = bool(value)
            self.Modified()

    @smproperty.intvector(name="Piece", default_values="-1")
    def SetPiece(self, value):
        # -1 reads every piece of a partitioned file (pvtu, pvd, VTKHDF, ...).
        if self._piece != int(value):
            self._piece = int(value)
            self.Modified()

    @smproperty.stringvector(name="GhostsInfo", information_only="1")
    def GetGhostPolicies(self):
        return list(_paraview.GHOST_POLICIES)

    @smproperty.stringvector(
        name="Ghosts", number_of_elements="1", default_values="keep"
    )
    @smdomain.xml(
        """
        <StringListDomain name="list">
            <RequiredProperties>
                <Property name="GhostsInfo" function="StringInfo"/>
            </RequiredProperties>
        </StringListDomain>
        """
    )
    def SetGhosts(self, value):
        if self._ghosts != value:
            self._ghosts = value
            self.Modified()

    @smproperty.doublevector(name="Deflection", default_values="0.0")
    def SetDeflection(self, value):
        # STEP and IGES: the chordal deflection, a length; 0 takes 1e-3 of
        # the bounding-box diagonal.
        if self._deflection != float(value):
            self._deflection = float(value)
            self.Modified()

    @smproperty.doublevector(name="AngularDeflection", default_values="0.5")
    def SetAngularDeflection(self, value):
        # STEP and IGES: the angular deflection, in radians.
        if self._angular_deflection != float(value):
            self._angular_deflection = float(value)
            self.Modified()

    @smproperty.intvector(name="RegionsAsBlocks", default_values="0")
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetRegionsAsBlocks(self, value):
        # One block per named region, in a data assembly, as ParaView's own
        # Exodus reader shows its blocks and sets.
        if self._regions_as_blocks != bool(value):
            self._regions_as_blocks = bool(value)
            self.Modified()

    @smproperty.stringvector(name="AddOns", information_only="1")
    def GetAddOns(self):
        # Which add-ons are loaded, e.g. "occt: no": a STEP read without
        # OpenCASCADE fails, and this says why.
        return [
            f"{name}: {'yes' if ok else 'no'}"
            for name, ok in _paraview.addon_status().items()
        ]

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

    def RequestDataObject(self, request, inInfoVec, outInfoVec):
        want = (
            vtkPartitionedDataSetCollection
            if self._regions_as_blocks
            else vtkUnstructuredGrid
        )
        info = outInfoVec.GetInformationObject(0)
        output = info.Get(vtkDataObject.DATA_OBJECT())
        if output is None or not output.IsA(want.__name__):
            info.Set(vtkDataObject.DATA_OBJECT(), want())
        return 1

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
        if self._filenames:
            meta = _paraview.metadata(self._filenames[0], self._file_format)
            for selection, key in (
                (self._point_arrays, "point_data_names"),
                (self._cell_arrays, "cell_data_names"),
            ):
                for name in meta.get(key, []):
                    if not selection.ArrayExists(name):
                        selection.AddArray(name)
        return 1

    def _read_kwargs(self):
        def ticked(selection):
            return [
                selection.GetArrayName(i)
                for i in range(selection.GetNumberOfArrays())
                if selection.ArrayIsEnabled(selection.GetArrayName(i))
            ]

        def names(selection):
            return [
                selection.GetArrayName(i) for i in range(selection.GetNumberOfArrays())
            ]

        return _paraview.read_kwargs(
            names(self._point_arrays),
            names(self._cell_arrays),
            ticked(self._point_arrays),
            ticked(self._cell_arrays),
            points_only=self._points_only,
            piece=self._piece,
            ghosts=self._ghosts,
        )

    def RequestData(self, request, inInfoVec, outInfoVec):
        executive = self.GetExecutive()
        info = outInfoVec.GetInformationObject(0)
        plan = self._plan()
        times = _paraview.time_values(plan)
        t = None
        if times and info.Has(executive.UPDATE_TIME_STEP()):
            t = info.Get(executive.UPDATE_TIME_STEP())
        index = _paraview.step_index(times, t)
        mesh = _paraview.read_step(
            plan,
            index,
            self._file_format,
            deflection=self._deflection,
            angular_deflection=self._angular_deflection,
            **self._read_kwargs(),
        )
        output = vtkDataObject.GetData(outInfoVec)
        if self._regions_as_blocks:
            _paraview.fill_collection(output, mesh)
        else:
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
@smdomain.datatype(
    dataTypes=["vtkDataSet", "vtkCompositeDataSet"], composite_data_supported=True
)
class MeshioPlusPlusWriter(VTKPythonAlgorithmBase):
    def __init__(self):
        # Any dataset (poly data, image data, ...) or a composite one, whose
        # leaves are merged with their block names as cell regions.
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=1, nOutputPorts=0, inputType="vtkDataObject"
        )
        self._file_format = None
        self._encoding = "default"
        self._codec = "none"
        self._float_format = ""
        self._provenance = True
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

    @smproperty.stringvector(name="WriteFormatInfo", information_only="1")
    def GetWriteFormats(self):
        return meshio_output_filetypes

    @smproperty.stringvector(name="FileFormat", number_of_elements="1")
    @smdomain.xml(
        """
        <StringListDomain name="list">
            <RequiredProperties>
                <Property name="WriteFormatInfo" function="StringInfo"/>
            </RequiredProperties>
        </StringListDomain>
        """
    )
    def SetFileFormat(self, file_format):
        file_format = None if file_format == "automatic" else file_format
        if self._file_format != file_format:
            self._file_format = file_format
            self.Modified()

    @smproperty.stringvector(name="EncodingInfo", information_only="1")
    def GetEncodings(self):
        return list(_paraview.ENCODINGS)

    @smproperty.stringvector(
        name="Encoding", number_of_elements="1", default_values="default"
    )
    @smdomain.xml(
        """
        <StringListDomain name="list">
            <RequiredProperties>
                <Property name="EncodingInfo" function="StringInfo"/>
            </RequiredProperties>
        </StringListDomain>
        """
    )
    def SetEncoding(self, value):
        if self._encoding != value:
            self._encoding = value
            self.Modified()

    @smproperty.stringvector(name="CodecInfo", information_only="1")
    def GetCodecs(self):
        return list(_paraview.CODECS)

    @smproperty.stringvector(
        name="Codec", number_of_elements="1", default_values="none"
    )
    @smdomain.xml(
        """
        <StringListDomain name="list">
            <RequiredProperties>
                <Property name="CodecInfo" function="StringInfo"/>
            </RequiredProperties>
        </StringListDomain>
        """
    )
    def SetCodec(self, value):
        if self._codec != value:
            self._codec = value
            self.Modified()

    @smproperty.stringvector(
        name="FloatFormat", number_of_elements="1", default_values=""
    )
    def SetFloatFormat(self, value):
        # A printf-style float format for text formats, e.g. ".6e".
        if self._float_format != (value or ""):
            self._float_format = value or ""
            self.Modified()

    @smproperty.intvector(name="Provenance", default_values="1")
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetProvenance(self, value):
        # Record ParaView as the source in the file's provenance block.
        if self._provenance != bool(value):
            self._provenance = bool(value)
            self.Modified()

    def _write_kwargs(self):
        return _paraview.write_kwargs(
            self._filename,
            self._file_format,
            self._encoding,
            self._codec,
            self._float_format,
        )

    def _provenance_scope(self):
        return _paraview.provenance_scope(
            f"ParaView {_paraview_version()}" if self._provenance else None,
            self._filename,
            self._file_format,
            self._encoding,
            self._codec,
            self._float_format,
        )

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
        mesh = _paraview.mesh_from_data_object(vtkDataObject.GetData(inInfoVec[0]))
        kwargs = self._write_kwargs()
        if len(self._times) < 2:
            with self._provenance_scope():
                meshioplusplus.write(self._filename, mesh, **kwargs)
            return 1
        # Every step, one pass of the pipeline each (CONTINUE_EXECUTING): a
        # one-file-per-step format is written as each step arrives; a format
        # holding a series in one file is written once the last one has.
        target = _paraview.series_path(self._filename, self._file_format)
        t, n = self._times[self._index], len(self._times)
        if target == str(self._filename):
            self._steps.append((t, mesh))
        else:
            with self._provenance_scope():
                _paraview.write_one_step(target, self._index, n, t, mesh, **kwargs)
        executive = self.GetExecutive()
        if self._index + 1 < n:
            self._index += 1
            request.Set(executive.CONTINUE_EXECUTING(), 1)
            return 1
        request.Remove(executive.CONTINUE_EXECUTING())
        if self._steps:
            with self._provenance_scope():
                meshioplusplus.write_sequence(self._filename, self._steps, **kwargs)
        self._index, self._steps = 0, []
        return 1

    def Write(self):
        self.Modified()
        self.Update()


@smproxy.filter(
    name="MeshioPlusPlusPipeline", label="meshio++ Pipeline", support_reload=False
)
@smhint.xml('<ShowInMenu category="meshio++"/>')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(
    dataTypes=["vtkDataSet", "vtkCompositeDataSet"], composite_data_supported=True
)
class MeshioPlusPlusPipeline(VTKPythonAlgorithmBase):
    """The settings pipeline's Operations, applied to the input
    (apply_operations): every operation the pipeline has, in one filter."""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=1,
            nOutputPorts=1,
            inputType="vtkDataObject",
            outputType="vtkUnstructuredGrid",
        )
        self._operations = "[]"

    @smproperty.stringvector(name="Operations", default_values="[]")
    @smhint.xml('<Widget type="multi_line" syntax="json"/>')
    def SetOperations(self, text):
        # The settings pipeline's Operations list as JSON, e.g.
        # [{"Op": "Refine"}, {"Op": "Smooth", "Iterations": 5}].
        if self._operations != text:
            self._operations = text
            self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        mesh = _paraview.mesh_from_data_object(
            vtkDataObject.GetData(inInfoVec[0]), _paraview.FILTER_OP
        )
        mesh, report = _paraview.run_operations(mesh, self._operations)
        _paraview.fill_grid_with_report(
            vtkUnstructuredGrid.GetData(outInfoVec), mesh, report
        )
        return 1


# --------------------------------------------------------------------------- #
# Dedicated filters (roadmap 6.9): thin proxies with typed properties.        #
# --------------------------------------------------------------------------- #
_DATASETS = ["vtkDataSet", "vtkCompositeDataSet"]


def _setter(key, kind):
    def setter(self, value):
        if self._values.get(key) != value:
            self._values[key] = value
            self.Modified()

    setter.__name__ = "Set" + key
    setter.__qualname__ = setter.__name__
    return setter


def _typed_property(key, kind, default):
    """The decorators of one step key's property."""
    if kind == "flag":
        return [
            smproperty.intvector(name=key, default_values=str(int(bool(default)))),
            smdomain.xml('<BooleanDomain name="bool"/>'),
        ]
    if kind == "number":
        return [smproperty.doublevector(name=key, default_values=str(float(default)))]
    # ParaView writes the attribute as given, so a default such as "<" (a
    # comparison) is escaped here.
    return [
        smproperty.stringvector(
            name=key,
            number_of_elements="1",
            default_values=html.escape(str(default), quote=True),
        )
    ]


def _step_filter(op, label):
    """A filter running the one pipeline step `op` on its input."""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=1,
            nOutputPorts=1,
            inputType="vtkDataObject",
            outputType="vtkUnstructuredGrid",
        )
        self._values = {}
        self._extra = ""

    def SetExtraKeys(self, text):
        # Any other key of the step, as a JSON object.
        if self._extra != text:
            self._extra = text
            self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        mesh = _paraview.mesh_from_data_object(
            vtkDataObject.GetData(inInfoVec[0]), _paraview.FILTER_OP
        )
        step = _paraview.step_from_properties(op, self._values, self._extra)
        mesh, report = _paraview.run_operations(mesh, [step])
        _paraview.fill_grid_with_report(
            vtkUnstructuredGrid.GetData(outInfoVec), mesh, report
        )
        return 1

    namespace = {"__init__": __init__, "RequestData": RequestData}
    for key, kind, default in _paraview.step_keys(op):
        method = _setter(key, kind)
        for decorate in reversed(_typed_property(key, kind, default)):
            method = decorate(method)
        namespace[method.__name__] = method
    namespace["SetExtraKeys"] = smproperty.stringvector(
        name="ExtraKeys", number_of_elements="1", default_values=""
    )(SetExtraKeys)
    cls = type("MeshioPlusPlus" + op, (VTKPythonAlgorithmBase,), namespace)
    cls.__doc__ = f"The settings pipeline's {op} step (doc/pipeline.md)."
    # Applied bottom-up, as the decorator stacks below are: smproxy.filter
    # renders the proxy's XML, so it comes last.
    decorated = smdomain.datatype(dataTypes=_DATASETS, composite_data_supported=True)(
        cls
    )
    decorated = smproperty.input(name="Input", port_index=0)(decorated)
    decorated = smhint.xml('<ShowInMenu category="meshio++"/>')(decorated)
    decorated = smproxy.filter(
        name=cls.__name__, label=f"meshio++ {label}", support_reload=False
    )(decorated)
    globals()[cls.__name__] = decorated


for _op, _label in _paraview.STEP_FILTERS:
    _step_filter(_op, _label)


def _mesh_in(inInfoVec, port):
    return _paraview.mesh_from_data_object(
        vtkDataObject.GetData(inInfoVec[port]), _paraview.FILTER_OP
    )


@smproxy.filter(
    name="MeshioPlusPlusPartition", label="meshio++ Partition", support_reload=False
)
@smhint.xml('<ShowInMenu category="meshio++"/>')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=_DATASETS, composite_data_supported=True)
class MeshioPlusPlusPartition(VTKPythonAlgorithmBase):
    """partition(): one partition per piece, with its ghost layers."""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=1,
            nOutputPorts=1,
            inputType="vtkDataObject",
            outputType="vtkPartitionedDataSet",
        )
        self._nparts, self._ghosts, self._method = 2, 0, "auto"

    @smproperty.intvector(name="Nparts", default_values="2")
    def SetNparts(self, value):
        if self._nparts != int(value):
            self._nparts = int(value)
            self.Modified()

    @smproperty.intvector(name="GhostLayers", default_values="0")
    def SetGhostLayers(self, value):
        if self._ghosts != int(value):
            self._ghosts = int(value)
            self.Modified()

    @smproperty.stringvector(
        name="Method", number_of_elements="1", default_values="auto"
    )
    def SetMethod(self, value):
        if self._method != value:
            self._method = value
            self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        pieces = _paraview.partition_pieces(
            _mesh_in(inInfoVec, 0), self._nparts, self._ghosts, method=self._method
        )
        _paraview.fill_partitioned(vtkPartitionedDataSet.GetData(outInfoVec), pieces)
        return 1


class _TwoInputs(VTKPythonAlgorithmBase):
    """A filter of two meshes, the second on port 1."""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=2,
            nOutputPorts=1,
            inputType="vtkDataObject",
            outputType="vtkUnstructuredGrid",
        )

    def _run(self, first, second):
        raise NotImplementedError

    def RequestData(self, request, inInfoVec, outInfoVec):
        mesh, report = self._run(_mesh_in(inInfoVec, 0), _mesh_in(inInfoVec, 1))
        _paraview.fill_grid_with_report(
            vtkUnstructuredGrid.GetData(outInfoVec), mesh, report
        )
        return 1


@smproxy.filter(
    name="MeshioPlusPlusConservativeInterpolate",
    label="meshio++ Conservative Interpolate",
    support_reload=False,
)
@smhint.xml('<ShowInMenu category="meshio++"/>')
@smproperty.input(name="Source", port_index=1)
@smdomain.datatype(dataTypes=_DATASETS, composite_data_supported=True)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=_DATASETS, composite_data_supported=True)
class MeshioPlusPlusConservativeInterpolate(_TwoInputs):
    """conservative_interpolate(): Source's data onto Input, mass-preserving."""

    def __init__(self):
        _TwoInputs.__init__(self)
        self._on_conflict = "error"

    @smproperty.stringvector(
        name="OnConflict", number_of_elements="1", default_values="error"
    )
    def SetOnConflict(self, value):
        # An array the target already has: error, overwrite or suffix.
        if self._on_conflict != value:
            self._on_conflict = value
            self.Modified()

    def _run(self, target, source):
        out = meshioplusplus.conservative_interpolate(
            source, target, on_conflict=self._on_conflict
        )
        return out, {}


@smproxy.filter(
    name="MeshioPlusPlusShrinkwrap", label="meshio++ Shrinkwrap", support_reload=False
)
@smhint.xml('<ShowInMenu category="meshio++"/>')
@smproperty.input(name="Target", port_index=1)
@smdomain.datatype(dataTypes=_DATASETS, composite_data_supported=True)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=_DATASETS, composite_data_supported=True)
class MeshioPlusPlusShrinkwrap(_TwoInputs):
    """shrinkwrap(): Input's points moved onto the Target surface."""

    def __init__(self):
        _TwoInputs.__init__(self)
        self._offset = 0.0

    @smproperty.doublevector(name="Offset", default_values="0.0")
    def SetOffset(self, value):
        if self._offset != float(value):
            self._offset = float(value)
            self.Modified()

    def _run(self, mesh, target):
        return meshioplusplus.shrinkwrap(
            mesh, target, offset=self._offset, return_report=True
        )


@smproxy.filter(
    name="MeshioPlusPlusUndoGreen", label="meshio++ Undo Green", support_reload=False
)
@smhint.xml('<ShowInMenu category="meshio++"/>')
@smproperty.input(name="Coarse", port_index=1)
@smdomain.datatype(dataTypes=_DATASETS, composite_data_supported=True)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=_DATASETS, composite_data_supported=True)
class MeshioPlusPlusUndoGreen(_TwoInputs):
    """undo_green(): Input (a refine output) with its green cells restored
    from Coarse, the mesh it was refined from."""

    def _run(self, fine, coarse):
        return meshioplusplus.undo_green(coarse, fine, return_report=True)


@smproxy.filter(
    name="MeshioPlusPlusPeriodicNodes",
    label="meshio++ Match Periodic Nodes",
    support_reload=False,
)
@smhint.xml('<ShowInMenu category="meshio++"/>')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=_DATASETS, composite_data_supported=True)
class MeshioPlusPlusPeriodicNodes(VTKPythonAlgorithmBase):
    """match_periodic_nodes(): a periodic:master point array."""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=1,
            nOutputPorts=1,
            inputType="vtkDataObject",
            outputType="vtkUnstructuredGrid",
        )
        self._slave, self._master = "", ""
        self._translate = (0.0, 0.0, 0.0)
        self._atol = 1e-8

    @smproperty.stringvector(name="Slave", number_of_elements="1", default_values="")
    def SetSlave(self, value):
        if self._slave != value:
            self._slave = value
            self.Modified()

    @smproperty.stringvector(name="Master", number_of_elements="1", default_values="")
    def SetMaster(self, value):
        if self._master != value:
            self._master = value
            self.Modified()

    @smproperty.doublevector(name="Translate", default_values=[0.0, 0.0, 0.0])
    def SetTranslate(self, x, y, z):
        if self._translate != (x, y, z):
            self._translate = (x, y, z)
            self.Modified()

    @smproperty.doublevector(name="Atol", default_values="1e-08")
    def SetAtol(self, value):
        if self._atol != float(value):
            self._atol = float(value)
            self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        mesh = _paraview.with_periodic_master(
            _mesh_in(inInfoVec, 0),
            self._slave,
            self._master,
            self._translate,
            self._atol,
        )
        _paraview.fill_grid(vtkUnstructuredGrid.GetData(outInfoVec), mesh)
        return 1


@smproxy.filter(
    name="MeshioPlusPlusContactPairs",
    label="meshio++ Contact Pairs",
    support_reload=False,
)
@smhint.xml('<ShowInMenu category="meshio++"/>')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=_DATASETS, composite_data_supported=True)
class MeshioPlusPlusContactPairs(VTKPythonAlgorithmBase):
    """contact_pairs(): the slave points, with their gap, normal and master."""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=1,
            nOutputPorts=1,
            inputType="vtkDataObject",
            outputType="vtkUnstructuredGrid",
        )
        self._slave, self._master, self._tolerance = "", "", 0.0

    @smproperty.stringvector(
        name="SlavePoints", number_of_elements="1", default_values=""
    )
    def SetSlavePoints(self, value):
        if self._slave != value:
            self._slave = value
            self.Modified()

    @smproperty.stringvector(
        name="MasterCells", number_of_elements="1", default_values=""
    )
    def SetMasterCells(self, value):
        if self._master != value:
            self._master = value
            self.Modified()

    @smproperty.doublevector(name="Tolerance", default_values="0.0")
    def SetTolerance(self, value):
        if self._tolerance != float(value):
            self._tolerance = float(value)
            self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        points = _paraview.contact_points(
            _mesh_in(inInfoVec, 0), self._slave, self._master, self._tolerance
        )
        _paraview.fill_grid(vtkUnstructuredGrid.GetData(outInfoVec), points)
        return 1


_alias(MeshioPlusPlusReader, "meshio++ reader", smproxy.source)
_alias(
    MeshioPlusPlusWriter,
    "meshio++ Writer",
    lambda **kw: smproxy.source(group="writers", type="WriterProxy", **kw),
)
