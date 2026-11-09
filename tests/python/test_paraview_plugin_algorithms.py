"""The ParaView plugin's algorithms, run on the PyPI ``vtk`` wheel.

The plugin's proxies are ``VTKPythonAlgorithmBase`` subclasses, which VTK
alone can run; only their ``smproxy``/``smproperty`` decorators come from
ParaView. Here those are stand-ins that leave each class and method as it is,
so the reader, the writer and the filters run their ``RequestData`` through a
real VTK pipeline in the default test matrix. ``test_paraview_plugin.py``
loads the same file under a real ``pvpython``.
"""

import importlib.util
import json
import pathlib
import sys
import types

import numpy as np
import pytest

import meshioplusplus

vtk = pytest.importorskip("vtk")
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase  # noqa: E402

PLUGIN = (
    pathlib.Path(__file__).resolve().parents[2]
    / "tools"
    / ("paraview-meshioplusplus-plugin.py")
)


class _Decorators:
    """Every ``smproxy.*``/``smproperty.*``/... factory: a decorator that
    returns its argument, marked as ParaView's are (``_pv_original_func``)."""

    def __getattr__(self, name):
        def factory(*args, **kwargs):
            def decorate(obj):
                try:
                    obj._pv_original_func = obj
                except (AttributeError, TypeError):
                    pass
                return obj

            return decorate

        return factory


@pytest.fixture(scope="module")
def plugin():
    stubs = {
        "paraview": types.ModuleType("paraview"),
        "paraview.util": types.ModuleType("paraview.util"),
        "paraview.util.vtkAlgorithm": types.ModuleType("paraview.util.vtkAlgorithm"),
    }
    algorithm = stubs["paraview.util.vtkAlgorithm"]
    algorithm.VTKPythonAlgorithmBase = VTKPythonAlgorithmBase
    for name in ("smproxy", "smproperty", "smdomain", "smhint"):
        setattr(algorithm, name, _Decorators())
    saved = {k: sys.modules.get(k) for k in stubs}
    sys.modules.update(stubs)
    try:
        spec = importlib.util.spec_from_file_location("mio_pv_plugin", PLUGIN)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        yield module
    finally:
        for key, value in saved.items():
            if value is None:
                sys.modules.pop(key, None)
            else:
                sys.modules[key] = value


def _tets():
    return meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]),
        [("tetra", [[0, 1, 2, 3], [1, 2, 3, 4]])],
        point_data={"u": np.arange(5.0), "v": np.ones(5)},
    )


def _reader(plugin, path, **properties):
    reader = plugin.MeshioPlusPlusReader()
    reader.AddFileName(str(path))
    for name, value in properties.items():
        getattr(reader, "Set" + name)(value)
    return reader


def test_reader(plugin, tmp_path):
    path = tmp_path / "m.vtu"
    meshioplusplus.write(path, _tets())
    reader = _reader(plugin, path)
    reader.Update()
    grid = reader.GetOutputDataObject(0)
    assert grid.GetNumberOfCells() == 2
    assert grid.GetPointData().GetArray("u") is not None


def test_reader_array_selection(plugin, tmp_path):
    path = tmp_path / "m.vtu"
    meshioplusplus.write(path, _tets())
    reader = _reader(plugin, path)
    reader.UpdateInformation()
    selection = reader.GetPointArraySelection()
    assert selection.GetNumberOfArrays() == 2
    selection.DisableArray("u")
    reader.Update()
    data = reader.GetOutputDataObject(0).GetPointData()
    assert data.GetArray("u") is None and data.GetArray("v") is not None


def test_reader_time_steps_and_writer(plugin, tmp_path):
    steps = []
    for t in (0.0, 0.5, 2.0):
        mesh = _tets()
        mesh.point_data["u"] = np.full(5, 10.0 * t)
        steps.append((t, mesh))
    src = tmp_path / "in.frd"
    meshioplusplus.write_sequence(src, steps)
    reader = _reader(plugin, src)
    reader.UpdateInformation()
    info = reader.GetOutputInformation(0)
    executive = reader.GetExecutive()
    assert list(info.Get(executive.TIME_STEPS())) == [0.0, 0.5, 2.0]
    reader.UpdateTimeStep(2.0)
    u = reader.GetOutputDataObject(0).GetPointData().GetArray("u")
    assert u.GetValue(0) == 20.0

    dst = tmp_path / "out.frd"
    writer = plugin.MeshioPlusPlusWriter()
    writer.SetInputConnection(reader.GetOutputPort())
    writer.SetFileName(str(dst))
    writer.SetWriteTimeSteps(1)
    writer.Write()
    assert meshioplusplus.read_metadata(dst)["time_values"] == [0.0, 0.5, 2.0]


def test_regions_as_blocks(plugin, tmp_path):
    from meshioplusplus._regions import Region

    mesh = _tets()
    mesh.regions = [Region("a", "cell", [0], 3, 1), Region("b", "cell", [1], 3, 2)]
    path = tmp_path / "r.vtu"
    meshioplusplus.write(path, mesh)
    reader = _reader(plugin, path, RegionsAsBlocks=1)
    reader.Update()
    data = reader.GetOutputDataObject(0)
    assert data.IsA("vtkPartitionedDataSetCollection")
    assert data.GetNumberOfPartitionedDataSets() == 2


def test_writer_options_and_provenance(plugin, tmp_path):
    sphere = vtk.vtkSphereSource()
    dst = tmp_path / "s.vtu"
    writer = plugin.MeshioPlusPlusWriter()
    writer.SetInputConnection(sphere.GetOutputPort())
    writer.SetFileName(str(dst))
    writer.SetEncoding("ascii")
    writer.Write()
    assert b'format="ascii"' in dst.read_bytes()
    assert any(
        "ParaView" in line for line in meshioplusplus.read_metadata(dst)["provenance"]
    )


def test_pipeline_filter(plugin):
    source = vtk.vtkCellTypeSource()
    source.SetCellType(vtk.VTK_TETRA)
    source.Update()
    cells = source.GetOutput().GetNumberOfCells()
    pipeline = plugin.MeshioPlusPlusPipeline()
    pipeline.SetInputConnection(source.GetOutputPort())
    pipeline.SetOperations('[{"Op": "Refine"}]')
    pipeline.Update()
    out = pipeline.GetOutputDataObject(0)
    assert out.GetNumberOfCells() == 8 * cells
    report = out.GetFieldData().GetAbstractArray("meshioplusplus:report").GetValue(0)
    assert json.loads(report)["steps"][0]["op"] == "Refine"


def test_dedicated_step_filters_exist(plugin):
    from meshioplusplus import _paraview

    for op, _ in _paraview.STEP_FILTERS:
        cls = getattr(plugin, "MeshioPlusPlus" + op)
        for key, _, _ in _paraview.step_keys(op):
            assert hasattr(cls, "Set" + key), (op, key)


def test_remesh_filter(plugin):
    """The roadmap's probe: remesh with 500 clusters on a sphere gives about
    500 points."""
    sphere = vtk.vtkSphereSource()
    sphere.SetThetaResolution(64)
    sphere.SetPhiResolution(64)
    remesh = plugin.MeshioPlusPlusRemesh()
    remesh.SetInputConnection(sphere.GetOutputPort())
    remesh.SetNumClusters(500.0)
    remesh.Update()
    points = remesh.GetOutputDataObject(0).GetNumberOfPoints()
    assert 400 <= points <= 600


def test_refine_filter_with_extra_keys(plugin):
    source = vtk.vtkCellTypeSource()
    source.SetCellType(vtk.VTK_TETRA)
    source.Update()
    refine = plugin.MeshioPlusPlusRefine()
    refine.SetInputConnection(source.GetOutputPort())
    refine.SetLevels(2.0)
    refine.SetExtraKeys('{"RecordLevels": true}')
    refine.Update()
    out = refine.GetOutputDataObject(0)
    assert out.GetNumberOfCells() == 64 * source.GetOutput().GetNumberOfCells()
    assert out.GetCellData().GetArray("refine:level") is not None


def test_partition_filter_with_ghosts(plugin):
    source = vtk.vtkCellTypeSource()
    source.SetCellType(vtk.VTK_HEXAHEDRON)
    source.SetBlocksDimensions(4, 4, 4)
    partition = plugin.MeshioPlusPlusPartition()
    partition.SetInputConnection(source.GetOutputPort())
    partition.SetNparts(2)
    partition.SetGhostLayers(1)
    partition.Update()
    out = partition.GetOutputDataObject(0)
    assert out.GetNumberOfPartitions() == 2
    total = sum(out.GetPartition(i).GetNumberOfCells() for i in range(2))
    assert total > 64  # the ghost layers are extra cells


def test_two_input_filter(plugin):
    """conservative_interpolate: the source's cell data on the target."""
    from meshioplusplus import _paraview

    source = _tets()
    source.cell_data["c"] = [np.array([1.0, 1.0])]
    target = _tets()
    f = plugin.MeshioPlusPlusConservativeInterpolate()
    f.SetOnConflict("overwrite")  # the target has the source's point arrays
    f.SetInputDataObject(0, _paraview.grid_from_mesh(target))
    f.SetInputDataObject(1, _paraview.grid_from_mesh(source))
    f.Update()
    out = f.GetOutputDataObject(0)
    np.testing.assert_allclose(
        [out.GetCellData().GetArray("c").GetValue(i) for i in range(2)], [1.0, 1.0]
    )
