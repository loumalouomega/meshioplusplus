"""VTKHDF against VTK's own ``vtkHDFReader`` / ``vtkHDFWriter``.

The round-trip tests in ``test_vtkhdf.py`` cannot catch a file that *we* read and
Kitware's reader does not -- a variable-length ``Type``, an untracked composite, a
missing ``Index``. These do, and they cover the roadmap's "done when": a transient
partitioned case opens with a working time slider and per-piece blocks, and a
VTK-written sample round-trips bit-exact on fields.

Skipped where ``vtk`` (pulled in by the ``pyvista`` extra) is not installed.
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus._regions import Region

h5py = pytest.importorskip("h5py")
pytest.importorskip("vtkmodules.vtkIOHDF")

import vtk  # noqa: E402
from vtkmodules.util import numpy_support as ns  # noqa: E402
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase  # noqa: E402
from vtkmodules.vtkCommonDataModel import (  # noqa: E402
    vtkCompositeDataSet,
    vtkPartitionedDataSet,
    vtkPolyData,
    vtkUnstructuredGrid,
)
from vtkmodules.vtkCommonExecutionModel import (  # noqa: E402
    vtkStreamingDemandDrivenPipeline as SDDP,
)
from vtkmodules.vtkIOHDF import vtkHDFReader, vtkHDFWriter  # noqa: E402

CUBE = [
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1],
]
CUBE_FACES = [
    [0, 3, 2, 1],
    [4, 5, 6, 7],
    [0, 1, 5, 4],
    [1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 0, 4, 7],
]


def _vtk_read(path):
    reader = vtkHDFReader()
    assert reader.CanReadFile(str(path)) == 1, "vtkHDFReader refuses the file"
    reader.SetFileName(str(path))
    reader.UpdateInformation()
    return reader


def _arr(a):
    return ns.vtk_to_numpy(a)


def _numpy_array(name, values):
    a = ns.numpy_to_vtk(np.ascontiguousarray(values), deep=1)
    a.SetName(name)
    return a


def _mesh(u=0.0):
    rng = np.random.default_rng(3)
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 1],
            [2, 0, 0],
            [2, 1, 0],
            [2, 0, 1],
            [3, 0, 0],
            [3, 1, 0],
            [3, 0, 1],
            [3, 1, 1],
        ],
        dtype=float,
    )
    return meshioplusplus.Mesh(
        pts,
        [("tetra", [[0, 1, 2, 3], [1, 2, 3, 4]]), ("wedge", [[5, 6, 8, 7, 9, 10]])],
        point_data={
            "u": rng.random(12) + u,
            "vel": rng.random((12, 3)),
            "id": np.arange(12, dtype=np.int32),
        },
        cell_data={
            "p": [rng.random(2), rng.random(1)],
            "tens": [rng.random((2, 9)), rng.random((1, 9))],
        },
        field_data={"gain": np.array([1.5, 2.5, 3.5])},
    )


@pytest.mark.parametrize("compression", [None, "gzip"])
def test_vtk_reads_an_unstructured_grid_bit_exact(compression, tmp_path):
    mesh = _mesh()
    p = tmp_path / "ug.vtkhdf"
    meshioplusplus.vtkhdf.write(p, mesh, compression=compression)
    reader = _vtk_read(p)
    reader.Update()
    g = reader.GetOutputDataObject(0)
    assert g.IsA("vtkUnstructuredGrid")
    assert (g.GetNumberOfPoints(), g.GetNumberOfCells()) == (12, 3)
    assert np.array_equal(_arr(g.GetPoints().GetData()), mesh.points)
    assert [g.GetCellType(i) for i in range(3)] == [
        vtk.VTK_TETRA,
        vtk.VTK_TETRA,
        vtk.VTK_WEDGE,
    ]
    ids = [g.GetCell(2).GetPointId(j) for j in range(6)]
    assert ids == [5, 8, 6, 7, 10, 9]  # meshio++'s wedge order permuted into VTK's
    pd, cd = g.GetPointData(), g.GetCellData()
    for name in ("u", "vel", "id"):
        assert np.array_equal(_arr(pd.GetArray(name)), mesh.point_data[name]), name
    assert _arr(pd.GetArray("id")).dtype == np.int32
    assert np.array_equal(_arr(cd.GetArray("p")), np.concatenate(mesh.cell_data["p"]))
    assert np.array_equal(
        _arr(cd.GetArray("tens")), np.concatenate(mesh.cell_data["tens"])
    )
    assert np.array_equal(
        _arr(g.GetFieldData().GetArray("gain")), mesh.field_data["gain"]
    )


def test_vtk_reads_a_version_1_0_file(tmp_path):
    p = tmp_path / "v10.vtkhdf"
    meshioplusplus.vtkhdf.write(p, _mesh(), version=(1, 0))
    reader = _vtk_read(p)
    reader.Update()
    assert reader.GetOutputDataObject(0).GetNumberOfCells() == 3


def test_vtk_reads_mixed_polyhedra(tmp_path):
    pts = np.array(CUBE + [[0, 0, 2]], dtype=float)
    mesh = meshioplusplus.Mesh(
        pts,
        [
            ("tetra", [[4, 5, 6, 8]]),
            ("polyhedron8", [[np.array(f) for f in CUBE_FACES]]),
            ("tetra", [[0, 1, 3, 4]]),
        ],
        point_data={"u": np.arange(9.0)},
        cell_data={"p": [np.array([1.0]), np.array([2.0]), np.array([3.0])]},
    )
    p = tmp_path / "poly.vtkhdf"
    meshioplusplus.vtkhdf.write(p, mesh)
    reader = _vtk_read(p)
    reader.Update()
    g = reader.GetOutputDataObject(0)
    assert [g.GetCellType(i) for i in range(3)] == [
        vtk.VTK_TETRA,
        vtk.VTK_POLYHEDRON,
        vtk.VTK_TETRA,
    ]
    assert g.GetCell(1).GetNumberOfFaces() == 6
    assert _arr(g.GetCellData().GetArray("p")).tolist() == [1.0, 2.0, 3.0]


def test_vtk_reads_polydata(tmp_path):
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [2, 2, 2], [3, 3, 3]], dtype=float
    )
    mesh = meshioplusplus.Mesh(
        pts,
        [
            ("triangle", [[0, 1, 2]]),
            ("vertex", [[4]]),
            ("quad", [[0, 1, 3, 2]]),
            ("line", [[0, 4]]),
            ("polygon", np.array([[0, 1, 3, 2, 5]])),
        ],
        cell_data={
            "p": [
                np.array([10.0]),
                np.array([20.0]),
                np.array([30.0]),
                np.array([40.0]),
                np.array([50.0]),
            ]
        },
    )
    p = tmp_path / "pd.vtkhdf"
    meshioplusplus.vtkhdf.write(p, mesh, dataset_type="PolyData")
    reader = _vtk_read(p)
    reader.Update()
    g = reader.GetOutputDataObject(0)
    assert g.IsA("vtkPolyData")
    assert (
        g.GetNumberOfVerts(),
        g.GetNumberOfLines(),
        g.GetNumberOfPolys(),
        g.GetNumberOfStrips(),
    ) == (1, 1, 3, 0)
    assert _arr(g.GetCellData().GetArray("p")).tolist() == [
        20.0,
        40.0,
        10.0,
        30.0,
        50.0,
    ]


def _composite_mesh(with_field_data):
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 1],
            [5, 0, 0],
            [6, 0, 0],
            [5, 1, 0],
            [5, 0, 1],
        ],
        dtype=float,
    )
    mesh = meshioplusplus.Mesh(
        pts,
        [("tetra", [[0, 1, 2, 3], [1, 2, 3, 4], [5, 6, 7, 8]])],
        point_data={"u": np.arange(9.0)},
        cell_data={"p": [np.array([1.0, 2.0, 3.0])]},
        field_data={"g": np.array([9.0])} if with_field_data else {},
    )
    mesh.regions = [
        Region("zeta", "cell", np.array([2]), tag=0),
        Region("alpha", "cell", np.array([0, 1]), tag=1),
    ]
    return mesh


@pytest.mark.parametrize("with_field_data", [False, True])
def test_vtk_reads_a_partitioned_dataset_collection(with_field_data, tmp_path):
    p = tmp_path / "pdc.vtkhdf"
    meshioplusplus.vtkhdf.write(
        p, _composite_mesh(with_field_data), dataset_type="PartitionedDataSetCollection"
    )
    reader = _vtk_read(p)
    reader.Update()
    g = reader.GetOutputDataObject(0)
    assert g.IsA("vtkPartitionedDataSetCollection")
    names = [
        g.GetMetaData(i).Get(vtkCompositeDataSet.NAME())
        for i in range(g.GetNumberOfPartitionedDataSets())
    ]
    assert names == [
        "zeta",
        "alpha",
    ], "block order must follow the (creation-ordered) Assembly, not the alphabet"
    zeta = g.GetPartitionedDataSet(0).GetPartition(0)
    alpha = g.GetPartitionedDataSet(1).GetPartition(0)
    assert (zeta.GetNumberOfCells(), alpha.GetNumberOfCells()) == (1, 2)
    assert _arr(zeta.GetCellData().GetArray("p")).tolist() == [3.0]
    assert _arr(alpha.GetCellData().GetArray("p")).tolist() == [1.0, 2.0]
    assert (
        zeta.GetNumberOfPoints() == 4 and zeta.GetPoint(0)[0] == 5.0
    )  # pruned to its own points


def test_vtk_reads_a_multiblock_dataset(tmp_path):
    p = tmp_path / "mb.vtkhdf"
    meshioplusplus.vtkhdf.write(
        p, _composite_mesh(False), dataset_type="MultiBlockDataSet"
    )
    reader = _vtk_read(p)
    reader.Update()
    g = reader.GetOutputDataObject(0)
    assert g.IsA("vtkMultiBlockDataSet")
    assert [
        g.GetMetaData(i).Get(vtkCompositeDataSet.NAME())
        for i in range(g.GetNumberOfBlocks())
    ] == ["zeta", "alpha"]


# ---- time: the on-disk form of the sequence engine ------------------------------ #
def _write_static_once(path, n_steps=3):
    """Geometry written once, fields appended per step: every step's offsets point at it."""
    from meshioplusplus.vtkhdf import _vtkhdf

    n, m = 5, 2
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float),
        [("tetra", [[0, 1, 2, 3], [1, 2, 3, 4]])],
    )
    meshioplusplus.vtkhdf.write(path, mesh)
    with h5py.File(path, "r+") as f:
        g = f["VTKHDF"]
        del g["PointData"], g["CellData"]
        g.create_group("PointData").create_dataset(
            "u", data=np.repeat(np.arange(n_steps, dtype=float), n)
        )
        g.create_group("CellData").create_dataset(
            "p", data=np.repeat(np.arange(n_steps) * 10.0, m)
        )
        s = g.create_group("Steps")
        s.attrs["NSteps"] = n_steps
        s["Values"] = np.arange(n_steps) * 0.5
        zeros = np.zeros(n_steps, dtype=np.int64)
        for name in (
            "PartOffsets",
            "PointOffsets",
            "CellOffsets",
            "ConnectivityIdOffsets",
        ):
            s[name] = zeros
        s["NumberOfParts"] = np.ones(n_steps, dtype=np.int64)
        s.create_group("PointDataOffsets")["u"] = np.arange(n_steps) * n
        s.create_group("CellDataOffsets")["p"] = np.arange(n_steps) * m
        s.create_group("FieldDataOffsets")
        s.create_group("FieldDataSizes")
    assert _vtkhdf.ROOT == "VTKHDF"


def test_vtk_time_slider_changes_fields_not_geometry(tmp_path):
    p = tmp_path / "once.vtkhdf"
    _write_static_once(p)
    reader = _vtk_read(p)
    info = reader.GetOutputInformation(0)
    assert list(info.Get(SDDP.TIME_STEPS())) == [0.0, 0.5, 1.0]
    seen = []
    for k, t in enumerate((0.0, 0.5, 1.0)):
        reader.UpdateTimeStep(t)
        g = reader.GetOutputDataObject(0)
        seen.append(
            (
                _arr(g.GetPoints().GetData()).copy(),
                _arr(g.GetPointData().GetArray("u")).tolist(),
                _arr(g.GetCellData().GetArray("p")).tolist(),
            )
        )
        assert seen[-1][1] == [float(k)] * 5 and seen[-1][2] == [10.0 * k] * 2
    assert all(
        np.array_equal(seen[0][0], s[0]) for s in seen
    )  # geometry identical at every step
    # and our own reader agrees step for step
    for k in range(3):
        m = meshioplusplus.read(p, time_step=k)
        assert m.point_data["u"].tolist() == [float(k)] * 5


# ---- VTK-written samples read by meshio++ --------------------------------------- #
def _tets(shift=0.0):
    g = vtkUnstructuredGrid()
    pts = vtk.vtkPoints()
    for q in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]:
        pts.InsertNextPoint(q[0] + shift, q[1], q[2])
    g.SetPoints(pts)
    g.InsertNextCell(vtk.VTK_TETRA, 4, [0, 1, 2, 3])
    g.InsertNextCell(vtk.VTK_TETRA, 4, [1, 2, 3, 4])
    return g


class _Temporal(VTKPythonAlgorithmBase):
    """A source producing ``builder(output, t)`` for t in 0, 0.5, 1."""

    def __init__(self, output_type, builder):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType=output_type
        )
        self.builder = builder
        self.times = [0.0, 0.5, 1.0]

    def RequestInformation(self, request, in_info, out_info):
        oi = out_info.GetInformationObject(0)
        oi.Set(SDDP.TIME_STEPS(), self.times, len(self.times))
        oi.Set(SDDP.TIME_RANGE(), [self.times[0], self.times[-1]], 2)
        return 1

    def RequestData(self, request, in_info, out_info):
        oi = out_info.GetInformationObject(0)
        t = oi.Get(SDDP.UPDATE_TIME_STEP()) if oi.Has(SDDP.UPDATE_TIME_STEP()) else 0.0
        self.builder(vtk.vtkDataObject.GetData(out_info), t)
        return 1


def _vtk_write(path, data=None, source=None):
    w = vtkHDFWriter()
    w.SetFileName(str(path))
    if source is not None:
        w.SetInputConnection(source.GetOutputPort())
        w.SetWriteAllTimeSteps(True)
    else:
        w.SetInputData(data)
    w.Write()


def test_vtk_written_static_sample_round_trips_bit_exact(tmp_path):
    rng = np.random.default_rng(5)
    g = _tets()
    u, vel, p = rng.random(5), rng.random((5, 3)), rng.random(2)
    g.GetPointData().AddArray(_numpy_array("u", u))
    g.GetPointData().AddArray(_numpy_array("vel", vel))
    g.GetCellData().AddArray(_numpy_array("p", p))
    g.GetFieldData().AddArray(_numpy_array("gain", np.array([1.5, 2.5, 3.5])))
    _vtk_write(tmp_path / "vtk.vtkhdf", g)
    m = meshioplusplus.read(tmp_path / "vtk.vtkhdf")
    assert m.cells[0].data.tolist() == [[0, 1, 2, 3], [1, 2, 3, 4]]
    assert np.array_equal(m.point_data["u"], u) and np.array_equal(
        m.point_data["vel"], vel
    )
    assert np.array_equal(m.cell_data["p"][0], p)
    assert np.array_equal(m.field_data["gain"], [1.5, 2.5, 3.5])
    assert m.points.dtype == np.float32  # the stored dtype is kept, not widened


def test_vtk_written_partitioned_transient_sample(tmp_path):
    def build(o, t):
        o.SetNumberOfPartitions(2)
        for i in range(2):
            g = _tets(5.0 * i)
            g.GetPointData().AddArray(_numpy_array("u", np.full(5, t + i)))
            g.GetCellData().AddArray(_numpy_array("p", np.full(2, 10 * (t + i))))
            o.SetPartition(i, g)

    _vtk_write(tmp_path / "pt.vtkhdf", source=_Temporal("vtkPartitionedDataSet", build))
    for k, t in enumerate((0.0, 0.5, 1.0)):
        m = meshioplusplus.read(tmp_path / "pt.vtkhdf", time_step=k)
        assert m.point_data["u"].tolist() == [t] * 5 + [t + 1] * 5
        assert m.cell_data["p"][0].tolist() == [10 * t] * 2 + [10 * (t + 1)] * 2
        assert m.field_data["meshio:time"].tolist() == [t]
        assert [(r.name, r.entries.tolist()) for r in m.regions] == [
            ("piece_0", [0, 1]),
            ("piece_1", [2, 3]),
        ]
        assert m.points[5, 0] == 5.0
    only = meshioplusplus.vtkhdf.read(tmp_path / "pt.vtkhdf", time_step=2, piece=1)
    assert only.point_data["u"].tolist() == [2.0] * 5


def test_vtk_written_polyhedra_and_polydata_samples(tmp_path):
    g = vtkUnstructuredGrid()
    pts = vtk.vtkPoints()
    for q in CUBE + [(0, 0, 2)]:
        pts.InsertNextPoint(*q)
    g.SetPoints(pts)
    g.InsertNextCell(vtk.VTK_TETRA, 4, [4, 5, 6, 8])
    faces = vtk.vtkIdList()
    faces.InsertNextId(6)
    for f in CUBE_FACES:
        faces.InsertNextId(4)
        for i in f:
            faces.InsertNextId(i)
    g.InsertNextCell(vtk.VTK_POLYHEDRON, faces)
    g.InsertNextCell(vtk.VTK_TETRA, 4, [0, 1, 3, 4])
    _vtk_write(tmp_path / "poly.vtkhdf", g)
    m = meshioplusplus.read(tmp_path / "poly.vtkhdf")
    assert [(c.type, len(c)) for c in m.cells] == [
        ("tetra", 1),
        ("polyhedron8", 1),
        ("tetra", 1),
    ]
    assert m.cells[0].data.tolist() == [[4, 5, 6, 8]] and m.cells[2].data.tolist() == [
        [0, 1, 3, 4]
    ]

    pd = vtkPolyData()
    p = vtk.vtkPoints()
    for q in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (2, 2, 2)]:
        p.InsertNextPoint(*q)
    pd.SetPoints(p)
    pd.Allocate(4)
    pd.InsertNextCell(vtk.VTK_TRIANGLE, 3, [0, 1, 2])
    pd.InsertNextCell(vtk.VTK_QUAD, 4, [0, 1, 3, 2])
    pd.InsertNextCell(vtk.VTK_LINE, 2, [0, 4])
    pd.InsertNextCell(vtk.VTK_VERTEX, 1, [4])
    pd.GetCellData().AddArray(_numpy_array("tag", np.arange(4.0)))
    _vtk_write(tmp_path / "pd.vtkhdf", pd)
    back = meshioplusplus.read(tmp_path / "pd.vtkhdf")
    assert [(c.type, len(c)) for c in back.cells] == [
        ("vertex", 1),
        ("line", 1),
        ("triangle", 1),
        ("quad", 1),
    ]
    assert [a.tolist() for a in back.cell_data["tag"]] == [[0.0], [1.0], [2.0], [3.0]]


def test_vtk_written_partitioned_polydata_keeps_cell_data_aligned(tmp_path):
    """VTK stores PolyData CellData piece-major while topology is category-major."""

    def piece(i):
        pd = vtkPolyData()
        p = vtk.vtkPoints()
        for q in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (2, 2, 2)]:
            p.InsertNextPoint(q[0] + 10 * i, q[1], q[2])
        pd.SetPoints(p)
        pd.Allocate(4)
        pd.InsertNextCell(vtk.VTK_TRIANGLE, 3, [0, 1, 2])
        pd.InsertNextCell(vtk.VTK_QUAD, 4, [0, 1, 3, 2])
        pd.InsertNextCell(vtk.VTK_LINE, 2, [0, 4])
        pd.InsertNextCell(vtk.VTK_VERTEX, 1, [4])
        pd.GetCellData().AddArray(
            _numpy_array("tag", np.array([100 * i + j for j in range(4)], dtype=float))
        )
        return pd

    pds = vtkPartitionedDataSet()
    pds.SetNumberOfPartitions(2)
    pds.SetPartition(0, piece(0))
    pds.SetPartition(1, piece(1))
    _vtk_write(tmp_path / "pdp.vtkhdf", pds)
    m = meshioplusplus.read(tmp_path / "pdp.vtkhdf")
    tags = np.concatenate(m.cell_data["tag"])
    kinds = [c.type for c in m.cells for _ in range(len(c))]
    assert kinds == ["vertex", "line", "triangle", "quad"] * 2
    assert tags.tolist() == [0, 1, 2, 3, 100, 101, 102, 103]
    assert [r.entries.tolist() for r in m.regions] == [[0, 1, 2, 3], [4, 5, 6, 7]]


def _series_via(writer, path):
    if writer == "cpp":
        from meshioplusplus import _core

        with _core.VtkhdfTimeSeriesWriter(str(path)) as w:
            grid = meshioplusplus.Mesh(
                np.array(
                    [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float
                ),
                [("tetra", np.array([[0, 1, 2, 3], [1, 2, 3, 4]]))],
            )
            w.write_points_cells(grid)
            for k in range(3):
                step = meshioplusplus.Mesh(
                    grid.points,
                    [("tetra", grid.cells[0].data)],
                    point_data={"u": np.full(5, float(k))},
                    cell_data={"p": [np.array([10.0 * k, 10.0 * k + 1])]},
                )
                w.write_data(0.5 * k, step)
    else:
        from meshioplusplus.vtkhdf import TimeSeriesWriter

        with TimeSeriesWriter(path) as w:
            w.write_points_cells(
                np.array(
                    [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float
                ),
                [("tetra", np.array([[0, 1, 2, 3], [1, 2, 3, 4]]))],
            )
            for k in range(3):
                w.write_data(
                    0.5 * k,
                    point_data={"u": np.full(5, float(k))},
                    cell_data={"p": [np.array([10.0 * k, 10.0 * k + 1])]},
                )


@pytest.mark.parametrize("writer", ["python", "cpp"])
def test_vtk_plays_a_series_written_by_either_engine(writer, tmp_path):
    """The roadmap's 'working time slider': TIME_STEPS, changing fields, one geometry."""
    from meshioplusplus import _core

    if writer == "cpp" and not hasattr(_core, "VtkhdfTimeSeriesWriter"):
        pytest.skip("_core predates VtkhdfTimeSeriesWriter")
    p = tmp_path / "series.vtkhdf"
    _series_via(writer, p)
    reader = _vtk_read(p)
    assert list(reader.GetOutputInformation(0).Get(SDDP.TIME_STEPS())) == [
        0.0,
        0.5,
        1.0,
    ]
    geometry = None
    for k, t in enumerate((0.0, 0.5, 1.0)):
        reader.UpdateTimeStep(t)
        g = reader.GetOutputDataObject(0)
        points = _arr(g.GetPoints().GetData()).copy()
        geometry = points if geometry is None else geometry
        assert np.array_equal(points, geometry)  # geometry x 1: identical at every step
        assert _arr(g.GetPointData().GetArray("u")).tolist() == [float(k)] * 5
        assert _arr(g.GetCellData().GetArray("p")).tolist() == [10.0 * k, 10.0 * k + 1]
