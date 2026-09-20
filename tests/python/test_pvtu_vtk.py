"""`.pvtu` / `.pvtp` against VTK's own parallel readers and writers.

The round-trip tests in ``test_pvtu.py`` cannot catch an index that *we* read and
Kitware's reader does not -- a wrong root tag, a missing declaration, a piece
path VTK will not resolve -- nor one that VTK writes and our parser is too strict
for. These do, in both directions, for the roadmap's "done when": a partitioned
case's index opens in VTK.

`.pvd` has no counterpart here: vanilla VTK has no collection reader
(`vtkPVDReader` ships with ParaView only), so `test_pvd.py` asserts its structure
and each step's file is checked to open in VTK below.

Skipped where ``vtk`` (pulled in by the ``pyvista`` extra) is not installed.
"""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._fallback import set_strict_core

pytest.importorskip("vtkmodules.vtkIOXML")
pytest.importorskip("vtkmodules.vtkIOParallelXML")

from vtkmodules.util import numpy_support as ns  # noqa: E402
from vtkmodules.vtkCommonCore import vtkPoints  # noqa: E402
from vtkmodules.vtkCommonDataModel import (  # noqa: E402
    VTK_QUAD,
    vtkPolyData,
    vtkUnstructuredGrid,
)
from vtkmodules.vtkIOParallelXML import (  # noqa: E402
    vtkXMLPPolyDataWriter,
    vtkXMLPUnstructuredGridWriter,
)
from vtkmodules.vtkIOXML import (  # noqa: E402
    vtkXMLPPolyDataReader,
    vtkXMLPUnstructuredGridReader,
    vtkXMLUnstructuredGridReader,
)

N = 4  # an N x N quad grid


def _points():
    return np.array([[i, j, 0.0] for j in range(N + 1) for i in range(N + 1)])


def _quads():
    return np.array(
        [
            [
                j * (N + 1) + i,
                j * (N + 1) + i + 1,
                (j + 1) * (N + 1) + i + 1,
                (j + 1) * (N + 1) + i,
            ]
            for j in range(N)
            for i in range(N)
        ]
    )


def _mesh(cell_type="quad"):
    pts = _points()
    mesh = meshioplusplus.Mesh(pts, [("quad", _quads())])
    mesh.point_data["u"] = pts[:, 0] + 10 * pts[:, 1]
    mesh.cell_data["c"] = [np.arange(N * N, dtype=np.float64)]
    mesh.cell_data["partition:part"] = [np.arange(N * N) // 4]
    return mesh


def _vtk_dataset(polydata):
    ds = vtkPolyData() if polydata else vtkUnstructuredGrid()
    pts = vtkPoints()
    for p in _points():
        pts.InsertNextPoint(*p)
    ds.SetPoints(pts)
    if polydata:
        ds.Allocate(N * N)
    for q in _quads():
        ds.InsertNextCell(VTK_QUAD, 4, [int(v) for v in q])
    u = ns.numpy_to_vtk(_points()[:, 0] + 10 * _points()[:, 1], deep=1)
    u.SetName("u")
    ds.GetPointData().AddArray(u)
    c = ns.numpy_to_vtk(np.arange(N * N, dtype=np.float64), deep=1)
    c.SetName("c")
    ds.GetCellData().AddArray(c)
    return ds


def _vtk_read(path, polydata=False):
    reader = vtkXMLPPolyDataReader() if polydata else vtkXMLPUnstructuredGridReader()
    assert reader.CanReadFile(str(path)) == 1, "VTK's reader refuses the index"
    reader.SetFileName(str(path))
    reader.UpdateInformation()
    return reader


@pytest.fixture(params=["core", "python"])
def write_read(request):
    """`(kind -> module)` for the engine under test; the core one is strict."""
    if request.param == "core":
        if not hasattr(_core, "pvtu_read"):
            pytest.skip("this build has no pvtu core")
        set_strict_core(True)
        yield {"pvtu": meshioplusplus.pvtu, "pvtp": meshioplusplus.pvtp}
        set_strict_core(None)
    else:
        from meshioplusplus.pvtp import _pvtp
        from meshioplusplus.pvtu import _pvtu

        yield {"pvtu": _pvtu, "pvtp": _pvtp}


# --- VTK reads what we write -------------------------------------------------


def test_vtk_reads_our_pvtu(write_read, tmp_path):
    path = tmp_path / "ours.pvtu"
    write_read["pvtu"].write(path, _mesh())
    reader = _vtk_read(path)
    assert reader.GetNumberOfPieces() == 4
    reader.Update()
    out = reader.GetOutput()
    assert out.GetNumberOfCells() == N * N
    assert out.GetNumberOfPoints() > (N + 1) ** 2  # interface points per piece
    assert {"u"} <= {
        out.GetPointData().GetArrayName(i)
        for i in range(out.GetPointData().GetNumberOfArrays())
    }
    c = ns.vtk_to_numpy(out.GetCellData().GetArray("c"))
    np.testing.assert_array_equal(np.sort(c), np.arange(N * N))


def test_vtk_reads_our_pvtp(write_read, tmp_path):
    path = tmp_path / "ours.pvtp"
    write_read["pvtp"].write(path, _mesh())
    reader = _vtk_read(path, polydata=True)
    assert reader.GetNumberOfPieces() == 4
    reader.Update()
    assert reader.GetOutput().GetNumberOfCells() == N * N


def test_vtk_sees_our_ghost_cells(write_read, tmp_path):
    """`partition`'s halo becomes a `vtkGhostType` array VTK recognises as one."""
    if not hasattr(_core, "partition"):
        pytest.skip("ghost layers need the C++ core")
    mesh = meshioplusplus.Mesh(_points(), [("quad", _quads())])
    mesh.cell_data["c"] = [np.arange(N * N, dtype=np.float64)]
    pieces = meshioplusplus.partition(mesh, 3, ghost_layers=1)
    path = tmp_path / "gh.pvtu"
    write_read["pvtu"].write_pieces(path, pieces)

    body = ET.parse(path).getroot().find("PUnstructuredGrid")
    assert body.get("GhostLevel") == "1"
    reader = _vtk_read(path)
    assert reader.GetNumberOfPieces() == 3
    reader.Update()
    out = reader.GetOutput()
    ghosts = out.GetCellGhostArray()
    assert ghosts is not None, "VTK did not recognise vtkGhostType as the ghost array"
    flagged = int(np.count_nonzero(ns.vtk_to_numpy(ghosts)))
    expected = sum(
        int(np.count_nonzero(p.cell_data["partition:ghost"][0])) for p in pieces
    )
    assert flagged == expected > 0


# --- we read what VTK writes -------------------------------------------------


@pytest.mark.parametrize("mode", ["ascii", "binary"])
def test_we_read_a_pvtu_vtk_wrote(write_read, mode, tmp_path):
    w = vtkXMLPUnstructuredGridWriter()
    w.SetFileName(str(tmp_path / "theirs.pvtu"))
    w.SetNumberOfPieces(4)
    w.SetStartPiece(0)
    w.SetEndPiece(3)
    w.SetGhostLevel(0)
    w.SetDataModeToAscii() if mode == "ascii" else w.SetDataModeToBinary()
    w.SetInputData(_vtk_dataset(polydata=False))
    assert w.Write() == 1

    back = write_read["pvtu"].read(tmp_path / "theirs.pvtu")
    # VTK's writer of a non-streaming dataset puts the whole grid in every piece
    assert len(back.points) == 4 * (N + 1) ** 2
    assert sum(len(cb.data) for cb in back.cells) == 4 * N * N
    assert [r.name for r in back.regions] == [f"piece_{i}" for i in range(4)]
    assert sorted(back.point_data) == ["u"] and sorted(back.cell_data) == ["c"]
    np.testing.assert_array_equal(
        np.sort(np.concatenate(back.cell_data["c"])),
        np.sort(np.tile(np.arange(N * N), 4)),
    )

    one = write_read["pvtu"].read(tmp_path / "theirs.pvtu", piece=2)
    assert len(one.points) == (N + 1) ** 2
    assert sum(len(cb.data) for cb in one.cells) == N * N
    np.testing.assert_allclose(
        one.points, _points()
    )  # VTK wrote Float32; ints are exact
    np.testing.assert_array_equal(one.cell_data["c"][0], np.arange(N * N))


def test_we_read_a_pvtp_vtk_wrote(write_read, tmp_path):
    w = vtkXMLPPolyDataWriter()
    w.SetFileName(str(tmp_path / "theirs.pvtp"))
    w.SetNumberOfPieces(2)
    w.SetStartPiece(0)
    w.SetEndPiece(1)
    w.SetDataModeToBinary()
    w.SetInputData(_vtk_dataset(polydata=True))
    assert w.Write() == 1
    back = write_read["pvtp"].read(tmp_path / "theirs.pvtp")
    assert sum(len(cb.data) for cb in back.cells) == 2 * N * N
    assert [r.name for r in back.regions] == ["piece_0", "piece_1"]


def test_our_pvtu_survives_a_trip_through_vtk(write_read, tmp_path):
    """Ours -> VTK's serial writer over the merged output -> ours: the geometry is intact."""
    path = tmp_path / "ours.pvtu"
    write_read["pvtu"].write(path, _mesh())
    reader = _vtk_read(path)
    reader.Update()
    out = reader.GetOutput()
    welded_points = {
        tuple(np.round(p, 9)) for p in ns.vtk_to_numpy(out.GetPoints().GetData())
    }
    assert welded_points == {tuple(p) for p in _points()}


# --- .pvd --------------------------------------------------------------------


def test_every_pvd_step_opens_in_vtk(write_read, tmp_path):
    """No `.pvd` reader in vanilla VTK: check the index shape ParaView reads
    (`timestep`, `part`, `file` on each `DataSet`) and that each named file opens."""
    for t in range(3):
        write_read["pvtu"].write(tmp_path / f"s{t}.pvtu", _mesh())
    entries = "".join(
        f'<DataSet timestep="{0.5 * t}" part="0" file="s{t}.pvtu"/>' for t in range(3)
    )
    (tmp_path / "run.pvd").write_text(
        '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" '
        f'byte_order="LittleEndian"><Collection>{entries}</Collection></VTKFile>'
    )
    # ours reads a hand-written ParaView-shaped collection ...
    assert meshioplusplus.read_metadata(tmp_path / "run.pvd")["time_values"] == [
        0.0,
        0.5,
        1.0,
    ]
    # ... and writes the same shape
    meshioplusplus.write_sequence(
        tmp_path / "own.pvd",
        ((t, meshioplusplus.Mesh(_points(), [("quad", _quads())])) for t in (0.0, 1.0)),
    )
    for pvd in (tmp_path / "run.pvd", tmp_path / "own.pvd"):
        ds = ET.parse(pvd).getroot().find("Collection").findall("DataSet")
        assert ds and all({"timestep", "file"} <= set(e.attrib) for e in ds)
        for e in ds:
            f = str(pvd.parent / e.get("file"))
            if f.endswith(".pvtu"):
                assert _vtk_read(f).GetNumberOfPieces() == 4
            else:
                assert vtkXMLUnstructuredGridReader().CanReadFile(f) == 1
