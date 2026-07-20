"""Optional VTK XML block codecs (zstd, lz4).

Guarded so the suite stays meaningful in a build with neither: the
availability, default-behaviour and error-message tests always run, and only
the round-trips are skipped.

The guarantee these exist to protect is the one that makes adding codecs safe
at all: **a pure build reads and writes exactly what it always did.** zlib
stays the default everywhere, and nothing about an existing file changes.
"""

from __future__ import annotations

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus.vtu import _vtu as _vtu_py

_core = pytest.importorskip("meshioplusplus._core")

HAS_ZLIB = getattr(_core, "__has_zlib__", False)
HAS_ZSTD = getattr(_core, "__has_zstd__", False)
HAS_LZ4 = getattr(_core, "__has_lz4__", False)

# The compressor= attribute each codec is recorded under. lz4 is a real VTK
# compressor; zstd is a meshio++ extension (VTK ships no ZSTD compressor).
ATTR = {
    "zlib": "vtkZLibDataCompressor",
    "lz4": "vtkLZ4DataCompressor",
    "zstd": "vtkZSTDDataCompressor",
}

# Unlike zstd/lz4 (default OFF), zlib defaults ON in CMake -- but it is still
# found via find_package, same as HDF5/netCDF, so a system without it (e.g.
# Windows CI, see CMakeLists.txt's MESHIOPLUSPLUS_WITH_ZLIB comment) compiles
# it out too. `meshioplusplus.write`/`.read` never notice (the Python fallback
# picks up zlib from the stdlib); only tests that call `_core` directly need
# this guard.
requires_zlib = pytest.mark.skipif(not HAS_ZLIB, reason="build has no zlib")
requires_zstd = pytest.mark.skipif(not HAS_ZSTD, reason="build has no zstd")
requires_lz4 = pytest.mark.skipif(not HAS_LZ4, reason="build has no lz4")


def _mesh():
    return meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
        point_data={"u": np.arange(4.0), "v": np.arange(4.0) * 3},
    )


def _assert_same(got, src):
    np.testing.assert_allclose(np.asarray(got.points)[:, :3], src.points)
    assert sorted(got.point_data) == sorted(src.point_data)
    for name, arr in src.point_data.items():
        np.testing.assert_allclose(np.asarray(got.point_data[name]), arr)


def test_zlib_capability_flag_exists():
    """Present as a boolean even when off (find_package can fail, e.g. Windows
    CI with no system zlib), matching the zstd/lz4 flags -- callers can branch
    without hasattr regardless of build configuration."""
    assert isinstance(getattr(_core, "__has_zlib__", None), bool)


def test_capability_flags_exist():
    # Present as booleans even when off, so callers can branch without hasattr.
    assert isinstance(_core.__has_zstd__, bool)
    assert isinstance(_core.__has_lz4__, bool)


@pytest.mark.parametrize("fmt", ["vtu", "vtp"])
def test_default_write_is_unchanged(tmp_path, fmt):
    """No codec asked for => zlib, exactly as before."""
    path = tmp_path / f"m.{fmt}"
    meshioplusplus.write(path, _mesh())
    text = path.read_text()
    assert ATTR["lz4"] not in text
    assert ATTR["zstd"] not in text
    _assert_same(meshioplusplus.read(path), _mesh())


@pytest.mark.parametrize(
    "codec",
    [
        "zlib",
        pytest.param("lz4", marks=requires_lz4),
        pytest.param("zstd", marks=requires_zstd),
    ],
)
@pytest.mark.parametrize("fmt", ["vtu", "vtp"])
def test_codec_round_trip(tmp_path, fmt, codec):
    path = tmp_path / f"m.{fmt}"
    src = _mesh()
    meshioplusplus.write(path, src, compression=codec)
    assert ATTR[codec] in path.read_text()
    _assert_same(meshioplusplus.read(path), src)


@pytest.mark.parametrize(
    "codec",
    [
        pytest.param("zlib", marks=requires_zlib),
        pytest.param("lz4", marks=requires_lz4),
        pytest.param("zstd", marks=requires_zstd),
    ],
)
def test_cpp_and_python_implementations_agree(tmp_path, codec):
    """Cross-compat both directions -- the real proof the framings match.

    A C++-written file must be readable by the pure-Python reference and vice
    versa; otherwise the two would silently be writing different formats under
    the same compressor= name.
    """
    src = _mesh()

    cpp_path = tmp_path / "cpp.vtu"
    _core.vtu_write_codec(str(cpp_path), src, True, codec)
    _assert_same(_vtu_py.read(cpp_path), src)

    py_path = tmp_path / "py.vtu"
    _vtu_py.write(py_path, src, binary=True, compression=codec)
    _assert_same(_core.vtu_read(str(py_path)), src)


def test_zlib_file_reads_in_a_build_without_optional_codecs(tmp_path):
    """A zlib file must never depend on zstd/lz4 being present."""
    path = tmp_path / "m.vtu"
    meshioplusplus.write(path, _mesh(), compression="zlib")
    assert ATTR["zlib"] in path.read_text()
    _assert_same(meshioplusplus.read(path), _mesh())


def test_unknown_codec_is_rejected(tmp_path):
    with pytest.raises(Exception, match="codec"):
        _core.vtu_write_codec(str(tmp_path / "m.vtu"), _mesh(), True, "bogus")


def _forge(tmp_path, codec):
    """A minimal VTU declaring `codec` without any real payload."""
    path = tmp_path / "forged.vtu"
    path.write_text(
        f'<VTKFile type="UnstructuredGrid" compressor="{ATTR[codec]}">'
        '<UnstructuredGrid><Piece NumberOfPoints="0"/>'
        "</UnstructuredGrid></VTKFile>"
    )
    return path


@pytest.mark.parametrize("codec", ["lz4", "zstd"])
def test_absent_codec_error_names_the_build_option(tmp_path, codec):
    """If the C++ core lacks a codec, the error must be actionable."""
    have = {"lz4": HAS_LZ4, "zstd": HAS_ZSTD}[codec]
    if have:
        pytest.skip(f"this build has {codec}")
    with pytest.raises(Exception) as excinfo:
        _core.vtu_read(str(_forge(tmp_path, codec)))
    assert "MESHIOPLUSPLUS_WITH_" in str(excinfo.value)


@pytest.mark.parametrize("codec", ["lz4", "zstd"])
def test_python_fallback_reports_the_missing_extra(codec):
    """With neither the C++ codec nor the Python module, fail by name.

    That combination is a genuinely new failure class -- a file readable by
    nothing in the stack -- so the message has to name the extra to install
    rather than surfacing an ImportError from deep inside the parser. Checked
    against the error path directly so it is exercised even when the modules
    *are* installed (as they are in CI).
    """
    shim = {"lz4": _vtu_py._lz4, "zstd": _vtu_py._zstd}[codec]
    message = str(shim._missing())
    assert "meshioplusplus[codecs]" in message
    assert f"MESHIOPLUSPLUS_WITH_{codec.upper()}=ON" in message
    assert {"lz4": "lz4", "zstd": "zstandard"}[codec] in message


def test_cli_compress_with_codec(tmp_path):
    from meshioplusplus._cli import main

    path = tmp_path / "m.vtu"
    meshioplusplus.write(path, _mesh())
    codec = "lz4" if HAS_LZ4 else "zlib"
    main(["compress", "--codec", codec, str(path)])
    assert ATTR[codec] in path.read_text()
    _assert_same(meshioplusplus.read(path), _mesh())


def test_cli_compress_rejects_codec_for_other_formats(tmp_path):
    """Silently ignoring --codec would be the worst outcome."""
    from meshioplusplus._cli import main

    path = tmp_path / "m.stl"
    meshioplusplus.write(path, _mesh())
    with pytest.raises(SystemExit):
        main(["compress", "--codec", "zstd", str(path)])


# ---------------------------------------------------------------------------
# Interop with a real VTK, which is the only thing that can actually settle it
# ---------------------------------------------------------------------------


def _vtk_grid_from(mesh):
    """Build the equivalent vtkUnstructuredGrid, for VTK-side writing."""
    import vtk

    points = vtk.vtkPoints()
    for p in mesh.points:
        points.InsertNextPoint(*p)
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    for block in mesh.cells:
        for row in block.data:
            cell = vtk.vtkTriangle()
            for i, node in enumerate(row):
                cell.GetPointIds().SetId(i, int(node))
            grid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
    return grid


def test_vtk_ships_lz4_but_not_zstd():
    """The fact that motivates treating the two codecs differently.

    lz4 is a real VTK compressor, so meshio++ writes it for interop. VTK has no
    ZSTD compressor, so `vtkZSTDDataCompressor` is explicitly a meshio++
    extension rather than a standard name we happen to support.
    """
    vtk = pytest.importorskip("vtk")
    assert hasattr(vtk, "vtkLZ4DataCompressor")
    assert not hasattr(vtk, "vtkZSTDDataCompressor")


@requires_lz4
def test_vtk_reads_our_lz4_output(tmp_path):
    """Our lz4 files must open in VTK/ParaView -- the whole point of using it."""
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import vtk_to_numpy

    src = _mesh()
    path = tmp_path / "ours.vtu"
    _core.vtu_write_codec(str(path), src, True, "lz4")
    assert ATTR["lz4"] in path.read_text()

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()

    assert grid.GetNumberOfPoints() == len(src.points)
    np.testing.assert_allclose(vtk_to_numpy(grid.GetPoints().GetData()), src.points)


@requires_lz4
def test_we_read_vtks_lz4_output(tmp_path):
    """And the converse: VTK's own raw-block framing must decode here."""
    vtk = pytest.importorskip("vtk")

    src = _mesh()
    path = tmp_path / "fromvtk.vtu"
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(_vtk_grid_from(src))
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToLZ4()
    writer.Write()
    assert ATTR["lz4"] in path.read_text(errors="ignore")

    back = _core.vtu_read(str(path))
    np.testing.assert_allclose(np.asarray(back.points)[:, :3], src.points)


@requires_zstd
def test_vtk_cleanly_refuses_our_zstd_output(tmp_path):
    """zstd must fail *visibly* in VTK, never be silently misread.

    That is the trade being made by writing a non-VTK compressor name: ParaView
    reports an unknown compressor rather than producing wrong geometry.
    """
    vtk = pytest.importorskip("vtk")

    path = tmp_path / "ours.vtu"
    _core.vtu_write_codec(str(path), _mesh(), True, "zstd")

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.SetGlobalWarningDisplay(0)  # the failure is expected; keep output clean
    reader.Update()
    # VTK cannot construct vtkZSTDDataCompressor, so no data comes through.
    assert reader.GetOutput().GetNumberOfPoints() == 0
