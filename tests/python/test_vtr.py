"""VTK XML RectilinearGrid (roadmap §1 tier B4, v11.6.0).

Unlike ``.vts``, the reader here is genuinely more capable than ``.vti``'s
uniform lattice: per-axis coordinates need not be evenly spaced at all, so a
graded grid (finer near a wall) reads correctly. The writer still requires a
UNIFORM dense lattice (no ``rectilinear_from_mesh`` detector exists yet --
documented in the format's own doc page).
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus._grid import lattice_from_mesh
from meshioplusplus.vtr import _vtr


def _lattice(n=3, origin=-0.5, spacing=0.25, point_data=True, cell_data=True):
    m = meshioplusplus.grid([n, n, n], (origin,) * 3, (spacing,) * 3)
    if point_data:
        m.point_data["f"] = np.arange(len(m.points), dtype=np.float64) * 0.5
    if cell_data:
        m.cell_data["tag"] = [np.arange(n**3, dtype=np.int64)]
    return m


@pytest.mark.parametrize(
    "binary, compression", [(False, None), (True, None), (True, "zlib")]
)
def test_roundtrip(binary, compression, tmp_path):
    mesh = _lattice()
    path = tmp_path / "g.vtr"
    meshioplusplus.vtr.write(path, mesh, binary=binary, compression=compression)
    back = meshioplusplus.vtr.read(path)
    assert np.allclose(back.points, mesh.points, atol=1e-12)
    assert len(back.cells) == 1
    assert np.array_equal(back.cells[0].data, mesh.cells[0].data)
    assert np.allclose(back.point_data["f"], mesh.point_data["f"], atol=1e-9)
    assert np.array_equal(back.cell_data["tag"][0], mesh.cell_data["tag"][0])


@pytest.mark.parametrize(
    "binary, compression", [(False, None), (True, None), (True, "zlib")]
)
def test_cross_compat(binary, compression, tmp_path):
    mesh = _lattice()
    cpp = tmp_path / "cpp.vtr"
    meshioplusplus.vtr.write(cpp, mesh, binary=binary, compression=compression)
    from_py = _vtr.read(cpp)
    assert np.allclose(from_py.points, mesh.points, atol=1e-12)
    assert np.allclose(from_py.point_data["f"], mesh.point_data["f"], atol=1e-9)

    py = tmp_path / "py.vtr"
    _vtr.write(py, mesh, binary=binary, compression=compression)
    from_cpp = meshioplusplus.vtr.read(py)
    assert np.allclose(from_cpp.points, mesh.points, atol=1e-12)
    assert np.allclose(from_cpp.point_data["f"], mesh.point_data["f"], atol=1e-9)
    assert np.array_equal(from_cpp.cell_data["tag"][0], mesh.cell_data["tag"][0])


def test_generic_io(tmp_path):
    mesh = _lattice(point_data=False, cell_data=False)
    path = tmp_path / "g.vtr"
    meshioplusplus.write(path, mesh)
    back = meshioplusplus.read(path)
    assert np.allclose(back.points, mesh.points)
    assert meshioplusplus.sniff_format(path) == "vtr"


def test_reads_a_graded_grid(tmp_path):
    """The identity `.vtr` exists for: per-axis spacing `.vti` could not
    express at all. Neither engine checks uniformity on read."""
    path = tmp_path / "graded.vtr"
    path.write_text(
        '<?xml version="1.0"?>\n'
        '<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">\n'
        '<RectilinearGrid WholeExtent="0 2 0 1 0 1">\n'
        '<Piece Extent="0 2 0 1 0 1">\n'
        "<Coordinates>\n"
        '<DataArray type="Float64" Name="x_coordinates" format="ascii">\n'
        "0 1 4\n</DataArray>\n"
        '<DataArray type="Float64" Name="y_coordinates" format="ascii">\n'
        "0 10\n</DataArray>\n"
        '<DataArray type="Float64" Name="z_coordinates" format="ascii">\n'
        "0 100\n</DataArray>\n"
        "</Coordinates>\n"
        "</Piece></RectilinearGrid></VTKFile>\n"
    )
    for reader in (meshioplusplus.vtr.read, _vtr.read):
        m = reader(path)
        assert len(m.points) == 12
        assert len(m.cells) == 1 and len(m.cells[0].data) == 2
        assert np.allclose(m.points[2], [4.0, 0.0, 0.0])
        assert np.allclose(m.points[3], [0.0, 10.0, 0.0])
        assert np.allclose(m.points[6], [0.0, 0.0, 100.0])


def test_refuses_a_mesh_that_is_not_a_uniform_lattice(tmp_path):
    path = tmp_path / "no.vtr"
    tetra = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    with pytest.raises(meshioplusplus.WriteError):
        _vtr.write(path, tetra)
    surf = meshioplusplus.extract_surface(meshioplusplus.grid([2, 2, 2]))
    partial = meshioplusplus.voxelize(surf, resolution=(4, 4, 4), fill="surface")
    assert lattice_from_mesh(partial) is None
    with pytest.raises(meshioplusplus.WriteError):
        _vtr.write(path, partial)


@pytest.mark.parametrize(
    "text",
    [
        '<VTKFile type="RectilinearGrid"><RectilinearGrid WholeExtent="0 2 0 1 0 1">'
        '<Piece Extent="0 1 0 1 0 1"/></RectilinearGrid></VTKFile>',
        '<VTKFile type="UnstructuredGrid"><UnstructuredGrid/></VTKFile>',
    ],
)
def test_declines_what_it_does_not_implement(text, tmp_path):
    path = tmp_path / "bad.vtr"
    path.write_text(text)
    with pytest.raises(meshioplusplus.ReadError):
        _vtr.read(path)


def test_metadata_agrees_with_a_real_read(tmp_path):
    mesh = _lattice(point_data=False, cell_data=False)
    path = tmp_path / "g.vtr"
    meshioplusplus.vtr.write(path, mesh, binary=True)
    meta = meshioplusplus.read_metadata(path)
    back = meshioplusplus.vtr.read(path)
    assert meta["num_points"] == len(back.points)
    assert meta["cell_blocks"][0]["type"] == "hexahedron"
    assert meta["cell_blocks"][0]["num_cells"] == len(back.cells[0].data)
