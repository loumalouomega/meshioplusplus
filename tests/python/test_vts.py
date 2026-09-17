"""VTK XML StructuredGrid (roadmap §1 tier B4, v11.6.0).

Unlike ``.vti``, the writer's dense-lattice requirement is not the reader's:
``.vts`` states connectivity implicitly (the same index formula) but points
explicitly, so a genuinely curved structured mesh reads correctly even
though it could never be written back (there is no lattice to recover).
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus._grid import lattice_from_mesh
from meshioplusplus.vts import _vts


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
    path = tmp_path / "g.vts"
    meshioplusplus.vts.write(path, mesh, binary=binary, compression=compression)
    back = meshioplusplus.vts.read(path)
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
    cpp = tmp_path / "cpp.vts"
    meshioplusplus.vts.write(cpp, mesh, binary=binary, compression=compression)
    from_py = _vts.read(cpp)
    assert np.allclose(from_py.points, mesh.points, atol=1e-12)
    assert np.allclose(from_py.point_data["f"], mesh.point_data["f"], atol=1e-9)

    py = tmp_path / "py.vts"
    _vts.write(py, mesh, binary=binary, compression=compression)
    from_cpp = meshioplusplus.vts.read(py)
    assert np.allclose(from_cpp.points, mesh.points, atol=1e-12)
    assert np.allclose(from_cpp.point_data["f"], mesh.point_data["f"], atol=1e-9)
    assert np.array_equal(from_cpp.cell_data["tag"][0], mesh.cell_data["tag"][0])


def test_generic_io(tmp_path):
    mesh = _lattice(point_data=False, cell_data=False)
    path = tmp_path / "g.vts"
    meshioplusplus.write(path, mesh)
    back = meshioplusplus.read(path)
    assert np.allclose(back.points, mesh.points)
    assert meshioplusplus.sniff_format(path) == "vts"


def test_reads_a_curved_structured_grid(tmp_path):
    """The identity `.vts` exists for: points `.vti` could not state at all."""
    path = tmp_path / "curved.vts"
    path.write_text(
        '<?xml version="1.0"?>\n'
        '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n'
        '<StructuredGrid WholeExtent="0 1 0 1 0 1">\n'
        '<Piece Extent="0 1 0 1 0 1">\n'
        "<Points><DataArray type=\"Float64\" Name=\"Points\" "
        'NumberOfComponents="3" format="ascii">\n'
        "0 0 0\n1 0 0\n0 1 0\n1 1 0\n0 0 1\n1 0 1\n0 1 1\n1.5 1.5 1.5\n"
        "</DataArray></Points>\n"
        "</Piece></StructuredGrid></VTKFile>\n"
    )
    for reader in (meshioplusplus.vts.read, _vts.read):
        m = reader(path)
        assert len(m.points) == 8
        assert len(m.cells) == 1 and len(m.cells[0].data) == 1
        assert np.allclose(m.points[7], [1.5, 1.5, 1.5])


def test_refuses_a_mesh_that_is_not_a_lattice(tmp_path):
    path = tmp_path / "no.vts"
    tetra = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    with pytest.raises(meshioplusplus.WriteError):
        _vts.write(path, tetra)
    surf = meshioplusplus.extract_surface(meshioplusplus.grid([2, 2, 2]))
    partial = meshioplusplus.voxelize(surf, resolution=(4, 4, 4), fill="surface")
    assert lattice_from_mesh(partial) is None
    with pytest.raises(meshioplusplus.WriteError):
        _vts.write(path, partial)


@pytest.mark.parametrize(
    "text",
    [
        '<VTKFile type="StructuredGrid"><StructuredGrid WholeExtent="0 2 0 1 0 1">'
        '<Piece Extent="0 1 0 1 0 1"/></StructuredGrid></VTKFile>',
        '<VTKFile type="UnstructuredGrid"><UnstructuredGrid/></VTKFile>',
    ],
)
def test_declines_what_it_does_not_implement(text, tmp_path):
    path = tmp_path / "bad.vts"
    path.write_text(text)
    with pytest.raises(meshioplusplus.ReadError):
        _vts.read(path)


def test_metadata_agrees_with_a_real_read(tmp_path):
    mesh = _lattice(point_data=False, cell_data=False)
    path = tmp_path / "g.vts"
    meshioplusplus.vts.write(path, mesh, binary=True)
    meta = meshioplusplus.read_metadata(path)
    back = meshioplusplus.vts.read(path)
    assert meta["num_points"] == len(back.points)
    assert meta["cell_blocks"][0]["type"] == "hexahedron"
    assert meta["cell_blocks"][0]["num_cells"] == len(back.cells[0].data)
