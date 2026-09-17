"""VTK XML MultiBlock (roadmap §1 tier B4, v11.6.0, part 3 of 3).

An index file plus one `.vtu` piece per cell block, combined on read via
`merge()` (no welding) into one mesh with one "cell" region per piece. Unlike
`.vts`/`.vtr`, there is no lattice restriction at all -- any mesh with one or
more cell blocks round-trips.
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus.vtm import _vtm


def _tri_quad_mesh():
    """Three blocks (triangle, quad, triangle) sharing points across
    blocks -- exercises point duplication across pieces and same-type block
    consolidation on the way back in, exactly like the C++ fixture."""
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 1, 0], [2, 1, 0], [1, 1, 0], [0, 1, 0]],
        dtype=np.float64,
    )
    cells = [
        ("triangle", np.array([[0, 1, 5], [0, 5, 6]])),
        ("quad", np.array([[1, 2, 4, 5]])),
        ("triangle", np.array([[2, 3, 4]])),
    ]
    return meshioplusplus.Mesh(points, cells)


@pytest.mark.parametrize(
    "binary, compression", [(False, None), (True, None), (True, "zlib")]
)
def test_roundtrip(binary, compression, tmp_path):
    mesh = _tri_quad_mesh()
    path = tmp_path / "g.vtm"
    meshioplusplus.vtm.write(path, mesh, binary=binary, compression=compression)
    back = meshioplusplus.vtm.read(path)

    assert len(back.points) == 11  # 4 + 4 + 3, points duplicated across pieces
    assert [c.type for c in back.cells] == ["triangle", "quad"]
    assert len(back.cells[0].data) == 3  # the two triangle pieces consolidate
    assert len(back.cells[1].data) == 1

    assert [r.name for r in back.regions] == ["block_0", "block_1", "block_2"]
    assert all(r.kind == "cell" for r in back.regions)
    assert sum(len(r.entries) for r in back.regions) == 4


@pytest.mark.parametrize(
    "binary, compression", [(False, None), (True, None), (True, "zlib")]
)
def test_cross_compat(binary, compression, tmp_path):
    mesh = _tri_quad_mesh()
    cpp = tmp_path / "cpp.vtm"
    meshioplusplus.vtm.write(cpp, mesh, binary=binary, compression=compression)
    from_py = _vtm.read(cpp)
    assert len(from_py.points) == 11
    assert [r.name for r in from_py.regions] == ["block_0", "block_1", "block_2"]

    py = tmp_path / "py.vtm"
    _vtm.write(py, mesh, binary=binary, compression=compression)
    from_cpp = meshioplusplus.vtm.read(py)
    assert len(from_cpp.points) == 11
    assert [c.type for c in from_cpp.cells] == ["triangle", "quad"]


def test_generic_io(tmp_path):
    mesh = _tri_quad_mesh()
    path = tmp_path / "g.vtm"
    meshioplusplus.write(path, mesh)
    back = meshioplusplus.read(path)
    assert len(back.points) == 11
    assert meshioplusplus.sniff_format(path) == "vtm"


def test_pieces_are_independently_readable_vtu_files(tmp_path):
    mesh = _tri_quad_mesh()
    path = tmp_path / "g.vtm"
    meshioplusplus.vtm.write(path, mesh, binary=True)

    stem = path.stem
    expected_cells = [2, 1, 1]
    for i in range(3):
        piece = path.parent / stem / f"{stem}_{i}.vtu"
        assert piece.exists()
        p = meshioplusplus.read(piece)
        assert len(p.cells) == 1
        assert len(p.cells[0].data) == expected_cells[i]


def test_metadata_agrees_with_a_real_read(tmp_path):
    mesh = _tri_quad_mesh()
    path = tmp_path / "g.vtm"
    meshioplusplus.vtm.write(path, mesh, binary=True)

    meta = meshioplusplus.read_metadata(path)
    back = meshioplusplus.vtm.read(path)
    assert meta["num_points"] == len(back.points)
    assert len(meta["cell_blocks"]) == len(back.cells)
    for mb, cb in zip(meta["cell_blocks"], back.cells):
        assert mb["type"] == cb.type
        assert mb["num_cells"] == len(cb.data)


def test_declines_a_piece_that_is_neither_vtu_nor_vtp(tmp_path):
    path = tmp_path / "bad.vtm"
    path.write_text(
        '<?xml version="1.0"?>\n'
        '<VTKFile type="vtkMultiBlockDataSet" version="1.0">\n'
        "<vtkMultiBlockDataSet>\n"
        '<Block index="0">\n'
        '<DataSet index="0" name="a" file="piece.obj"/>\n'
        "</Block>\n"
        "</vtkMultiBlockDataSet>\n"
        "</VTKFile>\n"
    )
    with pytest.raises(meshioplusplus.ReadError):
        _vtm.read(path)


def test_empty_index_reads_as_an_empty_mesh(tmp_path):
    path = tmp_path / "empty.vtm"
    path.write_text(
        '<?xml version="1.0"?>\n'
        '<VTKFile type="vtkMultiBlockDataSet" version="1.0">\n'
        "<vtkMultiBlockDataSet>\n"
        '<Block index="0"/>\n'
        "</vtkMultiBlockDataSet>\n"
        "</VTKFile>\n"
    )
    for reader in (meshioplusplus.vtm.read, _vtm.read):
        back = reader(path)
        assert len(back.points) == 0
        assert len(back.cells) == 0
