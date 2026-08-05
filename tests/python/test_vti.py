"""VTK XML ImageData.

The writer requires a lattice, so the fixture set is small on purpose: there is
exactly one shape this format can hold. What is worth testing instead is that the
two implementations agree on that shape's *geometry* -- ``Origin``, ``Spacing``
and ``WholeExtent`` are three attributes that a reader and a writer sharing a
convention can both get wrong together -- and that hand-written files whose
expected geometry is stated independently read correctly.
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus._grid import lattice_from_mesh
from meshioplusplus.vti import _vti


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
    path = tmp_path / "g.vti"
    meshioplusplus.vti.write(path, mesh, binary=binary, compression=compression)
    back = meshioplusplus.vti.read(path)
    assert np.allclose(back.points, mesh.points, atol=1e-12)
    assert len(back.cells) == 1
    assert np.array_equal(back.cells[0].data, mesh.cells[0].data)
    assert np.allclose(back.point_data["f"], mesh.point_data["f"], atol=1e-9)
    assert np.array_equal(back.cell_data["tag"][0], mesh.cell_data["tag"][0])


@pytest.mark.parametrize(
    "binary, compression", [(False, None), (True, None), (True, "zlib")]
)
def test_cross_compat(binary, compression, tmp_path):
    """Each writer's bytes read by the *other* implementation.

    The real cross-check on the base64 block framing: ``_vti._encode_binary`` is
    transcribed from VTU's writer closures, and this is what proves the
    transcription rather than a comment claiming it.
    """
    mesh = _lattice()
    cpp = tmp_path / "cpp.vti"
    meshioplusplus.vti.write(cpp, mesh, binary=binary, compression=compression)
    from_py = _vti.read(cpp)
    assert np.allclose(from_py.points, mesh.points, atol=1e-12)
    assert np.allclose(from_py.point_data["f"], mesh.point_data["f"], atol=1e-9)

    py = tmp_path / "py.vti"
    _vti.write(py, mesh, binary=binary, compression=compression)
    from_cpp = meshioplusplus.vti.read(py)
    assert np.allclose(from_cpp.points, mesh.points, atol=1e-12)
    assert np.allclose(from_cpp.point_data["f"], mesh.point_data["f"], atol=1e-9)
    assert np.array_equal(from_cpp.cell_data["tag"][0], mesh.cell_data["tag"][0])


def test_generic_io(tmp_path):
    """Registration: reachable through ``meshioplusplus.read``/``write``."""
    mesh = _lattice(point_data=False, cell_data=False)
    path = tmp_path / "g.vti"
    meshioplusplus.write(path, mesh)
    back = meshioplusplus.read(path)
    assert np.allclose(back.points, mesh.points)
    assert meshioplusplus.sniff_format(path) == "vti"


def test_reads_a_foreign_file_with_a_shifted_extent(tmp_path):
    """A hand-written file whose expected geometry is stated here, not produced.

    An extent that does not start at zero puts the lo corner at
    ``Origin + start * Spacing``; dropping that offset translates the whole grid
    by one cell, silently and plausibly.
    """
    path = tmp_path / "foreign.vti"
    path.write_text(
        '<?xml version="1.0"?>\n'
        '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n'
        '<ImageData WholeExtent="2 4 0 1 0 1" Origin="1 2 3" Spacing="0.5 2 4">\n'
        '<Piece Extent="2 4 0 1 0 1">\n'
        '<PointData><DataArray type="Int32" Name="i" format="ascii">\n'
        "0 1 2 3 4 5 6 7 8 9 10 11\n</DataArray></PointData>\n"
        "</Piece></ImageData></VTKFile>\n"
    )
    for reader in (meshioplusplus.vti.read, _vti.read):
        m = reader(path)
        assert len(m.points) == 12
        assert len(m.cells) == 1 and len(m.cells[0].data) == 2
        assert np.allclose(m.points[0], [2.0, 2.0, 3.0])
        assert np.allclose(m.points[1], [2.5, 2.0, 3.0])  # x fastest
        assert np.allclose(m.points[3], [2.0, 4.0, 3.0])  # next y row
        assert np.array_equal(m.point_data["i"], np.arange(12))


def test_refuses_a_mesh_that_is_not_a_lattice(tmp_path):
    path = tmp_path / "no.vti"
    tetra = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    with pytest.raises(meshioplusplus.WriteError):
        _vti.write(path, tetra)
    # A PARTIAL lattice is the interesting refusal: right cells, right places,
    # and ImageData still cannot express its holes.
    surf = meshioplusplus.extract_surface(meshioplusplus.grid([2, 2, 2]))
    partial = meshioplusplus.voxelize(surf, resolution=(4, 4, 4), fill="surface")
    assert lattice_from_mesh(partial) is None
    with pytest.raises(meshioplusplus.WriteError):
        _vti.write(path, partial)


@pytest.mark.parametrize(
    "text",
    [
        # a piece that is not the whole extent
        '<VTKFile type="ImageData"><ImageData WholeExtent="0 2 0 1 0 1">'
        '<Piece Extent="0 1 0 1 0 1"/></ImageData></VTKFile>',
        # a rotated lattice
        '<VTKFile type="ImageData"><ImageData WholeExtent="0 1 0 1 0 1" '
        'Direction="0 1 0 -1 0 0 0 0 1"><Piece Extent="0 1 0 1 0 1"/>'
        "</ImageData></VTKFile>",
        # the wrong dataset type
        '<VTKFile type="UnstructuredGrid"><UnstructuredGrid/></VTKFile>',
    ],
)
def test_declines_what_it_does_not_implement(text, tmp_path):
    path = tmp_path / "bad.vti"
    path.write_text(text)
    with pytest.raises(meshioplusplus.ReadError):
        _vti.read(path)


def test_the_recovered_header_equals_the_written_one(tmp_path):
    """The identity the format exists for.

    No format persists arbitrary ``field_data``, so a grid written anywhere else
    loses its ``sdf:*`` header. ``.vti`` stores the geometry itself, and
    ``lattice_from_mesh`` reads it back off the points -- so the round trip is
    exact rather than approximate.
    """
    mesh = meshioplusplus.grid([5, 5, 5], (1.25,) * 3, (0.125,) * 3)
    before = lattice_from_mesh(mesh)
    assert before is not None
    path = tmp_path / "g.vti"
    meshioplusplus.vti.write(path, mesh)
    after = lattice_from_mesh(meshioplusplus.vti.read(path))
    assert after is not None
    for a, b in zip(after, before):
        assert np.allclose(a, b, rtol=0, atol=1e-15)


def test_lattice_from_mesh_cpp_matches_python():
    """The recovery agrees with ``detail::lattice_from_mesh``.

    Checked through the .vti writer, which is the C++ side's only consumer: if
    the two disagreed about a mesh's origin or spacing, the file this writes
    would place the grid somewhere the Python twin does not.
    """
    import os
    import tempfile

    mesh = meshioplusplus.grid([3, 4, 5], (-1.5, 0.25, 2.0), (0.5, 0.25, 0.125))
    dims, origin, spacing = lattice_from_mesh(mesh)
    fd, path = tempfile.mkstemp(suffix=".vti")
    os.close(fd)
    try:
        meshioplusplus.vti.write(path, mesh, binary=False)
        text = open(path).read()
    finally:
        os.unlink(path)
    assert f'WholeExtent="0 {dims[0]} 0 {dims[1]} 0 {dims[2]}"' in text

    def num(v):
        return f"{float(v):.17g}"

    assert f'Origin="{num(origin[0])} {num(origin[1])} {num(origin[2])}"' in text
    assert f'Spacing="{num(spacing[0])} {num(spacing[1])} {num(spacing[2])}"' in text

    # And the Python writer writes the identical attribute line.
    fd, path = tempfile.mkstemp(suffix=".vti")
    os.close(fd)
    try:
        _vti.write(path, mesh, binary=False)
        py_text = open(path).read()
    finally:
        os.unlink(path)
    marker = "<ImageData "
    assert (
        text[text.index(marker) : text.index(">", text.index(marker))]
        == py_text[py_text.index(marker) : py_text.index(">", py_text.index(marker))]
    )


def test_lattice_from_mesh_rejects_what_is_not_one():
    tetra = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    assert lattice_from_mesh(tetra) is None

    # A permuted lattice: identical plane sets, different mesh. A recovery that
    # only looked at the distinct coordinate values would accept this.
    m = meshioplusplus.grid([2, 2, 2])
    permuted = meshioplusplus.Mesh(
        m.points[::-1].copy(), [("hexahedron", m.cells[0].data)]
    )
    assert lattice_from_mesh(permuted) is None

    # A graded grid: uniform plane counts, non-uniform spacing.
    xs = np.array([0.0, 1.0, 3.0])
    pts = np.array(
        [
            [xs[i], float(j), float(k)]
            for k in range(3)
            for j in range(3)
            for i in range(3)
        ]
    )
    graded = meshioplusplus.Mesh(pts, [("hexahedron", m.cells[0].data)])
    assert lattice_from_mesh(graded) is None
