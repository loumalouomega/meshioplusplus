import pytest

import meshioplusplus

from . import helpers

test_set = [
    helpers.empty_mesh,
    helpers.line_mesh,
    helpers.tri_mesh,
    helpers.tri_mesh_2d,
    helpers.quad_mesh,
]


def _drawable_cell_count(mesh):
    return sum(
        len(cb.data) for cb in mesh.cells if cb.type in ("line", "triangle", "quad")
    )


@pytest.mark.parametrize("mesh", test_set)
def test(mesh, tmp_path):
    filepath = tmp_path / "out.tikz"
    meshioplusplus.write_points_cells(filepath, mesh.points, mesh.cells)

    content = filepath.read_text()
    assert "\\documentclass{standalone}" in content
    assert "\\begin{tikzpicture}" in content
    assert "\\end{tikzpicture}" in content
    assert content.count("\\draw") == _drawable_cell_count(mesh)


@pytest.mark.parametrize("mesh", test_set)
def test_cpp_matches_python(mesh, tmp_path):
    # TikZ is plain text: the C++ core writer must be byte-identical to the
    # pure-Python reference.
    cpp = tmp_path / "cpp.tikz"
    py = tmp_path / "py.tikz"
    meshioplusplus._core.tikz_write(str(cpp), mesh)
    meshioplusplus.tikz._tikz.write(str(py), mesh)
    assert cpp.read_text() == py.read_text()


def test_standalone_false(tmp_path):
    mesh = helpers.tri_mesh_2d
    filepath = tmp_path / "snippet.tikz"
    meshioplusplus.tikz.write(filepath, mesh, standalone=False)

    content = filepath.read_text()
    assert "\\documentclass" not in content
    assert "\\begin{tikzpicture}" in content
    assert "\\draw" in content


def test_non_flat_3d_raises(tmp_path):
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        [("triangle", [[0, 1, 2]])],
    )
    filepath = tmp_path / "out.tikz"
    with pytest.raises(meshioplusplus.WriteError):
        meshioplusplus.tikz.write(filepath, mesh)
