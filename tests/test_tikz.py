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


@pytest.mark.parametrize("mesh", test_set)
def test(mesh, tmp_path):
    filepath = tmp_path / "out.tikz"
    meshioplusplus.write_points_cells(filepath, mesh.points, mesh.cells)

    content = filepath.read_text()
    # Default is a full, compilable standalone document.
    assert "\\documentclass{standalone}" in content
    assert "\\begin{tikzpicture}" in content
    assert "\\end{tikzpicture}" in content


def test_standalone_false(tmp_path):
    mesh = helpers.tri_mesh_2d
    filepath = tmp_path / "snippet.tikz"
    meshioplusplus.tikz.write(filepath, mesh, standalone=False)

    content = filepath.read_text()
    # Snippet-only: no document wrapper, just the tikzpicture environment.
    assert "\\documentclass" not in content
    assert "\\begin{tikzpicture}" in content
    assert "\\draw" in content


def test_non_flat_3d_raises(tmp_path):
    # A genuinely non-flat 3D mesh (nonzero z) cannot be drawn in 2D.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        [("triangle", [[0, 1, 2]])],
    )
    filepath = tmp_path / "out.tikz"
    with pytest.raises(meshioplusplus.WriteError):
        meshioplusplus.tikz.write(filepath, mesh)
