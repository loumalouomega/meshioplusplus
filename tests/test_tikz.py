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


def test_non_flat_shell_is_projected(tmp_path):
    # A non-flat 3D surface mesh is projected (isometric camera) and drawn.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        [("triangle", [[0, 1, 2]])],
    )
    filepath = tmp_path / "out.tikz"
    meshioplusplus.tikz.write(filepath, mesh)
    content = filepath.read_text()
    assert "\\documentclass{standalone}" in content
    assert content.count("\\draw") == 1


@pytest.mark.parametrize(
    "mesh, num_faces",
    [
        (helpers.tet_mesh, 6),
        (helpers.hex_mesh, 6),
        (helpers.wedge_mesh, 5),
        (helpers.pyramid_mesh, 5),
    ],
)
def test_volume_mesh_draws_skin(mesh, num_faces, tmp_path):
    import copy

    filepath = tmp_path / "skin.tikz"
    meshioplusplus.tikz.write(filepath, copy.deepcopy(mesh))
    content = filepath.read_text()
    assert content.count("\\draw") == num_faces
    assert content.count("-- cycle;") == num_faces


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tet_mesh,
        helpers.hex_mesh,
        helpers.wedge_mesh,
        helpers.pyramid_mesh,
    ],
)
def test_cpp_matches_python_3d(mesh, tmp_path):
    import copy

    # The byte-identity guarantee extends to the 3D projected path: skin
    # extraction, camera arithmetic, depth sort, and formatting must all
    # agree between the C++ core and the Python reference.
    cpp = tmp_path / "cpp.tikz"
    py = tmp_path / "py.tikz"
    meshioplusplus._core.tikz_write(str(cpp), copy.deepcopy(mesh))
    meshioplusplus.tikz._tikz.write(str(py), copy.deepcopy(mesh))
    assert cpp.read_text() == py.read_text()


def test_camera_angles_change_projection(tmp_path):
    import copy

    p1 = tmp_path / "iso.tikz"
    p2 = tmp_path / "other.tikz"
    meshioplusplus.tikz.write(p1, copy.deepcopy(helpers.hex_mesh))
    meshioplusplus.tikz.write(
        p2, copy.deepcopy(helpers.hex_mesh), azimuth=10.0, elevation=60.0
    )
    assert p1.read_text() != p2.read_text()
    assert p1.read_text().count("\\draw") == p2.read_text().count("\\draw") == 6
