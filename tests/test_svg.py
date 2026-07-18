from xml.etree import ElementTree as ET

import pytest

import meshioplusplus

from . import helpers

SVG_NS = "{http://www.w3.org/2000/svg}"

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
    filepath = tmp_path / "out.svg"
    meshioplusplus.write_points_cells(filepath, mesh.points, mesh.cells)
    # Output is valid SVG with one <path> per drawable cell.
    paths = ET.parse(filepath).getroot().findall(f"{SVG_NS}path")
    assert len(paths) == _drawable_cell_count(mesh)


@pytest.mark.parametrize("mesh", test_set)
def test_cpp_matches_python(mesh, tmp_path):
    # The C++ core writer and the pure-Python reference agree on path count.
    cpp = tmp_path / "cpp.svg"
    py = tmp_path / "py.svg"
    meshioplusplus._core.svg_write(str(cpp), mesh)
    meshioplusplus.svg._svg.write(str(py), mesh)

    n_cpp = len(ET.parse(cpp).getroot().findall(f"{SVG_NS}path"))
    n_py = len(ET.parse(py).getroot().findall(f"{SVG_NS}path"))
    assert n_cpp == n_py == _drawable_cell_count(mesh)


def _count_paths(filepath):
    return len(ET.parse(filepath).getroot().findall(f"{SVG_NS}path"))


def test_non_flat_shell_is_projected(tmp_path):
    # A non-flat 3D surface mesh is projected (isometric camera) and drawn.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        [("triangle", [[0, 1, 2]])],
    )
    filepath = tmp_path / "out.svg"
    meshioplusplus.svg.write(filepath, mesh)
    assert _count_paths(filepath) == 1


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

    filepath = tmp_path / "skin.svg"
    meshioplusplus.svg.write(filepath, copy.deepcopy(mesh))
    assert _count_paths(filepath) == num_faces


def test_camera_angles_change_projection(tmp_path):
    import copy

    p1 = tmp_path / "iso.svg"
    p2 = tmp_path / "other.svg"
    meshioplusplus.svg.write(p1, copy.deepcopy(helpers.hex_mesh))
    meshioplusplus.svg.write(
        p2, copy.deepcopy(helpers.hex_mesh), azimuth=10.0, elevation=60.0
    )
    assert _count_paths(p1) == _count_paths(p2) == 6
    assert p1.read_bytes() != p2.read_bytes()


def test_cpp_matches_python_3d(tmp_path):
    import copy

    # Path data must agree between the C++ core and the Python reference on
    # the 3D path (the stroke-width default formatting is allowed to differ:
    # C++ %g vs Python str(), as on the flat path).
    cpp = tmp_path / "cpp.svg"
    py = tmp_path / "py.svg"
    meshioplusplus._core.svg_write(str(cpp), copy.deepcopy(helpers.hex_mesh))
    meshioplusplus.svg._svg.write(str(py), copy.deepcopy(helpers.hex_mesh))
    d_cpp = [p.get("d") for p in ET.parse(cpp).getroot().findall(f"{SVG_NS}path")]
    d_py = [p.get("d") for p in ET.parse(py).getroot().findall(f"{SVG_NS}path")]
    assert d_cpp == d_py
