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


def test_non_flat_3d_raises(tmp_path):
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        [("triangle", [[0, 1, 2]])],
    )
    filepath = tmp_path / "out.svg"
    with pytest.raises(meshioplusplus.WriteError):
        meshioplusplus.svg.write(filepath, mesh)
