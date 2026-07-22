import copy

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus.triangle import _triangle

from . import helpers

triangle6_2d_mesh = meshioplusplus.Mesh(
    [
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.5, 0.0],
        [1.0, 0.5],
        [0.5, 0.5],
    ],
    [("triangle6", [[0, 1, 2, 3, 4, 5]])],
)

data_mesh = meshioplusplus.Mesh(
    helpers.tri_mesh_2d.points,
    helpers.tri_mesh_2d.cells,
    point_data={
        "triangle:attr1": np.array([0.5, 1.25, -2.75, 3.0]),
        "triangle:ref": np.array([1.0, 2.0, 3.0, 4.0]),
    },
    cell_data={"triangle:ref": [np.array([11, 22])]},
)

line_mesh_2d = meshioplusplus.Mesh(
    [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
    [("line", [[0, 1], [1, 2], [2, 3], [3, 0]])],
    cell_data={"triangle:ref": [np.array([1, 2, 3, 4])]},
)


@pytest.mark.parametrize("mesh", [helpers.tri_mesh_2d, triangle6_2d_mesh, data_mesh])
@pytest.mark.parametrize("extension", [".node", ".ele"])
def test_node_ele(mesh, extension, tmp_path):
    helpers.write_read(
        tmp_path,
        meshioplusplus.triangle.write,
        meshioplusplus.triangle.read,
        mesh,
        1.0e-15,
        extension=extension,
    )


def test_poly(tmp_path):
    helpers.write_read(
        tmp_path,
        meshioplusplus.triangle.write,
        meshioplusplus.triangle.read,
        line_mesh_2d,
        1.0e-15,
        extension=".poly",
    )


@pytest.mark.parametrize("mesh", [helpers.tri_mesh_2d, triangle6_2d_mesh])
def test_cross_compat(mesh, tmp_path):
    # shim (C++ where available) write -> pure-Python read
    p = tmp_path / "cross.node"
    meshioplusplus.triangle.write(p, copy.deepcopy(mesh))
    out = _triangle.read(p)
    assert np.allclose(mesh.points, out.points, atol=1.0e-15, rtol=0.0)
    assert np.array_equal(mesh.cells[0].data, out.cells[0].data)

    # pure-Python write -> shim (C++ where available) read
    _triangle.write(p, copy.deepcopy(mesh))
    out = meshioplusplus.triangle.read(p)
    assert np.allclose(mesh.points, out.points, atol=1.0e-15, rtol=0.0)
    assert np.array_equal(mesh.cells[0].data, out.cells[0].data)


def test_lone_node_file(tmp_path):
    # A .node file without an .ele sibling is a valid point cloud.
    p = tmp_path / "cloud.node"
    p.write_text("3 2 0 0\n0 0.0 0.0\n1 1.0 0.0\n2 0.0 1.0\n")
    for reader in (meshioplusplus.triangle.read, _triangle.read):
        mesh = reader(p)
        assert mesh.points.shape == (3, 2)
        assert len(mesh.cells) == 0


def test_extension_dispatch(tmp_path):
    # A 2D .node/.ele pair dispatches to triangle (tetgen raises and the
    # extension loop falls through); a 3D pair still dispatches to tetgen.
    p2 = tmp_path / "two.node"
    meshioplusplus.triangle.write(p2, copy.deepcopy(helpers.tri_mesh_2d))
    mesh = meshioplusplus.read(p2)
    assert mesh.points.shape[1] == 2
    assert mesh.cells[0].type == "triangle"

    p3 = tmp_path / "three.node"
    meshioplusplus.tetgen.write(p3, copy.deepcopy(helpers.tet_mesh))
    mesh = meshioplusplus.read(p3)
    assert mesh.points.shape[1] == 3
    assert mesh.cells[0].type == "tetra"


@pytest.mark.parametrize("writer", [meshioplusplus.triangle.write, _triangle.write])
def test_3d_points_raise(writer, tmp_path):
    with pytest.raises(meshioplusplus.WriteError):
        writer(tmp_path / "bad.node", copy.deepcopy(helpers.tri_mesh))


def test_mixed_triangle_orders_raise(tmp_path):
    mesh = meshioplusplus.Mesh(
        triangle6_2d_mesh.points,
        helpers.tri_mesh_2d.cells + triangle6_2d_mesh.cells,
    )
    with pytest.raises(meshioplusplus.WriteError):
        _triangle.write(tmp_path / "mixed.node", mesh)
