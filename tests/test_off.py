import pathlib

import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        # helpers.empty_mesh,
        helpers.tri_mesh,
        helpers.quad_mesh,
        helpers.polygon_mesh_one_cell,
    ],
)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.off.write, meshioplusplus.off.read, mesh, 1.0e-15
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.off")
    # With additional, insignificant suffix:
    helpers.generic_io(tmp_path / "test.0.off")


def test_quad_faces():
    # https://github.com/loumalouomega/meshioplusplus/issues/35 — OFF files
    # with non-triangular (here: quad) faces used to be rejected outright.
    this_dir = pathlib.Path(__file__).resolve().parent
    mesh = meshioplusplus.read(this_dir / "meshes" / "off" / "cube_example.off")
    assert mesh.points.shape == (8, 3)
    assert len(mesh.cells) == 1
    assert mesh.cells[0].type == "quad"
    assert mesh.cells[0].data.shape == (6, 4)


def test_quad_and_triangulated_cube_agree():
    # Same cube, once as 6 quads and once pre-triangulated into 12 triangles
    # (2 per quad) — both should describe the same geometry.
    this_dir = pathlib.Path(__file__).resolve().parent
    quad = meshioplusplus.read(this_dir / "meshes" / "off" / "cube_example.off")
    tri = meshioplusplus.read(
        this_dir / "meshes" / "off" / "cube_example_as_triangs.off"
    )

    assert np.array_equal(quad.points, tri.points)
    assert quad.cells[0].type == "quad"
    assert quad.cells[0].data.shape == (6, 4)
    assert tri.cells[0].type == "triangle"
    assert tri.cells[0].data.shape == (12, 3)


def test_mixed_face_blocks(tmp_path):
    # Faces are grouped by vertex count into separate triangle/quad/polygon
    # runs; a run boundary is drawn whenever the count changes.
    mesh = meshioplusplus.Mesh(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 1.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ],
        [
            ("triangle", [[0, 1, 2], [0, 1, 4]]),
            ("polygon", [[0, 1, 5, 6, 3]]),
        ],
    )
    helpers.write_read(
        tmp_path, meshioplusplus.off.write, meshioplusplus.off.read, mesh, 1.0e-15
    )
