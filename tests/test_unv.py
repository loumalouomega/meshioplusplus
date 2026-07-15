import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.line_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.triangle6_mesh,
        helpers.quad_mesh,
        helpers.quad8_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
        helpers.wedge_mesh,
    ],
)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.unv.write, meshioplusplus.unv.read, mesh, 1.0e-12
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.unv")
    helpers.generic_io(tmp_path / "test.0.unv")


def test_groups(tmp_path):
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    mesh.point_sets = {"corners": np.array([0, 2])}
    mesh.cell_sets = {"all": [np.array([0, 1])]}
    p = tmp_path / "g.unv"
    meshioplusplus.unv.write(p, mesh)
    out = meshioplusplus.unv.read(p)
    assert np.array_equal(out.point_sets["corners"], [0, 2])
    assert np.array_equal(out.cell_sets["all"][0], [0, 1])
