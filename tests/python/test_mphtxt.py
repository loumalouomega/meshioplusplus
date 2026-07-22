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
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.wedge_mesh,
        helpers.tri_quad_mesh,
    ],
)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.mphtxt.write, meshioplusplus.mphtxt.read, mesh, 1.0e-12
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.mphtxt")
    helpers.generic_io(tmp_path / "test.0.mphtxt")
