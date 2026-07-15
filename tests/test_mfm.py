import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.line_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.quad_mesh,
        helpers.tet_mesh,
        helpers.hex_mesh,
        helpers.wedge_mesh,
    ],
)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.mfm.write, meshioplusplus.mfm.read, mesh, 1.0e-12
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.mfm")
    helpers.generic_io(tmp_path / "test.0.mfm")


def test_reject_mixed(tmp_path):
    with pytest.raises(meshioplusplus.WriteError):
        meshioplusplus.mfm.write(tmp_path / "x.mfm", helpers.tri_quad_mesh)
