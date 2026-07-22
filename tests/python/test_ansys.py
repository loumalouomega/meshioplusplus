import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.empty_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.quad_mesh,
        helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.hex_mesh,
        helpers.pyramid_mesh,
        helpers.wedge_mesh,
    ],
)
@pytest.mark.parametrize("binary", [False, True])
def test(mesh, binary, tmp_path):
    def writer(*args, **kwargs):
        return meshioplusplus.ansys.write(*args, binary=binary, **kwargs)

    helpers.write_read(tmp_path, writer, meshioplusplus.ansys.read, mesh, 1.0e-15)


# --- malformed-input / error-path coverage ---


def test_ansys_garbage_raises(tmp_path):
    p = tmp_path / "bad.msh"
    p.write_text("this is not an ansys mesh file\n")
    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.read(p, file_format="ansys")
