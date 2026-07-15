import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tri_mesh_2d,
        helpers.tri_mesh,
        helpers.tet_mesh,
    ],
)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path,
        meshioplusplus.freefem.write,
        meshioplusplus.freefem.read,
        mesh,
        1.0e-12,
        ".msh",
    )


def test_explicit_file_format(tmp_path):
    p = tmp_path / "test.msh"
    meshioplusplus.freefem.write(p, helpers.tri_mesh_2d)
    mesh = meshioplusplus.read(p, file_format="freefem")
    assert mesh.cells[0].type == "triangle"
    assert len(mesh.cells[0].data) == 2
