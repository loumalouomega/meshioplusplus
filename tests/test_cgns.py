import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        # helpers.empty_mesh,
        helpers.tet_mesh
    ],
)
def test(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.cgns.write, meshioplusplus.cgns.read, mesh, 1.0e-15
    )
