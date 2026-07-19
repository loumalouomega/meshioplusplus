import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        # helpers.empty_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.tet_mesh,
        helpers.add_cell_data(
            helpers.tri_mesh, [("a", (), float), ("b", (), np.int64)]
        ),
    ],
)
def test_dolfin(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.dolfin.write, meshioplusplus.dolfin.read, mesh, 1.0e-15
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.xml")
    # With additional, insignificant suffix:
    helpers.generic_io(tmp_path / "test.0.xml")


# --- malformed-input / error-path coverage (Python reference reader) ---
from meshioplusplus.dolfin._dolfin import read as _dolfin_py_read  # noqa: E402


def test_dolfin_vertices_before_mesh_raises(tmp_path):
    p = tmp_path / "bad.xml"
    p.write_text('<dolfin><vertices size="1"/></dolfin>')
    with pytest.raises(meshioplusplus.ReadError):
        _dolfin_py_read(str(p))
