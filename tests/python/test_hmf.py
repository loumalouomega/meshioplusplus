import numpy as np
import pytest

import meshioplusplus

from . import helpers

test_set_full = [
    helpers.empty_mesh,
    helpers.line_mesh,
    helpers.tri_mesh,
    helpers.line_tri_mesh,
    helpers.tri_mesh_2d,
    helpers.triangle6_mesh,
    helpers.quad_mesh,
    helpers.quad8_mesh,
    helpers.tri_quad_mesh,
    helpers.tet_mesh,
    helpers.tet10_mesh,
    helpers.hex_mesh,
    helpers.hex20_mesh,
    helpers.add_point_data(helpers.tri_mesh, 1),
    helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
]


@pytest.mark.parametrize("mesh", test_set_full)
@pytest.mark.parametrize("compression", [None, "gzip"])
def test_xdmf3(mesh, compression, tmp_path):
    def write(*args, **kwargs):
        return meshioplusplus.xdmf.write(*args, compression=compression, **kwargs)

    helpers.write_read(tmp_path, write, meshioplusplus.xdmf.read, mesh, 1.0e-14)


@pytest.mark.skip
def test_generic_io(tmp_path):
    with pytest.warns(UserWarning):
        helpers.generic_io(tmp_path / "test.hmf")

    with pytest.warns(UserWarning):
        # With additional, insignificant suffix:
        helpers.generic_io(tmp_path / "test.0.hmf")


def test_python_reader_splits_cell_data_by_cell_count(tmp_path):
    # The h5py reader (what a build without HDF5 uses) split the raw cell data
    # by the length of each block's type *name* -- 8 for "triangle" -- rather
    # than its cell count, so any multi-block mesh with cell data failed.
    pytest.importorskip("h5py")
    from meshioplusplus.hmf import _hmf

    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]),
        [("triangle", [[0, 1, 2], [1, 3, 2]]), ("line", [[0, 1]])],
        cell_data={"c": [np.array([1.5, 2.5]), np.array([3.5])]},
    )
    path = tmp_path / "blocks.hmf"
    _hmf.write(path, mesh)
    back = _hmf.read(path)
    assert [b.type for b in back.cells] == ["triangle", "line"]
    assert [list(d) for d in back.cell_data["c"]] == [[1.5, 2.5], [3.5]]
