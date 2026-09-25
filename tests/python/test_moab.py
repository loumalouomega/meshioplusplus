import pytest

import meshioplusplus

from . import helpers

h5py = pytest.importorskip("h5py")


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.empty_mesh,
        helpers.line_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.tet_mesh,
    ],
)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.h5m.write, meshioplusplus.h5m.read, mesh, 1.0e-15
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.h5m")
    # With additional, insignificant suffix:
    helpers.generic_io(tmp_path / "test.0.h5m")


def test_python_writer_writes_cell_data_per_block(tmp_path, capsys):
    # The Python writer (the core declines a mesh with types h5m cannot hold)
    # treated `cell_data` as {cell type: {name: array}} and crashed on the
    # {name: [one array per block]} it actually is.
    import numpy as np

    from meshioplusplus.h5m._h5m import write as py_write

    mesh = meshioplusplus.Mesh(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
        [
            ("triangle", [[0, 1, 2]]),
            ("quad", [[0, 1, 2, 3]]),
            ("tetra", [[0, 1, 2, 3]]),
        ],
        cell_data={"c": [np.array([1.0]), np.array([2.0]), np.array([3.0])]},
    )
    path = tmp_path / "m.h5m"
    py_write(path, mesh)
    assert "'quad'" in capsys.readouterr().err  # named, not a literal '%s'
    with h5py.File(path) as f:
        elements = f["tstt"]["elements"]
        assert elements["Tri3"]["tags"]["c"][()].tolist() == [1.0]
        assert elements["Tet4"]["tags"]["c"][()].tolist() == [3.0]
