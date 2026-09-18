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
        # v9.9.0: point data, as `dim="0"` mesh functions. This used to be
        # dropped outright while cell data round-tripped.
        helpers.add_point_data(helpers.tri_mesh, 1, dtype=float),
        helpers.add_point_data(helpers.tet_mesh, 1, dtype=np.int32),
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


def test_point_data_is_written_as_a_dim_zero_mesh_function(tmp_path):
    """`dim` is the entity dimension, so vertices are 0 -- the discriminator."""
    mesh = helpers.add_point_data(helpers.tri_mesh, 1, dtype=float)
    out = tmp_path / "pd.xml"
    meshioplusplus.dolfin.write(out, mesh)
    name = next(iter(mesh.point_data))
    sibling = tmp_path / f"pd_{name}.xml"
    assert sibling.exists(), "no sibling mesh_function file was written"
    assert 'dim="0"' in sibling.read_text()


def test_point_and_cell_data_of_the_same_name_do_not_clobber(tmp_path):
    """Both want the same sibling file; cell data owns it, the point one warns."""
    mesh = meshioplusplus.Mesh(
        helpers.tri_mesh.points.copy(),
        [(cb.type, cb.data.copy()) for cb in helpers.tri_mesh.cells],
    )
    ncells = sum(len(cb.data) for cb in mesh.cells)
    mesh.cell_data["x"] = [np.arange(ncells, dtype=float)]
    mesh.point_data["x"] = np.zeros(len(mesh.points))
    out = tmp_path / "clash.xml"
    meshioplusplus.dolfin.write(out, mesh)

    back = meshioplusplus.dolfin.read(out)
    assert "x" in back.cell_data, "the cell array must be the one that survives"
    assert np.array_equal(back.cell_data["x"][0], np.arange(ncells, dtype=float))
    assert "x" not in back.point_data


# --- roadmap §1 "DOLFIN writes every cell-data block to the same sibling
# file" ---


def _two_triangle_block_mesh():
    # Two distinct triangle cell blocks, each with its own cell-data values,
    # so a write that truncates the sibling file per block is caught: block
    # 0's values (10, 11) would be overwritten by block 1's (20) if the
    # writer reopened "<stem>_region.xml" once per block.
    mesh = meshioplusplus.Mesh(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0], [2, 1, 0]],
        [("triangle", [[0, 1, 2], [0, 2, 3]]), ("triangle", [[1, 4, 5]])],
        cell_data={"region": [np.array([10.0, 11.0]), np.array([20.0])]},
    )
    return mesh


def test_two_cell_data_blocks_both_survive_in_the_sibling_file(tmp_path):
    mesh = _two_triangle_block_mesh()
    out = tmp_path / "two_blocks.xml"
    meshioplusplus.dolfin.write(out, mesh)

    back = meshioplusplus.dolfin.read(out)
    assert "region" in back.cell_data
    assert (
        back.cell_data["region"][0].shape[0] == 3
    ), "the second block's row was overwritten or lost"
    assert np.array_equal(back.cell_data["region"][0], [10.0, 11.0, 20.0])


def test_two_cell_data_blocks_python_engine(tmp_path):
    # The pure-Python engine (_dolfin.py) has its own independent copy of
    # this bug; exercise it directly rather than only through the
    # C++-accelerated default dispatch.
    from meshioplusplus.dolfin._dolfin import read as _py_read
    from meshioplusplus.dolfin._dolfin import write as _py_write

    mesh = _two_triangle_block_mesh()
    out = tmp_path / "two_blocks_py.xml"
    _py_write(str(out), mesh)

    back = _py_read(str(out))
    assert "region" in back.cell_data
    assert back.cell_data["region"][0].shape[0] == 3
    assert np.array_equal(back.cell_data["region"][0], [10.0, 11.0, 20.0])
