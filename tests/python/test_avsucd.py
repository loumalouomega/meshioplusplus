import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.empty_mesh,
        helpers.tri_mesh,
        helpers.quad_mesh,
        helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.hex_mesh,
        helpers.add_cell_data(
            helpers.tri_mesh,
            [("avsucd:material", (), int), ("a", (), float), ("b", (3,), float)],
        ),
        helpers.add_point_data(helpers.add_point_data(helpers.tri_mesh, 1), 3),
    ],
)
def test(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.avsucd.write, meshioplusplus.avsucd.read, mesh, 1.0e-13
    )


# --- roadmap §1 "writers that drop data without a warning" ---


def _two_int_cell_data_mesh():
    ncells = len(helpers.tri_mesh.cells[0].data)
    mesh = meshioplusplus.Mesh(
        helpers.tri_mesh.points.copy(),
        [(cb.type, cb.data.copy()) for cb in helpers.tri_mesh.cells],
        cell_data={
            "first_material": [np.arange(ncells, dtype=np.int32)],
            "second_material": [np.arange(ncells, dtype=np.int32) * 2],
        },
    )
    return mesh


def test_warns_when_demoting_an_extra_integer_cell_data_array(tmp_path, capfd):
    # AVS-UCD can only carry one integer "material id" column; the other
    # integer array is not dropped -- it's demoted to real-valued cell data
    # -- but that used to happen with no diagnostic at all.
    mesh = _two_int_cell_data_mesh()
    meshioplusplus.avsucd.write(tmp_path / "two_int.avs", mesh)
    err = capfd.readouterr().err
    assert "first_material" in err
    assert "second_material" in err
    # The array is demoted (still written, as real-valued data), not
    # skipped -- the old message's "skipping" wording was itself wrong.
    assert "skip" not in err.lower()


def test_warns_when_demoting_python_engine(tmp_path, capsys):
    # meshio++ format writers emit these warnings via rich to stderr (not
    # warnings.warn), so they are captured with capsys rather than
    # pytest.warns -- see tests/python/test_mesh.py's same convention.
    from meshioplusplus.avsucd._avsucd import write as _py_write

    mesh = _two_int_cell_data_mesh()
    _py_write(str(tmp_path / "two_int_py.avs"), mesh)
    err = capsys.readouterr().err
    assert "first_material" in err
    assert "second_material" in err
    assert "skip" not in err.lower()
