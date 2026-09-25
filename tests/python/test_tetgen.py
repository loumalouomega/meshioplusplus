import pathlib

import pytest

import meshioplusplus

from . import helpers

test_set = [
    # helpers.empty_mesh,
    helpers.tet_mesh
]


@pytest.mark.parametrize("mesh", test_set)
def test(mesh, tmp_path):
    helpers.write_read(
        tmp_path,
        meshioplusplus.tetgen.write,
        meshioplusplus.tetgen.read,
        mesh,
        1.0e-15,
        extension=".node",
    )


@pytest.mark.parametrize(
    "filename, point_ref_sum, cell_ref_sum", [("mesh.ele", 12, 373)]
)
def test_point_cell_refs(filename, point_ref_sum, cell_ref_sum):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "tetgen" / filename

    mesh = meshioplusplus.read(filename)
    assert mesh.point_data["tetgen:ref"].sum() == point_ref_sum
    assert mesh.cell_data["tetgen:ref"][0].sum() == cell_ref_sum


@pytest.mark.parametrize("empty", ["node", "ele"])
def test_python_reader_refuses_an_empty_file_instead_of_looping(tmp_path, empty):
    # The header scan used to loop for ever at EOF (readline() returns "").
    from meshioplusplus.tetgen._tetgen import read as py_read

    (tmp_path / "m.node").write_text("" if empty == "node" else "1 3 0 0\n0 0 0 0\n")
    (tmp_path / "m.ele").write_text("" if empty == "ele" else "1 4 0\n0 0 0 0 0\n")
    with pytest.raises(meshioplusplus.ReadError, match="no header line"):
        py_read(tmp_path / "m.node")
