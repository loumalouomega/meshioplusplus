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


@pytest.mark.parametrize("engine", ["core", "python"])
def test_smesh_surface_round_trips(engine, tmp_path):
    """Roadmap 7.7.2: the ``.smesh`` surface TetGen meshes, in both engines,
    each reading what the other wrote."""
    import numpy as np

    from meshioplusplus.tetgen import _tetgen

    points = (
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float)
        + np.pi * 1e-3
    )
    mesh = meshioplusplus.Mesh(
        points,
        [
            ("triangle", [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]]),
            ("quad", [[1, 2, 4, 3]]),
        ],
        cell_data={"tetgen:ref": [np.array([1, 2, 3, 4]), np.array([5])]},
    )
    path = tmp_path / "surface.smesh"
    if engine == "core":
        meshioplusplus.write(path, mesh)
        back = _tetgen.read(path)
    else:
        _tetgen.write(path, mesh)
        back = meshioplusplus.read(path)
    assert np.array_equal(back.points, points)
    assert [c.type for c in back.cells] == ["triangle", "quad"]
    assert np.array_equal(back.cells[1].data, [[1, 2, 4, 3]])
    assert np.array_equal(back.cell_data["tetgen:ref"][0], [1, 2, 3, 4])
    assert np.array_equal(back.cell_data["tetgen:ref"][1], [5])


def test_smesh_refuses_volume_cells_and_bad_facets(tmp_path):
    import numpy as np

    from meshioplusplus.tetgen import _tetgen

    tet = meshioplusplus.Mesh(np.eye(4, 3), [("tetra", [[0, 1, 2, 3]])])
    for write in (meshioplusplus.write, _tetgen.write):
        with pytest.raises(meshioplusplus.WriteError, match="tetra"):
            write(tmp_path / "bad.smesh", tet)
    bad = tmp_path / "range.smesh"
    bad.write_text("3 3 0 0\n0 0 0 0\n1 1 0 0\n2 0 1 0\n1 0\n3 0 1 9\n0\n0\n")
    for read in (meshioplusplus.read, _tetgen.read):
        with pytest.raises(meshioplusplus.ReadError):
            read(bad)
