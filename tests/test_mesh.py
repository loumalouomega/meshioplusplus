import copy

import numpy as np
import pytest
from numpy.testing import assert_equal

import meshioplusplus

from . import helpers


def test_cells_dict():
    mesh = copy.deepcopy(helpers.tri_mesh)
    assert len(mesh.cells_dict) == 1
    assert np.array_equal(mesh.cells_dict["triangle"], [[0, 1, 2], [0, 2, 3]])

    # two cells groups
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        [("triangle", [[0, 1, 2]]), ("triangle", [[0, 2, 3]])],
        cell_data={"a": [[0.5], [1.3]]},
    )
    assert len(mesh.cells_dict) == 1
    assert_equal(mesh.cells_dict, {"triangle": [[0, 1, 2], [0, 2, 3]]})
    assert_equal(mesh.cell_data_dict, {"a": {"triangle": [0.5, 1.3]}})


def test_sets_to_int_data():
    mesh = helpers.tri_mesh_5
    mesh = helpers.add_point_sets(mesh)
    mesh = helpers.add_cell_sets(mesh)

    mesh.point_sets_to_data()
    mesh.cell_sets_to_data()

    assert mesh.cell_sets == {}
    assert_equal(mesh.cell_data, {"grain0-grain1": [[0, 0, 1, 1, 1]]})

    assert mesh.point_sets == {}
    assert_equal(mesh.point_data, {"fixed-loose": [0, 0, 0, 1, 1, 1, 1]})

    # now back to set data
    mesh.cell_data_to_sets("grain0-grain1")
    mesh.point_data_to_sets("fixed-loose")

    assert mesh.cell_data == {}
    assert_equal(mesh.cell_sets, {"grain0": [[0, 1]], "grain1": [[2, 3, 4]]})

    assert mesh.point_data == {}
    assert_equal(mesh.point_sets, {"fixed": [0, 1, 2], "loose": [3, 4, 5, 6]})


def test_sets_to_int_data_warning(capsys):
    # meshio++ emits these warnings via rich to stderr (not warnings.warn), so
    # they are captured with capsys rather than pytest.warns.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]],
        {"triangle": [[0, 1, 2], [1, 2, 3]]},
        cell_sets={"tag": [[0]]},
    )
    mesh.cell_sets_to_data()
    assert "not part of any cell set" in capsys.readouterr().err
    assert np.all(mesh.cell_data["tag"] == np.array([[0, -1]]))

    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]],
        {"triangle": [[0, 1, 2], [1, 2, 3]]},
        point_sets={"tag": [[0, 1, 3]]},
    )
    mesh.point_sets_to_data()
    assert "Not all points are part of a point set" in capsys.readouterr().err
    assert np.all(mesh.point_data["tag"] == np.array([[0, 0, -1, 0]]))


def test_int_data_to_sets():
    # Deep-copy the shared fixture: this test mutates cell_data/cell_sets in
    # place, and leaking that onto the module-level helpers.tri_mesh corrupts
    # later round-trip tests (surfaces on the Python-writer path, e.g. Windows).
    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.cell_data = {"grain0-grain1": [np.array([0, 1])]}

    mesh.cell_data_to_sets("grain0-grain1")

    assert_equal(mesh.cell_sets, {"grain0": [[0]], "grain1": [[1]]})


def test_gh_1165():
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        {
            "triangle": [[0, 1, 2], [1, 2, 3]],
            "line": [[0, 1], [0, 2], [1, 3], [2, 3]],
        },
        cell_sets={
            "test": [[], [1]],
            "sets": [[0, 1], [0, 2, 3]],
        },
    )

    mesh.cell_sets_to_data()
    mesh.cell_data_to_sets("test-sets")

    assert_equal(mesh.cell_sets, {"test": [[], [1]], "sets": [[0, 1], [0, 2, 3]]})


def test_copy():
    mesh = helpers.tri_mesh
    mesh2 = mesh.copy()

    assert np.all(mesh.points == mesh2.points)
    assert not np.may_share_memory(mesh.points, mesh2.points)


def test_get_cell_data():
    # Two triangle blocks + a quad block, with per-block cell_data.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [2.0, 0.0]],
        [
            ("triangle", [[0, 1, 2]]),
            ("quad", [[0, 1, 2, 3]]),
            ("triangle", [[1, 4, 2]]),
        ],
        cell_data={"a": [np.array([10.0]), np.array([20.0]), np.array([30.0])]},
    )
    # get_cell_data concatenates the data of all blocks of the given type.
    assert_equal(mesh.get_cell_data("a", "triangle"), [10.0, 30.0])
    assert_equal(mesh.get_cell_data("a", "quad"), [20.0])


def test_init_point_data_length_mismatch():
    with pytest.raises(ValueError):
        meshioplusplus.Mesh(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            [("triangle", [[0, 1, 2]])],
            point_data={"p": [1.0, 2.0, 3.0]},  # 3 != 4 points
        )


def test_init_cell_data_block_count_mismatch():
    with pytest.raises(ValueError):
        meshioplusplus.Mesh(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            [("triangle", [[0, 1, 2]])],  # one cell block
            cell_data={"a": [[1.0], [2.0]]},  # two data blocks
        )


def test_init_cell_data_block_length_mismatch():
    with pytest.raises(ValueError):
        meshioplusplus.Mesh(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            [("triangle", [[0, 1, 2], [0, 2, 3]])],  # 2 cells
            cell_data={"a": [[1.0]]},  # block has 1 value, not 2
        )


def test_topological_dimension():
    from meshioplusplus import topological_dimension

    assert topological_dimension["vertex"] == 0
    assert topological_dimension["line"] == 1
    assert topological_dimension["triangle"] == 2
    assert topological_dimension["quad"] == 2
    assert topological_dimension["tetra"] == 3
    assert topological_dimension["hexahedron"] == 3


def test_cellblock_polygon_dim():
    # A polygon block reports topological dimension 2.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        [("polygon", [[0, 1, 2, 3]])],
    )
    assert mesh.cells[0].dim == 2
