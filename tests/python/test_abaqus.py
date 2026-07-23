import pathlib

import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.empty_mesh,
        helpers.tri_mesh,
        helpers.triangle6_mesh,
        helpers.quad_mesh,
        helpers.quad8_mesh,
        helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
    ],
)
def test(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.abaqus.write, meshioplusplus.abaqus.read, mesh, 1.0e-15
    )


@pytest.mark.parametrize(
    "filename, ref_sum, ref_num_cells, ref_num_cell_sets",
    [
        ("UUea.inp", 4950.0, 50, 10),
        ("nle1xf3c.inp", 32.215275528, 12, 3),
        ("element_elset.inp", 6.0, 2, 3),
        # The `*ELSET` lives in the `*INCLUDE`d file. It used to be dropped —
        # `_abaqus.py`'s `merge()` carried points and cells across but left
        # cell sets behind under a TODO. Both readers now carry it, so the
        # expected count is 1 rather than 0.
        ("wInclude_main.inp", 1.5, 2, 1),
    ],
)
def test_reference_file(filename, ref_sum, ref_num_cells, ref_num_cell_sets):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "abaqus" / filename

    mesh = meshioplusplus.read(filename)

    assert np.isclose(np.sum(mesh.points), ref_sum)
    assert sum(len(cells.data) for cells in mesh.cells) == ref_num_cells
    assert len(mesh.cell_sets) == ref_num_cell_sets


def test_elset(tmp_path):
    points = np.array(
        [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.5, 0.0], [0.0, 0.5, 0.0]]
    )
    cells = [
        ("triangle", np.array([[0, 1, 2]])),
        ("triangle", np.array([[0, 1, 3]])),
    ]
    cell_sets = {
        "right": [np.array([0]), np.array([])],
        "left": [np.array([]), np.array([1])],
    }
    mesh_ref = meshioplusplus.Mesh(points, cells, cell_sets=cell_sets)

    filepath = tmp_path / "test.inp"
    meshioplusplus.abaqus.write(filepath, mesh_ref)
    mesh = meshioplusplus.abaqus.read(filepath)

    assert np.allclose(mesh_ref.points, mesh.points)

    assert len(mesh_ref.cells) == len(mesh.cells)
    for ic, cell in enumerate(mesh_ref.cells):
        assert cell.type == mesh.cells[ic].type
        assert np.allclose(cell.data, mesh.cells[ic].data)

    assert sorted(mesh_ref.cell_sets.keys()) == sorted(mesh.cell_sets.keys())
    for k, v in mesh_ref.cell_sets.items():
        for ic in range(len(mesh_ref.cells)):
            assert np.allclose(v[ic], mesh.cell_sets[k][ic])
