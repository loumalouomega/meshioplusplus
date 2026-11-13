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
        helpers.quad_mesh,
        helpers.tri_quad_mesh,
        helpers.polygon_mesh,
    ],
)
def test_obj(mesh, tmp_path):
    for k, c in enumerate(mesh.cells):
        mesh.cells[k] = meshioplusplus.CellBlock(c.type, c.data.astype(np.int32))

    helpers.write_read(
        tmp_path, meshioplusplus.obj.write, meshioplusplus.obj.read, mesh, 1.0e-12
    )


@pytest.mark.skip("Fails point data consistency check.")
@pytest.mark.parametrize(
    "filename, ref_sum, ref_num_cells", [("elephav.obj", 3.678372172450000e05, 1148)]
)
def test_reference_file(filename, ref_sum, ref_num_cells):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "obj" / filename

    mesh = meshioplusplus.read(filename)
    tol = 1.0e-5
    s = np.sum(mesh.points)
    assert abs(s - ref_sum) < tol * abs(ref_sum)
    assert mesh.cells[0].type == "triangle"
    assert len(mesh.cells[0].data) == ref_num_cells


_CUBE_FACE_NORMALS = """\
v 0 0 0
v 1 0 0
v 1 1 0
v 0 1 0
v 0 0 1
v 1 0 1
vn 0 0 -1
vn 0 -1 0
f 1//1 4//1 3//1 2//1
f 1//2 2//2 6//2 5//2
"""


def test_per_corner_normals_split_points(tmp_path):
    # vertices 1 and 2 are shared by faces with different normals: they split.
    path = tmp_path / "n.obj"
    path.write_text(_CUBE_FACE_NORMALS)
    mesh = meshioplusplus.read(path)
    nrm = mesh.point_data["obj:vn"]
    assert len(mesh.points) == len(nrm) == 8
    parent = mesh.point_data["obj:parent_point"]
    assert parent[:6].tolist() == [0, 1, 2, 3, 4, 5]
    for block in mesh.cells:
        for cell in block.data:
            n = {tuple(nrm[i]) for i in cell}
            assert len(n) == 1  # one normal per face
    assert {tuple(n) for n in nrm[~np.isnan(nrm).any(axis=1)]} == {
        (0.0, 0.0, -1.0),
        (0.0, -1.0, 0.0),
    }


def test_one_normal_per_vertex_is_unchanged(tmp_path):
    path = tmp_path / "n.obj"
    path.write_text(
        "v 0 0 0\nv 1 0 0\nv 0 1 0\nvn 0 0 1\nvn 0 0 1\nvn 0 0 1\n" "f 1//1 2//2 3//3\n"
    )
    mesh = meshioplusplus.read(path)
    assert len(mesh.points) == 3
    assert "obj:parent_point" not in mesh.point_data
    assert mesh.point_data["obj:vn"].shape == (3, 3)
