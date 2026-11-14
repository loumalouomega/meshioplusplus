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


@pytest.mark.parametrize("engine", ["core", "python"])
def test_materials_and_a_texture_round_trip(tmp_path, engine, monkeypatch):
    from meshioplusplus import _fallback, _png
    from meshioplusplus._regions import Region
    from meshioplusplus.obj import _obj

    image = np.full((2, 2, 4), 200, dtype=np.uint8)
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
        point_data={"uv": np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)},
        field_data={
            "texture:png": np.frombuffer(_png.encode(image), np.uint8),
            "obj:material:red:Kd": np.array([1.0, 0.0, 0.0]),
        },
        regions=[Region("red", "cell", [0])],
    )
    path = tmp_path / "part.obj"
    if engine == "core":
        monkeypatch.setattr(_fallback, "_strict", True)
        meshioplusplus.obj.write(path, mesh)
    else:
        _obj.write(path, mesh)
    text = path.read_text()
    assert (
        "mtllib part.mtl" in text
        and "usemtl red" in text
        and "usemtl material0" in text
    )
    assert (tmp_path / "part.png").read_bytes() == _png.encode(image)
    for read in (meshioplusplus.obj.read, _obj.read):
        back = read(path)
        regions = {r.name: np.asarray(r.entries).tolist() for r in back.regions}
        assert regions == {"red": [0], "material0": [1]}
        assert bytes(back.field_data["texture:png"]) == _png.encode(image)
        assert np.allclose(back.field_data["obj:material:red:Kd"], [1, 0, 0])
