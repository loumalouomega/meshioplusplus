import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [helpers.empty_mesh, helpers.tri_mesh],
)
@pytest.mark.parametrize(
    "binary, tol",
    [
        (False, 1.0e-15),
        # binary STL only operates in single precision
        (True, 1.0e-8),
    ],
)
def test_stl(mesh, binary, tol, tmp_path):
    def writer(*args, **kwargs):
        return meshioplusplus.stl.write(*args, binary=binary, **kwargs)

    helpers.write_read(tmp_path, writer, meshioplusplus.stl.read, mesh, tol)


@pytest.mark.parametrize("binary", [False, True])
@pytest.mark.parametrize(
    "mesh, num_triangles",
    [
        (helpers.tet_mesh, 6),  # 2 tets, shared face removed
        (helpers.hex_mesh, 12),  # 6 boundary quads, triangulated
        (helpers.wedge_mesh, 8),  # 2 tris + 3 quads -> 2 + 6
        (helpers.pyramid_mesh, 6),  # 4 tris + 1 quad -> 4 + 2
    ],
)
def test_stl_volume_mesh_writes_skin(mesh, num_triangles, binary, tmp_path):
    import copy

    filepath = tmp_path / "skin.stl"
    meshioplusplus.stl.write(filepath, copy.deepcopy(mesh), binary=binary)
    out = meshioplusplus.stl.read(filepath)
    assert len(out.get_cells_type("triangle")) == num_triangles


def test_stl_skin_false_legacy(tmp_path):
    import copy

    # skin=False keeps the legacy behavior: volume cells are discarded with a
    # warning and the STL comes out empty.
    filepath = tmp_path / "legacy.stl"
    meshioplusplus.stl.write(filepath, copy.deepcopy(helpers.tet_mesh), skin=False)
    out = meshioplusplus.stl.read(filepath)
    assert len(out.cells) == 0


def test_stl_python_fallback_matches_skin(tmp_path):
    import copy

    # The pure-Python reference writer performs the same skin extraction.
    from meshioplusplus.stl import _stl

    filepath = tmp_path / "py_skin.stl"
    _stl.write(filepath, copy.deepcopy(helpers.hex_mesh))
    out = meshioplusplus.stl.read(filepath)
    assert len(out.get_cells_type("triangle")) == 12


def test_binary_attribute_word_round_trip(tmp_path):
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], float)
    attr = np.array([0, 7], dtype=np.uint16)
    mesh = meshioplusplus.Mesh(
        pts,
        [("triangle", np.array([[0, 1, 2], [1, 3, 2]]))],
        cell_data={"stl:attribute": [attr]},
    )
    path = tmp_path / "a.stl"
    meshioplusplus.write(path, mesh, file_format="stl", binary=True)
    back = meshioplusplus.read(path)
    assert np.array_equal(back.cell_data["stl:attribute"][0], attr)
    # all-zero attributes are not kept
    plain = meshioplusplus.Mesh(pts, mesh.cells)
    meshioplusplus.write(path, plain, file_format="stl", binary=True)
    assert "stl:attribute" not in meshioplusplus.read(path).cell_data


_TWO_SOLIDS = """solid left
facet normal 0 0 1
 outer loop
  vertex 0 0 0
  vertex 1 0 0
  vertex 0 1 0
 endloop
endfacet
endsolid left
solid right
facet normal 0 0 1
 outer loop
  vertex 1 0 0
  vertex 1 1 0
  vertex 0 1 0
 endloop
endfacet
facet normal 0 0 1
 outer loop
  vertex 1 1 0
  vertex 2 1 0
  vertex 1 2 0
 endloop
endfacet
endsolid right
"""


def test_several_solids_are_cell_regions(tmp_path):
    path = tmp_path / "two.stl"
    path.write_text(_TWO_SOLIDS)
    mesh = meshioplusplus.read(path)
    by_name = {r.name: r.entries.tolist() for r in mesh.regions}
    assert by_name == {"left": [0], "right": [1, 2]}
    out = tmp_path / "out.stl"
    meshioplusplus.write(out, mesh, file_format="stl")
    text = out.read_text()
    assert "solid left" in text and "endsolid right" in text
    back = meshioplusplus.read(out)
    assert {r.name: r.entries.tolist() for r in back.regions} == by_name


def test_one_solid_adds_no_region(tmp_path):
    path = tmp_path / "one.stl"
    path.write_text(_TWO_SOLIDS.split("solid right")[0])
    assert not meshioplusplus.read(path).regions
