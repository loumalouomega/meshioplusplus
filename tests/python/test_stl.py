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
