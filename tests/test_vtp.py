import copy

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus.vtp import _vtp

from . import helpers

test_set = [
    helpers.line_mesh,
    helpers.tri_mesh,
    helpers.tri_mesh_2d,
    helpers.tri_mesh_one_cell,
    helpers.quad_mesh,
    helpers.tri_quad_mesh,
    helpers.polygon_mesh,
    helpers.polygon_mesh_one_cell,
    # NOT polygon2_mesh: PolyData has no cell-type array, so its 3-/4-noded
    # polygon cells legitimately read back as triangle/quad.
    helpers.add_point_data(helpers.tri_mesh, 1),
    helpers.add_point_data(helpers.tri_mesh, 3),
    helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
    helpers.add_cell_data(helpers.tri_quad_mesh, [("a", (), np.float64)]),
]


@pytest.mark.parametrize("mesh", test_set)
@pytest.mark.parametrize(
    "binary, compression", [(False, None), (True, None), (True, "zlib")]
)
def test(mesh, binary, compression, tmp_path):
    def writer(filename, mesh):
        meshioplusplus.vtp.write(filename, mesh, binary=binary, compression=compression)

    tol = 1.0e-15 if binary else 1.0e-10
    helpers.write_read(
        tmp_path, writer, meshioplusplus.vtp.read, mesh, tol, extension=".vtp"
    )


@pytest.mark.parametrize("mesh", [helpers.tri_mesh, helpers.polygon_mesh])
@pytest.mark.parametrize(
    "binary, compression", [(False, None), (True, None), (True, "zlib")]
)
def test_cross_compat(mesh, binary, compression, tmp_path):
    # shim (C++ where available) write -> pure-Python read
    p = tmp_path / "cross.vtp"
    meshioplusplus.vtp.write(
        p, copy.deepcopy(mesh), binary=binary, compression=compression
    )
    out = _vtp.read(p)
    assert np.allclose(mesh.points, out.points, atol=1.0e-13, rtol=0.0)

    # pure-Python write -> shim (C++ where available) read
    _vtp.write(p, copy.deepcopy(mesh), binary=binary, compression=compression)
    out = meshioplusplus.vtp.read(p)
    assert np.allclose(mesh.points, out.points, atol=1.0e-13, rtol=0.0)


def test_lzma_roundtrip(tmp_path):
    # lzma is Python-only (the C++ path raises and the shim falls back).
    def writer(filename, mesh):
        meshioplusplus.vtp.write(filename, mesh, binary=True, compression="lzma")

    helpers.write_read(
        tmp_path,
        writer,
        meshioplusplus.vtp.read,
        helpers.tri_mesh,
        1.0e-15,
        extension=".vtp",
    )


@pytest.mark.parametrize("writer", [meshioplusplus.vtp.write, _vtp.write])
def test_volume_cells_raise(writer, tmp_path):
    with pytest.raises(meshioplusplus.WriteError):
        writer(tmp_path / "vol.vtp", copy.deepcopy(helpers.tet_mesh))


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.vtp")
    # With additional, insignificant suffix:
    helpers.generic_io(tmp_path / "test.0.vtp")
