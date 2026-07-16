import numpy as np
import pytest

import meshioplusplus
from meshioplusplus.dex import _dex as dex_py


def _vector_mesh():
    pts = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
    return meshioplusplus.Mesh(
        pts, [], point_data={"mGradT": np.arange(9.0).reshape(3, 3)}
    )


def _scalar_mesh():
    pts = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
    return meshioplusplus.Mesh(
        pts, [], point_data={"temp": np.array([10.0, 20.0, 30.0])}
    )


@pytest.mark.parametrize("mesh_fn", [_vector_mesh, _scalar_mesh])
def test_roundtrip(mesh_fn, tmp_path):
    """DEX round-trips both node coordinates and field values."""
    mesh = mesh_fn()
    p = tmp_path / "a.dex"
    meshioplusplus.dex.write(p, mesh)
    out = meshioplusplus.dex.read(p)
    assert np.allclose(out.points, mesh.points)
    (name,) = mesh.point_data
    assert np.allclose(out.point_data[name], mesh.point_data[name])


def test_cpp_python_parity(tmp_path):
    mesh = _vector_mesh()
    p_cpp = str(tmp_path / "cpp.dex")
    p_py = str(tmp_path / "py.dex")
    meshioplusplus.dex.write(p_cpp, mesh)  # C++ path
    dex_py.write(p_py, mesh)  # Python reference
    vals = mesh.point_data["mGradT"]
    for reader, path in (
        (meshioplusplus.dex.read, p_py),
        (dex_py.read, p_cpp),
    ):
        out = reader(path)
        assert np.allclose(out.points, mesh.points)
        assert np.allclose(out.point_data["mGradT"], vals)


def test_generic_io(tmp_path):
    mesh = _scalar_mesh()
    p = tmp_path / "b.dex"
    meshioplusplus.write(p, mesh)
    out = meshioplusplus.read(p)
    assert np.allclose(out.point_data["temp"], [10.0, 20.0, 30.0])
