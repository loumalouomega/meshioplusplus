import numpy as np

import meshioplusplus
from meshioplusplus.mff import _mff as mff_py


def _mesh():
    return meshioplusplus.Mesh(
        np.empty((5, 0)),
        [],
        point_data={"mff:field": np.array([1.5, -2.25, 3.0, 4.0, 5.0])},
    )


def test_roundtrip(tmp_path):
    """MFF preserves field values (geometry/shape are not recoverable)."""
    mesh = _mesh()
    p = tmp_path / "a.mff"
    meshioplusplus.mff.write(p, mesh)
    out = meshioplusplus.mff.read(p)
    assert np.allclose(out.point_data["mff:field"], [1.5, -2.25, 3.0, 4.0, 5.0])


def test_cpp_python_parity(tmp_path):
    mesh = _mesh()
    p_cpp = str(tmp_path / "cpp.mff")
    p_py = str(tmp_path / "py.mff")
    meshioplusplus.mff.write(p_cpp, mesh)  # C++ path
    mff_py.write(p_py, mesh)  # Python reference
    expected = [1.5, -2.25, 3.0, 4.0, 5.0]
    assert np.allclose(meshioplusplus.mff.read(p_py).point_data["mff:field"], expected)
    assert np.allclose(mff_py.read(p_cpp).point_data["mff:field"], expected)


def test_generic_io(tmp_path):
    """Extension-based dispatch round-trips a field through the top-level API."""
    mesh = _mesh()
    p = tmp_path / "b.mff"
    meshioplusplus.write(p, mesh)
    out = meshioplusplus.read(p)
    assert np.allclose(out.point_data["mff:field"], [1.5, -2.25, 3.0, 4.0, 5.0])
