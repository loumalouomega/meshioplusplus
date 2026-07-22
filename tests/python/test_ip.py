import numpy as np

import meshioplusplus
from meshioplusplus.ip import _ip as ip_py


def _mesh():
    pts = np.array([[0.0, 0], [1, 0], [0, 1], [1, 1]])  # 2D
    return meshioplusplus.Mesh(
        pts,
        [],
        point_data={
            "x-velocity": np.array([1.0, 2.0, 3.0, 4.0]),
            "pressure": np.array([10.0, 20.0, 30.0, 40.0]),
        },
    )


def test_roundtrip(tmp_path):
    """IP round-trips node coordinates and every field component."""
    mesh = _mesh()
    p = tmp_path / "a.ip"
    meshioplusplus.ip.write(p, mesh)
    out = meshioplusplus.ip.read(p)
    assert np.allclose(out.points, mesh.points)
    assert np.allclose(out.point_data["x-velocity"], [1, 2, 3, 4])
    assert np.allclose(out.point_data["pressure"], [10, 20, 30, 40])


def test_cpp_python_parity(tmp_path):
    mesh = _mesh()
    p_cpp = str(tmp_path / "cpp.ip")
    p_py = str(tmp_path / "py.ip")
    meshioplusplus.ip.write(p_cpp, mesh)  # C++ path
    ip_py.write(p_py, mesh)  # Python reference
    for reader, path in ((meshioplusplus.ip.read, p_py), (ip_py.read, p_cpp)):
        out = reader(path)
        assert np.allclose(out.points, mesh.points)
        assert np.allclose(out.point_data["pressure"], [10, 20, 30, 40])


def test_read_version2(tmp_path):
    """A hand-written version-2 IP file (no parentheses) reads correctly."""
    text = "2\n2\n3\n1\ntemp\n1.0\n2.0\n3.0\n4.0\n5.0\n6.0\n100.0\n200.0\n300.0\n"
    p = tmp_path / "v2.ip"
    p.write_text(text)
    out = meshioplusplus.ip.read(p)
    assert np.allclose(out.points, [[1, 4], [2, 5], [3, 6]])
    assert np.allclose(out.point_data["temp"], [100, 200, 300])


def test_generic_io(tmp_path):
    mesh = _mesh()
    p = tmp_path / "b.ip"
    meshioplusplus.write(p, mesh)
    out = meshioplusplus.read(p)
    assert np.allclose(out.point_data["x-velocity"], [1, 2, 3, 4])
