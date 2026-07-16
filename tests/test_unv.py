import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core

from . import helpers

# wedge15 (Solid Parabolic Wedge, UNV descriptor 113) is not one of the shared
# helper fixtures, so build one locally.
wedge15_mesh = meshioplusplus.Mesh(
    np.arange(45.0).reshape(15, 3),
    [("wedge15", np.arange(15).reshape(1, 15))],
)


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.line_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.triangle6_mesh,
        helpers.quad_mesh,
        helpers.quad8_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
        helpers.wedge_mesh,
        wedge15_mesh,
    ],
)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.unv.write, meshioplusplus.unv.read, mesh, 1.0e-12
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.unv")
    helpers.generic_io(tmp_path / "test.0.unv")


def test_groups(tmp_path):
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    mesh.point_sets = {"corners": np.array([0, 2])}
    mesh.cell_sets = {"all": [np.array([0, 1])]}
    p = tmp_path / "g.unv"
    meshioplusplus.unv.write(p, mesh)
    out = meshioplusplus.unv.read(p)
    assert np.array_equal(out.point_sets["corners"], [0, 2])
    assert np.array_equal(out.cell_sets["all"][0], [0, 1])


def _field_mesh():
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    mesh.point_data = {
        "temp": np.array([1.0, 2, 3, 4, 5]),
        "disp": np.arange(15.0).reshape(5, 3),
    }
    mesh.cell_data = {
        "stress": [np.arange(12.0).reshape(2, 6)],
        "unv:pid": [np.array([3, 4])],
    }
    return mesh


def _assert_fields(out):
    assert np.allclose(out.point_data["temp"], [1, 2, 3, 4, 5])
    assert np.allclose(out.point_data["disp"], np.arange(15.0).reshape(5, 3))
    assert np.allclose(out.cell_data["stress"][0], np.arange(12.0).reshape(2, 6))
    assert np.array_equal(out.cell_data["unv:pid"][0], [3, 4])


def test_fields_roundtrip(tmp_path):
    """Dataset 2414 fields (scalar/vector/tensor) at nodes and elements."""
    mesh = _field_mesh()
    p = tmp_path / "f.unv"
    meshioplusplus.unv.write(p, mesh)
    _assert_fields(meshioplusplus.unv.read(p))


def test_fields_cpp_python_parity(tmp_path):
    """C++ and Python field paths produce mutually readable output."""
    mesh = _field_mesh()
    p_cpp = str(tmp_path / "cpp.unv")
    p_py = str(tmp_path / "py.unv")
    _core.unv_write(p_cpp, mesh, {}, {})
    meshioplusplus.unv._unv.write(p_py, mesh)
    # every reader reads every writer's output
    _assert_fields(_core.unv_read(p_cpp))
    _assert_fields(_core.unv_read(p_py))
    _assert_fields(meshioplusplus.unv._unv.read(p_cpp))
    _assert_fields(meshioplusplus.unv._unv.read(p_py))


def test_fields_code_aster(tmp_path):
    """Code-Aster mode emits legacy datasets 55/57 and round-trips."""
    mesh = _field_mesh()
    p = tmp_path / "ca.unv"
    meshioplusplus.unv.write(p, mesh, code_aster=True)
    text = p.read_text()
    assert "\n    55\n" in text and "\n    57\n" in text
    _assert_fields(meshioplusplus.unv.read(p))


def test_node_dataset_781(tmp_path):
    mesh = _field_mesh()
    p = tmp_path / "n.unv"
    meshioplusplus.unv.write(p, mesh, node_dataset=781)
    assert "\n   781\n" in p.read_text()
    _assert_fields(meshioplusplus.unv.read(p))


def test_unknown_descriptor_skipped(tmp_path):
    """An unsupported FE descriptor warns and is skipped, not fatal."""
    text = (
        "    -1\n  2411\n"
        "         1         1         1        11\n   0.0   0.0   0.0\n"
        "         2         1         1        11\n   1.0   0.0   0.0\n"
        "    -1\n"
        "    -1\n  2412\n"
        "         1       999         1         1        11         2\n"
        "         1         2\n"
        "    -1\n"
    )
    p = tmp_path / "bad.unv"
    p.write_text(text)
    # neither the C++ nor the Python reader should raise
    out = meshioplusplus.unv.read(p)
    assert len(out.points) == 2 and len(out.cells) == 0
    out_py = meshioplusplus.unv._unv.read(p)
    assert len(out_py.points) == 2 and len(out_py.cells) == 0


def test_groups_cpp_path(tmp_path):
    """Force the C++ group path and confirm point_sets/cell_sets survive."""
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    mesh.point_sets = {"corners": np.array([0, 2])}
    mesh.cell_sets = {"all": [np.array([0, 1])], "first": [np.array([0])]}
    p = str(tmp_path / "g.unv")
    _core.unv_write(
        p, mesh, dict(mesh.point_sets), {k: list(v) for k, v in mesh.cell_sets.items()}
    )
    out = _core.unv_read(p)
    assert np.array_equal(out.point_sets["corners"], [0, 2])
    assert np.array_equal(out.cell_sets["all"][0], [0, 1])
    assert np.array_equal(out.cell_sets["first"][0], [0])
