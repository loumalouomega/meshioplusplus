"""Tests for the affine transform operation."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import transform
from meshioplusplus._transform import _build_matrix, _transform_py


def _cube():
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )
    return meshioplusplus.Mesh(
        pts, [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))]
    )


def test_translate_moves_points_only():
    mesh = _cube()
    out = transform(mesh, translate=[10, 20, 30])
    assert np.allclose(out.points, mesh.points + [10, 20, 30])
    # connectivity unchanged
    assert np.array_equal(out.cells[0].data, mesh.cells[0].data)


def test_scale_per_axis_and_uniform():
    mesh = _cube()
    out = transform(mesh, scale=[2, 3, 4])
    assert np.allclose(out.points[6], [2, 3, 4])
    out2 = transform(mesh, scale=5)
    assert np.allclose(out2.points[6], [5, 5, 5])


def test_rotate_90_about_z_maps_known_point():
    mesh = _cube()
    out = transform(mesh, rotate=("z", 90))
    # (1,0,0) -> (0,1,0)
    assert np.allclose(out.points[1], [0, 1, 0], atol=1e-12)
    # (0,1,0) -> (-1,0,0)
    assert np.allclose(out.points[3], [-1, 0, 0], atol=1e-12)


def test_scale_units_multiplies_coordinates():
    mesh = _cube()
    out = transform(mesh, scale_units=0.001)
    assert np.allclose(out.points, mesh.points * 0.001)


def test_matrix_equivalent_to_translate():
    mesh = _cube()
    mat = np.eye(4)
    mat[:3, 3] = [1, 2, 3]
    out = transform(mesh, matrix=mat)
    assert np.allclose(out.points, mesh.points + [1, 2, 3])


def test_rotate_vector_point_data():
    mesh = _cube()
    mesh.point_data["v"] = np.tile([1.0, 0.0, 0.0], (8, 1))
    out = transform(mesh, rotate=("z", 90), rotate_vector_data=True)
    # each vector (1,0,0) rotates to (0,1,0)
    assert np.allclose(
        out.point_data["v"], np.tile([0.0, 1.0, 0.0], (8, 1)), atol=1e-12
    )
    # without the flag, vector data is untouched
    out2 = transform(mesh, rotate=("z", 90))
    assert np.allclose(out2.point_data["v"], np.tile([1.0, 0.0, 0.0], (8, 1)))


def test_rotate_vector_cell_data():
    # A vector living on a CELL rotates exactly as one living on a point does.
    # Before v10.33.0 the flag reached point data only and cell data rode
    # through in the old frame -- silently, since nothing about an array says
    # which frame it is in.
    mesh = _cube()
    mesh.cell_data["vel"] = [np.array([[1.0, 0.0, 0.0]])]
    mesh.cell_data["stress"] = [np.diag([1.0, 2.0, 3.0]).reshape(1, 9)]
    mesh.cell_data["mat"] = [np.array([7], dtype=np.int64)]
    out = transform(mesh, rotate=("z", 90), rotate_vector_data=True)
    assert np.allclose(out.cell_data["vel"][0][0], [0.0, 1.0, 0.0], atol=1e-12)
    assert np.allclose(
        out.cell_data["stress"][0][0].reshape(3, 3),
        np.diag([2.0, 1.0, 3.0]),
        atol=1e-12,
    )
    # an integer tag is never a vector
    assert np.array_equal(out.cell_data["mat"][0], [7])
    # and without the flag cell data is untouched
    plain = transform(mesh, rotate=("z", 90))
    assert np.allclose(plain.cell_data["vel"][0][0], [1.0, 0.0, 0.0])


def test_cpp_matches_python():
    core = pytest.importorskip("meshioplusplus._core")
    mesh = _cube()
    mesh.point_data["v"] = np.tile([1.0, 2.0, 3.0], (8, 1))
    mesh.cell_data["cv"] = [np.array([[1.0, 2.0, 3.0]])]
    mesh.cell_data["ct"] = [np.arange(9, dtype=np.float64).reshape(1, 9)]
    mat = _build_matrix(None, None, ("z", 37.0), None, None)
    got = core.transform(mesh, mat.reshape(-1).tolist(), True)
    ref = _transform_py(mesh, mat, True)
    assert np.allclose(got.points, ref.points, atol=1e-12)
    assert np.allclose(got.point_data["v"], ref.point_data["v"], atol=1e-12)
    # The two engines must agree on cell data as well, or a Windows CI run
    # (which takes the numpy path) would disagree with a Linux one.
    for name in ("cv", "ct"):
        assert np.allclose(got.cell_data[name][0], ref.cell_data[name][0], atol=1e-12)


def test_roundtrip_write_read(tmp_path):
    mesh = _cube()
    out = transform(mesh, translate=[1, 2, 3])
    p = tmp_path / "t.vtu"
    meshioplusplus.write(p, out)
    back = meshioplusplus.read(p)
    assert np.allclose(back.points, out.points)


def test_sets_pass_through():
    mesh = _cube()
    mesh.point_sets = {"a": np.array([0, 1, 2])}
    out = transform(mesh, translate=[1, 0, 0])
    assert np.array_equal(out.point_sets["a"], [0, 1, 2])


def _ragged():
    """A ragged polygon block (a 3-gon and a 4-gon) and a polyhedron block."""
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [2, 0, 0],
            [0.5, 0.5, 1],
        ],
        dtype=float,
    )
    polygons = [np.array([0, 1, 4]), np.array([0, 1, 2, 3])]
    pyramid = [
        [
            np.array([0, 3, 2, 1]),
            np.array([0, 1, 5]),
            np.array([1, 2, 5]),
            np.array([2, 3, 5]),
            np.array([3, 0, 5]),
        ]
    ]
    return meshioplusplus.Mesh(
        pts,
        [("polygon", polygons), ("polyhedron5", pyramid)],
        point_data={"v": np.tile([1.0, 0.0, 0.0], (6, 1))},
        cell_data={"c": [np.array([1.0, 2.0]), np.array([3.0])]},
    )


def _assert_same_cells(a, b):
    assert [cb.type for cb in a.cells] == [cb.type for cb in b.cells]
    for ca, cb in zip(a.cells, b.cells):
        assert len(ca.data) == len(cb.data)
        for ra, rb in zip(ca.data, cb.data):
            if isinstance(ra, list) or isinstance(rb, list):
                assert len(ra) == len(rb)
                for fa, fb in zip(ra, rb):
                    assert np.array_equal(np.asarray(fa), np.asarray(fb))
            else:
                assert np.array_equal(np.asarray(ra), np.asarray(rb))


@pytest.mark.parametrize("rotate_vector_data", [False, True])
def test_ragged_mesh_core_matches_twin(rotate_vector_data):
    # transform used to fail on every ragged mesh: the binding refused it and
    # the twin's np.array(list) raised on the inhomogeneous rows.
    mesh = _ragged()
    m = _build_matrix([1, 2, 3], None, ("z", 90), None, None)
    out = transform(mesh, matrix=m, rotate_vector_data=rotate_vector_data)
    ref = _transform_py(mesh, m, rotate_vector_data)
    assert np.allclose(out.points, ref.points)
    _assert_same_cells(out, mesh)
    _assert_same_cells(ref, mesh)
    assert np.allclose(out.point_data["v"], ref.point_data["v"])
    for a, b in zip(out.cell_data["c"], ref.cell_data["c"]):
        assert np.array_equal(a, b)


def test_ragged_twin_does_not_alias_input():
    mesh = _ragged()
    out = _transform_py(mesh, np.eye(4), False)
    out.cells[0].data[1][0] = 99
    out.cells[1].data[0][0][0] = 99
    assert mesh.cells[0].data[1][0] == 0
    assert mesh.cells[1].data[0][0][0] == 0


def test_ragged_mesh_runs_in_core():
    pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._fallback import set_strict_core

    set_strict_core(True)  # a decline is now an error, not a twin run
    try:
        out = transform(_ragged(), translate=[0, 0, 1])
    finally:
        set_strict_core(None)
    assert np.allclose(out.points[:, 2], _ragged().points[:, 2] + 1)
