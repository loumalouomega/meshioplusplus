"""Tests for point <-> cell data averaging (location conversion)."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._data_average import _cell_measures, _to_cell_py, _to_point_py

from .helpers_data import assert_same_geometry, data_mesh


def test_point_to_cell_equals_the_hand_computed_cell_means():
    m = data_mesh()
    out = mp.point_data_to_cell_data(m, keys=["T"])
    # T = {0, 1, 11, 10, 2, 12}; triangles {0,1,2} and {0,2,3}; quad {1,4,5,2}.
    assert out.cell_data["T"][0] == pytest.approx([(0 + 1 + 11) / 3, (0 + 11 + 10) / 3])
    assert out.cell_data["T"][1] == pytest.approx([(1 + 2 + 12 + 11) / 4])
    assert len(out.cell_data["T"]) == len(m.cells)
    assert_same_geometry(m, out)


def test_point_to_cell_of_a_linear_field():
    # T = x + 10*y, so each cell's mean is the field at its centroid.
    m = data_mesh()
    out = mp.point_data_to_cell_data(m, keys=["T"])
    for block, values in zip(m.cells, out.cell_data["T"]):
        conn = np.asarray(block.data)
        centroids = m.points[conn].mean(axis=1)
        expected = centroids[:, 0] + 10 * centroids[:, 1]
        assert values == pytest.approx(expected)


@pytest.mark.parametrize(
    ("cell_type", "points", "expected"),
    [
        ("tetra", [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], 1.0 / 6.0),
        (
            "hexahedron",
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
            1.0,
        ),
        (
            "wedge",
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]],
            0.5,
        ),
        (
            "pyramid",
            [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 1]],
            1.0 / 3.0,
        ),
    ],
)
def test_cell_measures_are_translation_invariant(cell_type, points, expected):
    """`_cell_measures`' divergence-theorem volume decomposes a closed surface
    into tetrahedra from the ORIGIN, so it gives the true (translation-
    invariant) volume only when every face table entry is wound outward the
    same way -- one inward face silently makes the "volume" grow with
    distance from the origin instead of depending only on the cell's shape.

    All four `_FACES` tables (tetra/hexahedron/wedge/pyramid) had exactly
    this defect until fixed to match `detail/cell_faces.cpp`'s
    Newell-normal-gtested windings verbatim: the same reference-element cell
    translated far from the origin used to report a wildly different
    "volume" purely from where it sat in space.
    """
    p0 = np.asarray(points, dtype=float)
    p1 = p0 + np.array([37.0, -19.0, 53.0])
    pts = np.vstack([p0, p1])
    conn = np.array(
        [list(range(len(p0))), list(range(len(p0), 2 * len(p0)))], dtype=np.int64
    )
    mesh = mp.Mesh(pts, [(cell_type, conn)])
    measures = _cell_measures(mesh, mesh.cells[0])
    assert measures == pytest.approx([expected, expected])


def test_cell_to_point_unweighted_is_the_incident_cell_mean():
    m = data_mesh()
    out = mp.cell_data_to_point_data(m, keys=["mat"])
    # mat: triangles = {1, 2}, quad = {3}
    got = out.point_data["mat"]
    assert got[0] == pytest.approx(1.5)  # triangles 0 and 1
    assert got[2] == pytest.approx(2.0)  # both triangles and the quad
    assert got[3] == pytest.approx(2.0)  # triangle 1 only
    assert got[4] == pytest.approx(3.0)  # quad only
    assert_same_geometry(m, out)


def test_weighted_differs_from_unweighted_on_a_non_uniform_mesh():
    # Each triangle has area 1/2; the quad has area 1.
    m = data_mesh()
    plain = mp.cell_data_to_point_data(m, keys=["mat"])
    weighted = mp.cell_data_to_point_data(m, keys=["mat"], weighted=True)
    assert plain.point_data["mat"][2] == pytest.approx((1 + 2 + 3) / 3)
    # (0.5*1 + 0.5*2 + 1*3) / 2 = 2.25
    assert weighted.point_data["mat"][2] == pytest.approx(2.25)
    assert plain.point_data["mat"][2] != weighted.point_data["mat"][2]


def test_constant_field_survives_the_round_trip():
    m = data_mesh()
    m.point_data["T"] = np.full(6, 7.0)
    mid = mp.point_data_to_cell_data(m, keys=["T"])
    back = mp.cell_data_to_point_data(mid, keys=["T"])
    assert back.point_data["T"] == pytest.approx(np.full(6, 7.0))


def test_vector_arrays_are_handled_component_wise():
    m = data_mesh()
    out = mp.point_data_to_cell_data(m, keys=["v"])
    assert out.cell_data["v"][0].shape == (2, 3)
    # triangle {0,1,2}: (1,0,0)+(0,1,0)+(0,0,1) -> (1/3, 1/3, 1/3)
    assert out.cell_data["v"][0][0] == pytest.approx([1 / 3, 1 / 3, 1 / 3])


def test_integer_input_yields_float_output():
    m = data_mesh()
    out = mp.cell_data_to_point_data(m, keys=["tag"])
    assert np.issubdtype(out.point_data["tag"].dtype, np.floating)


def test_isolated_point_yields_nan():
    m = mp.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [5, 5, 0]]),
        [("triangle", np.array([[0, 1, 2]]))],
    )
    m.cell_data["mat"] = [np.array([4.0])]
    out = mp.cell_data_to_point_data(m)
    assert out.point_data["mat"][0] == pytest.approx(4.0)
    assert np.isnan(out.point_data["mat"][3])


def test_suffix_leaves_the_source_alone():
    m = data_mesh()
    out = mp.point_data_to_cell_data(m, keys=["T"], suffix="_c")
    assert "T_c" in out.cell_data
    assert "T" in out.point_data
    assert "T" not in out.cell_data


def test_unknown_name_raises():
    m = data_mesh()
    with pytest.raises(ValueError):
        mp.point_data_to_cell_data(m, keys=["nope"])


def test_geometry_unchanged_by_both_directions():
    m = data_mesh()
    assert_same_geometry(m, mp.point_data_to_cell_data(m))
    assert_same_geometry(m, mp.cell_data_to_point_data(m))


def test_cpp_matches_python():
    pytest.importorskip("meshioplusplus._core")
    m = data_mesh()

    cpp = mp.point_data_to_cell_data(m, keys=["T", "v"])
    py = _to_cell_py(m, ["T", "v"], "", "", True, "ignore", 0.0)
    for name in ("T", "v"):
        for a, b in zip(cpp.cell_data[name], py.cell_data[name]):
            assert np.allclose(a, b, equal_nan=True)

    for weighted in (False, True):
        cpp = mp.cell_data_to_point_data(m, keys=["mat"], weighted=weighted)
        py = _to_point_py(
            m,
            ["mat"],
            "measure" if weighted else "uniform",
            "",
            "",
            True,
            "ignore",
            0.0,
        )
        assert np.allclose(cpp.point_data["mat"], py.point_data["mat"], equal_nan=True)


def test_roundtrip_write_read(tmp_path):
    m = data_mesh()
    out = mp.point_data_to_cell_data(m, keys=["T"], suffix="_c")
    path = tmp_path / "out.vtu"
    mp.write(path, out)
    back = mp.read(path)
    assert np.allclose(
        np.concatenate(back.cell_data["T_c"]), np.concatenate(out.cell_data["T_c"])
    )
