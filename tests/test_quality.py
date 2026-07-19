import numpy as np
import pytest

import meshioplusplus
from meshioplusplus._quality import _compute_quality_py

from . import helpers

SQRT2 = np.sqrt(2.0)


def _val(report, name, block=0, cell=0):
    return report["cell_arrays"][name][block].ravel()[cell]


def test_unit_cube_hex_is_ideal():
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
    mesh = meshioplusplus.Mesh(
        pts, [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))]
    )
    r = meshioplusplus.compute_quality(mesh)
    assert _val(r, "quality:volume") == pytest.approx(1.0)
    assert _val(r, "quality:scaled_jacobian") == pytest.approx(1.0)
    assert _val(r, "quality:skewness") == pytest.approx(0.0, abs=1e-12)
    assert r["num_inverted"] == 0 and r["num_degenerate"] == 0


def test_regular_tetra_is_ideal():
    s = SQRT2
    pts = np.array(
        [
            [0, 0, 0],
            [s, 0, 0],
            [s / 2, s * np.sqrt(3) / 2, 0],
            [s / 2, s * np.sqrt(3) / 6, s * np.sqrt(2 / 3)],
        ],
        dtype=float,
    )
    mesh = meshioplusplus.Mesh(pts, [("tetra", np.array([[0, 1, 2, 3]]))])
    r = meshioplusplus.compute_quality(mesh)
    assert _val(r, "quality:scaled_jacobian") == pytest.approx(1.0)
    assert _val(r, "quality:aspect_ratio") == pytest.approx(1.0)
    ideal = np.degrees(np.arccos(1 / 3))
    assert _val(r, "quality:min_dihedral") == pytest.approx(ideal)
    assert _val(r, "quality:max_dihedral") == pytest.approx(ideal)


def test_equilateral_triangle_2d():
    pts = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]], dtype=float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2]]))])
    r = meshioplusplus.compute_quality(mesh)
    assert _val(r, "quality:aspect_ratio") == pytest.approx(1.0)
    assert _val(r, "quality:min_angle") == pytest.approx(60.0)
    assert _val(r, "quality:max_angle") == pytest.approx(60.0)


def test_unit_square_2d():
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
    mesh = meshioplusplus.Mesh(pts, [("quad", np.array([[0, 1, 2, 3]]))])
    r = meshioplusplus.compute_quality(mesh)
    assert _val(r, "quality:volume") == pytest.approx(1.0)
    assert _val(r, "quality:skewness") == pytest.approx(0.0, abs=1e-12)
    assert _val(r, "quality:aspect_ratio") == pytest.approx(1.0)
    assert _val(r, "quality:warpage") == pytest.approx(0.0, abs=1e-12)


def test_inverted_tetra_flagged():
    s = SQRT2
    pts = np.array(
        [
            [0, 0, 0],
            [s, 0, 0],
            [s / 2, s * np.sqrt(3) / 2, 0],
            [s / 2, s * np.sqrt(3) / 6, s * np.sqrt(2 / 3)],
        ],
        dtype=float,
    )
    # swap two nodes -> negative signed volume
    mesh = meshioplusplus.Mesh(pts, [("tetra", np.array([[0, 2, 1, 3]]))])
    r = meshioplusplus.compute_quality(mesh)
    assert _val(r, "quality:volume") < 0.0
    assert _val(r, "quality:inverted") == 1.0
    assert r["num_inverted"] == 1


def test_degenerate_triangle_flagged():
    pts = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2]]))])
    r = meshioplusplus.compute_quality(mesh)
    assert _val(r, "quality:degenerate") == 1.0
    assert r["num_degenerate"] == 1
    assert np.isnan(_val(r, "quality:aspect_ratio"))


def test_attach_quality_adds_cell_data():
    mesh = helpers.tet_mesh.copy()
    out = meshioplusplus.attach_quality(mesh)
    assert "quality:volume" in out.cell_data
    assert "quality:inverted" in out.cell_data
    assert len(out.cell_data["quality:volume"]) == len(mesh.cells)


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tri_mesh,
        helpers.quad_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
        helpers.wedge_mesh,
        helpers.pyramid_mesh,
        helpers.tri_quad_mesh,
    ],
)
def test_cpp_matches_python(mesh):
    from meshioplusplus import _core

    a = _core.compute_quality(mesh)
    b = _compute_quality_py(mesh)
    assert a["num_cells"] == b["num_cells"]
    assert a["num_inverted"] == b["num_inverted"]
    assert a["num_degenerate"] == b["num_degenerate"]
    for name in b["cell_arrays"]:
        for x, y in zip(a["cell_arrays"][name], b["cell_arrays"][name]):
            np.testing.assert_allclose(x, y, rtol=1e-9, atol=1e-10, equal_nan=True)
