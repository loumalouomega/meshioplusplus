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


def test_quad_scaled_jacobian():
    # Verdict's quad scaled Jacobian (roadmap 6.8): the least corner
    # (e1 x e2) . n / (|e1| |e2|). A square scores 1, a 60-degree rhombus
    # sin 60, a reflex (arrowhead) quad below 0, in both engines and in 3-D.
    c = np.sqrt(3.0) / 2.0
    cases = [
        ([[0, 0], [1, 0], [1, 1], [0, 1]], 1.0),
        ([[0, 0], [1, 0], [1.5, c], [0.5, c]], c),
        ([[0, 0], [2, 1], [0, 2], [0.5, 1]], None),
    ]
    for pts, want in cases:
        for points in (np.array(pts, float), np.column_stack([pts, np.zeros(4)])):
            mesh = meshioplusplus.Mesh(points, [("quad", np.array([[0, 1, 2, 3]]))])
            for r in (meshioplusplus.compute_quality(mesh), _compute_quality_py(mesh)):
                sj = _val(r, "quality:scaled_jacobian")
                if want is None:
                    assert sj < 0.0
                else:
                    assert sj == pytest.approx(want)


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


def test_cell_measures_of_reference_cells():
    s = SQRT2
    hexa = meshioplusplus.Mesh(
        np.array(
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
        ),
        [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))],
    )
    r = meshioplusplus.compute_quality(hexa)
    assert _val(r, "quality:size") == pytest.approx(1.0)
    assert _val(r, "quality:diagonal") == pytest.approx(np.sqrt(3.0))
    assert _val(r, "quality:min_edge") == pytest.approx(1.0)
    assert _val(r, "quality:max_edge") == pytest.approx(1.0)

    rect = meshioplusplus.Mesh(
        np.array([[0, 0], [2, 0], [2, 1], [0, 1]], dtype=float),
        [("quad", np.array([[0, 1, 2, 3]]))],
    )
    r = meshioplusplus.compute_quality(rect)
    assert _val(r, "quality:size") == pytest.approx(2.0)
    assert _val(r, "quality:diagonal") == pytest.approx(np.sqrt(5.0))
    assert _val(r, "quality:min_edge") == pytest.approx(1.0)
    assert _val(r, "quality:max_edge") == pytest.approx(2.0)

    tet = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    r = meshioplusplus.compute_quality(tet)
    assert _val(r, "quality:size") == pytest.approx(1.0 / 6.0)
    assert _val(r, "quality:min_edge") == pytest.approx(1.0)
    assert _val(r, "quality:max_edge") == pytest.approx(s)

    h = np.sqrt(3.0) / 2
    wedge = meshioplusplus.Mesh(
        np.array(
            [[0, 0, 0], [1, 0, 0], [0.5, h, 0], [0, 0, 1], [1, 0, 1], [0.5, h, 1]],
            dtype=float,
        ),
        [("wedge", np.array([[0, 1, 2, 3, 4, 5]]))],
    )
    r = meshioplusplus.compute_quality(wedge)
    assert _val(r, "quality:size") == pytest.approx(h / 2)
    assert _val(r, "quality:max_edge") == pytest.approx(1.0)
    assert _val(r, "quality:diagonal") == pytest.approx(s)

    pyr = meshioplusplus.Mesh(
        np.array(
            [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 1 / s]], dtype=float
        ),
        [("pyramid", np.array([[0, 1, 2, 3, 4]]))],
    )
    r = meshioplusplus.compute_quality(pyr)
    assert _val(r, "quality:size") == pytest.approx(1.0 / (3.0 * s))
    assert _val(r, "quality:min_edge") == pytest.approx(1.0)
    assert _val(r, "quality:max_edge") == pytest.approx(1.0)


def test_cell_measures_survive_a_degenerate_cell():
    pts = np.array([[0, 0, 0], [1, 0, 0], [3, 0, 0]], dtype=float)
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2]]))])
    r = meshioplusplus.compute_quality(mesh)
    assert _val(r, "quality:degenerate") == 1.0
    assert _val(r, "quality:size") == pytest.approx(0.0, abs=1e-12)
    assert _val(r, "quality:min_edge") == pytest.approx(1.0)
    assert _val(r, "quality:max_edge") == pytest.approx(3.0)
    assert _val(r, "quality:diagonal") == pytest.approx(3.0)


def test_size_is_unsigned_and_sums_to_the_stats_totals():
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
    inv = meshioplusplus.Mesh(pts, [("tetra", np.array([[0, 2, 1, 3]]))])
    r = meshioplusplus.compute_quality(inv)
    assert _val(r, "quality:volume") < 0.0
    assert _val(r, "quality:size") == pytest.approx(-_val(r, "quality:volume"))

    mesh = meshioplusplus.Mesh(pts, [("tetra", np.array([[0, 1, 2, 3]]))])
    r = meshioplusplus.compute_quality(mesh)
    total = meshioplusplus.compute_stats(mesh)["unsigned_volume"]
    assert r["cell_arrays"]["quality:size"][0].sum() == pytest.approx(total)


_S3, _S6 = np.sqrt(3.0), np.sqrt(6.0)
_REG_TRI = [[0, 0, 0], [1, 0, 0], [0.5, _S3 / 2, 0]]
_REG_TET = [[0, 0, 0], [1, 0, 0], [0.5, _S3 / 2, 0], [0.5, _S3 / 6, _S6 / 3]]
_UNIT_HEX = [
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1],
]


def _verdict(points, cell_type, n):
    mesh = meshioplusplus.Mesh(
        np.array(points, dtype=float), [(cell_type, np.array([list(range(n))]))]
    )
    r = meshioplusplus.compute_quality(mesh)
    return lambda name: _val(r, "quality:" + name)


def test_verdict_metrics_of_the_ideal_cells():
    t = _verdict(_REG_TRI, "triangle", 3)
    assert t("condition") == pytest.approx(1.0)
    assert t("shape") == pytest.approx(1.0)
    assert t("edge_ratio") == pytest.approx(1.0)
    assert np.isnan(t("stretch")) and np.isnan(t("taper"))
    q = _verdict([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], "quad", 4)
    assert q("condition") == pytest.approx(1.0)
    assert q("stretch") == pytest.approx(1.0)
    assert q("taper") == pytest.approx(0.0, abs=1e-12)
    te = _verdict(_REG_TET, "tetra", 4)
    for name in ("condition", "shape", "aspect_gamma", "edge_ratio"):
        assert te(name) == pytest.approx(1.0), name
    h = _verdict(_UNIT_HEX, "hexahedron", 8)
    for name in ("condition", "stretch", "diagonal_ratio", "edge_ratio"):
        assert h(name) == pytest.approx(1.0), name


def test_verdict_metrics_of_stretched_cells():
    # A 2 x 1 x 1 box: J = diag(2, 1, 1), |J| = sqrt 6, |J^-1| = 1.5.
    pts = [[x * 2, y, z] for x, y, z in _UNIT_HEX]
    h = _verdict(pts, "hexahedron", 8)
    assert h("edge_ratio") == pytest.approx(2.0)
    assert h("condition") == pytest.approx(np.sqrt(6.0) * 1.5 / 3.0)
    assert h("stretch") == pytest.approx(np.sqrt(3.0) / np.sqrt(6.0))
    assert h("diagonal_ratio") == pytest.approx(1.0)
    # An isosceles trapezoid: X1 = (3, 0), X2 = (0, 2), X12 = (-1, 0) -> taper 1/2.
    q = _verdict([[0, 0, 0], [2, 0, 0], [1.5, 1, 0], [0.5, 1, 0]], "quad", 4)
    assert q("taper") == pytest.approx(0.5)
    # a 2 x 1 rectangle: every corner (4 + 1) / (2 * 2) = 5/4; stretch sqrt2 / sqrt5.
    r = _verdict([[0, 0, 0], [2, 0, 0], [2, 1, 0], [0, 1, 0]], "quad", 4)
    assert r("condition") == pytest.approx(1.25)
    assert r("stretch") == pytest.approx(np.sqrt(2.0) / np.sqrt(5.0))


def test_verdict_tetra_metrics_match_an_independent_numpy_evaluation():
    rng = np.random.default_rng(5)
    for _ in range(5):
        pts = rng.normal(size=(4, 3))
        if (
            np.linalg.det(
                np.stack([pts[1] - pts[0], pts[2] - pts[0], pts[3] - pts[0]], 1)
            )
            < 0
        ):
            pts[[1, 2]] = pts[[2, 1]]
        t = _verdict(pts, "tetra", 4)
        e = np.stack([pts[1] - pts[0], pts[2] - pts[0], pts[3] - pts[0]], 1)
        winv = np.array(
            [[1, -1 / _S3, -1 / _S6], [0, 2 / _S3, -1 / _S6], [0, 0, 3 / _S6]]
        )
        f = e @ winv
        cond = np.linalg.norm(f) * np.linalg.norm(np.linalg.inv(f)) / 3
        assert t("condition") == pytest.approx(cond, rel=1e-10)
        shape = 3 * np.linalg.det(f) ** (2 / 3) / np.linalg.norm(f) ** 2
        assert t("shape") == pytest.approx(shape, rel=1e-10)
        edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        srms = np.sqrt(sum(np.sum((pts[i] - pts[j]) ** 2) for i, j in edges) / 6)
        gamma = srms**3 / (6 * np.sqrt(2) * abs(np.linalg.det(e)) / 6)
        assert t("aspect_gamma") == pytest.approx(gamma, rel=1e-10)


def test_a_triangles_condition_is_its_frobenius_aspect():
    rng = np.random.default_rng(2)
    for _ in range(4):
        p = rng.normal(size=(3, 3))
        t = _verdict(p, "triangle", 3)
        a = np.linalg.norm(p[1] - p[0])
        b = np.linalg.norm(p[2] - p[1])
        c = np.linalg.norm(p[0] - p[2])
        area = 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
        assert t("condition") == pytest.approx(
            (a * a + b * b + c * c) / (4 * _S3 * area), rel=1e-10
        )
        assert t("shape") == pytest.approx(1.0 / t("condition"))


def test_verdict_metrics_are_nan_for_a_degenerate_cell():
    t = _verdict([[0, 0, 0], [1, 0, 0], [2, 0, 0]], "triangle", 3)
    assert np.isnan(t("condition")) and np.isnan(t("shape"))
    assert t("edge_ratio") == pytest.approx(2.0)  # the lengths survive
    flat = _verdict([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], "tetra", 4)
    assert np.isnan(flat("condition")) and np.isnan(flat("aspect_gamma"))


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
