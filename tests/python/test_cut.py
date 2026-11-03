"""``cut_level_set`` and ``cut_quadrature`` (doc/cut.md): the part of a mesh
that level sets select and quadrature on it, through CutCells, an optional
dependency of the C++ core. The measures are exact for planar cuts; the
refusals run everywhere, and without CutCells the operations refuse naming the
CMake option."""

import contextlib
import io
import json

import numpy as np
import pytest

import meshioplusplus as mio

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - the numpy-only install
    _core = None

needs_core = pytest.mark.skipif(
    not hasattr(_core, "cut_level_set"), reason="core predates cut_level_set"
)
needs_cutcells = pytest.mark.skipif(
    not mio.cut_has_cutcells(), reason="this build has no CutCells"
)
PLANE = [((0.3, 0.0, 0.0), (1.0, 0.0, 0.0))]


def _measure(mesh, dim):
    """The summed measure of the linear cells of dimension ``dim``."""
    total = 0.0
    for block in mesh.cells:
        pts = mesh.points[block.data]
        if dim == 3 and block.type == "tetra":
            a, b, c, d = (pts[:, k] for k in range(4))
            total += (
                np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))).sum() / 6
            )
        elif dim == 2 and block.type == "triangle":
            a, b, c = (pts[:, k] for k in range(3))
            total += np.linalg.norm(np.cross(b - a, c - a), axis=1).sum() / 2
        elif dim == 1 and block.type == "line":
            total += np.linalg.norm(pts[:, 1] - pts[:, 0], axis=1).sum()
    return total


@needs_core
def test_without_cutcells_the_operations_refuse_by_name():
    if mio.cut_has_cutcells():
        pytest.skip("CutCells is compiled in")
    box = mio.box((2, 2, 2))
    with pytest.raises(NotImplementedError, match="MESHIOPLUSPLUS_WITH_CUTCELLS"):
        mio.cut_level_set(box, planes=PLANE)
    with pytest.raises(NotImplementedError, match="MESHIOPLUSPLUS_WITH_CUTCELLS"):
        mio.cut_quadrature(box, planes=PLANE)


def test_the_python_checks_run_before_the_core():
    box = mio.box((1, 1, 1))
    with pytest.raises(ValueError, match="give a level set"):
        mio.cut_level_set(box)
    with pytest.raises(ValueError, match="backend"):
        mio.cut_level_set(box, planes=PLANE, backend="marching")
    with pytest.raises(ValueError, match="plane 0"):
        mio.cut_level_set(box, planes=[((0, 0, 0), (0, 0, 0))])
    with pytest.raises(ValueError, match="order"):
        mio.cut_quadrature(box, planes=PLANE, order=-1)


@needs_cutcells
@pytest.mark.parametrize("cell_type", ["hexahedron", "tetra", "wedge"])
def test_a_plane_cut_is_exact_and_carries_the_data(cell_type):
    box = mio.box((3, 3, 3), cell_type=cell_type)
    box.point_data["f"] = 2 * box.points[:, 0] - box.points[:, 1] + 1
    box.cell_data["id"] = [np.arange(len(box.cells[0].data)) + 10]
    part, report = mio.cut_level_set(box, planes=PLANE, return_report=True)
    assert report["measure"] == pytest.approx(0.3, abs=1e-12)
    assert report["dimension"] == 3 and report["backend"] == "lut"
    assert report["selection"] == "plane0 < 0"
    assert np.all(part.points[:, 0] <= 0.3 + 1e-12)
    np.testing.assert_allclose(
        part.point_data["f"], 2 * part.points[:, 0] - part.points[:, 1] + 1, atol=1e-10
    )
    for k, block in enumerate(part.cells):
        parents = part.cell_data["cut:parent_cell"][k]
        np.testing.assert_array_equal(part.cell_data["id"][k], parents + 10)
    face = mio.cut_level_set(box, planes=PLANE, selection="plane0 = 0")
    assert {b.type for b in face.cells} <= {"triangle", "quad"}
    simplices = mio.cut_level_set(box, planes=PLANE, triangulate=True)
    assert {b.type for b in simplices.cells} <= {"tetra"}
    assert _measure(simplices, 3) == pytest.approx(0.3, abs=1e-12)


@needs_cutcells
def test_two_level_sets_and_the_quadrature():
    box = mio.box((2, 2, 2))
    box.point_data["a"] = box.points[:, 0] - 0.5
    box.point_data["b"] = box.points[:, 1] - 0.5
    _, report = mio.cut_level_set(
        box, ["a", "b"], selection="a < 0 and b < 0", return_report=True
    )
    assert report["measure"] == pytest.approx(0.25, abs=1e-12)
    cloud, rep = mio.cut_quadrature(
        box, ["a", "b"], selection="a = 0 and b < 0", order=3, return_report=True
    )
    assert rep["dimension"] == 2
    assert rep["sum_weights"] == pytest.approx(0.5, abs=1e-12)
    w = cloud.point_data["quadrature:weight"]
    assert np.all(np.abs(cloud.points[:, 0] - 0.5) < 1e-12)
    # A degree-3 rule integrates z^3 over the face x = 0.5, y < 0.5 exactly.
    assert float(w @ cloud.points[:, 2] ** 3) == pytest.approx(0.5 / 4, abs=1e-12)
    assert cloud.point_data["quadrature:reference"].shape == (len(w), 3)
    with pytest.raises(ValueError, match="'c'"):
        mio.cut_level_set(box, ["a"], selection="c < 0")


@needs_cutcells
def test_a_quadratic_level_set_on_tetra10_is_more_accurate():
    r, exact = 0.3, 4 / 3 * np.pi * 0.3**3

    def ball(m):
        m.point_data["ball"] = np.linalg.norm(m.points - 0.5, axis=1) - r
        return m

    lin = ball(mio.box((4, 4, 4), cell_type="tetra"))
    quad = ball(mio.box((4, 4, 4), cell_type="tetra", order=2))
    e1 = abs(
        mio.cut_quadrature(lin, "ball", return_report=True)[1]["sum_weights"] - exact
    )
    _, rep = mio.cut_quadrature(quad, "ball", order=4, return_report=True)
    assert rep["backend"] == "quadrays"
    assert abs(rep["sum_weights"] - exact) < e1 / 5


@needs_cutcells
def test_pipeline_and_cli(tmp_path):
    box = mio.box((2, 2, 2))
    out, report = mio.apply_operations(
        box, [{"Op": "CutLevelSet", "Planes": [0.3, 0, 0, 1, 0, 0]}]
    )
    step = report["steps"][0]
    assert step["Measure"] == pytest.approx(0.3, abs=1e-12)
    assert step["Dimension"] == 3
    src = tmp_path / "box.vtu"
    mio.write(src, box)
    import meshioplusplus._cli  # noqa: F401  (mio._cli.main below)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        mio._cli.main(
            [
                "cut",
                str(src),
                str(tmp_path / "c.vtu"),
                "--plane",
                "0.3,0,0,1,0,0",
                "--json",
            ]
        )
    assert json.loads(buf.getvalue())["measure"] == pytest.approx(0.3, abs=1e-12)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        mio._cli.main(
            [
                "cut-quadrature",
                str(src),
                str(tmp_path / "q.vtu"),
                "--plane",
                "0.3,0,0,1,0,0",
                "--selection",
                "plane0 = 0",
                "--json",
            ]
        )
    assert json.loads(buf.getvalue())["sum_weights"] == pytest.approx(1.0, abs=1e-12)
