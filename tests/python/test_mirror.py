"""mirror (roadmap 3.2.5) through the Python API.

A reflecting ``transform`` leaves every cell inverted -- the probe that made
``mirror`` an operation rather than documentation -- and ``mirror`` keeps
them right side out. The flip tables are pinned in tests/cpp/test_mirror.cpp.
"""

import numpy as np
import pytest

import meshioplusplus as mio

_core = pytest.importorskip("meshioplusplus._core")
if not hasattr(_core, "mirror"):  # pragma: no cover - an older core
    pytest.skip("core predates mirror", allow_module_level=True)


def test_a_reflecting_transform_inverts_and_mirror_does_not():
    box = mio.box(dims=(2, 2, 2), cell_type="tetra")
    flipped = mio.transform(box, matrix=np.diag([-1.0, 1.0, 1.0, 1.0]))
    assert mio.check_quality(flipped)["num_inverted"] == 48
    out, report = mio.mirror(box, return_report=True)
    assert mio.check_quality(out)["num_inverted"] == 0
    assert report["num_flipped_cells"] == 48
    np.testing.assert_allclose(out.points[:, 0], -box.points[:, 0])


@pytest.mark.parametrize("cell_type", ["tetra", "hexahedron", "wedge"])
@pytest.mark.parametrize("order", [1, 2])
def test_every_box_comes_out_right_side_out(cell_type, order):
    box = mio.box(dims=(2, 3, 2), cell_type=cell_type, order=order)
    out = mio.mirror(box, origin=(0.1, 0.2, 0.3), normal=(1.0, 2.0, -0.5))
    assert mio.check_quality(out)["num_inverted"] == 0


def test_keep_original_and_weld_make_a_whole_model():
    half = mio.box(dims=(2, 2, 2), cell_type="tetra")
    out, report = mio.mirror(half, keep_original=True, return_report=True)
    assert report["num_points_on_plane"] == 9
    assert len(out.points) == 2 * len(half.points) - 9
    assert report["num_dropped_sides"] == 8 and report["warnings"]
    _, conf = mio.check_conformity(out, return_report=True)
    assert conf["conforming"]
    xmin = [r for r in out.regions if r.name == "xmin"]
    assert len(xmin) == 1 and len(xmin[0].entries) == 0
    unwelded = mio.mirror(half, keep_original=True, weld=False)
    assert len(unwelded.points) == 2 * len(half.points)


def test_vector_data_is_reflected_on_the_copies():
    m = mio.box(dims=(1, 1, 1), origin=(1.0, 0.0, 0.0))
    m.point_data["v"] = np.tile([1.0, 2.0, 3.0], (len(m.points), 1))
    out = mio.mirror(m, keep_original=True, reflect_vector_data=True)
    n = len(m.points)
    np.testing.assert_array_equal(out.point_data["v"][:n], m.point_data["v"])
    np.testing.assert_array_equal(out.point_data["v"][n:, 0], -1.0)


def test_refusals():
    box = mio.box()
    with pytest.raises(ValueError, match="normal"):
        mio.mirror(box, normal=(0, 0, 0))
    with pytest.raises(ValueError, match="normal"):
        mio.mirror(box, normal=(1, 0))
    planar = mio.Mesh(
        np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]), [("triangle", [[0, 1, 2]])]
    )
    with pytest.raises(ValueError, match="component 2"):
        mio.mirror(planar, normal=(0, 0, 1))
    assert mio.mirror(planar, normal=(0, 1, 0)).points.shape == (3, 2)
