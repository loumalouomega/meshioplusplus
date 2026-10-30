"""find_intersections (roadmap 3.1.2): self-intersection and interpenetration."""

import numpy as np
import pytest

import meshioplusplus as mio

_core = pytest.importorskip("meshioplusplus._core")
if not hasattr(_core, "find_intersections"):
    pytest.skip("core built without find_intersections", allow_module_level=True)


def _sphere(cx=0.0, n=3):
    return mio.icosa_sphere(1.0, n, center=(cx, 0.0, 0.0))


def _counts(mesh):
    return np.concatenate(mesh.cell_data["intersections:count"])


def test_a_clean_sphere_has_no_pairs():
    out, rep = mio.find_intersections(_sphere(), return_report=True)
    assert rep["clean"]
    assert rep["num_pairs"] == 0
    assert rep["pairs"].shape == (0, 2)
    assert rep["other"] is None
    assert not _counts(out).any()
    # Bit-identical but for the flags.
    np.testing.assert_array_equal(out.points, _sphere().points)


def test_two_overlapping_spheres():
    a, b = _sphere(), _sphere(1.0)
    flagged, rep = mio.find_intersections(a, b, return_report=True)
    assert not rep["clean"]
    assert rep["num_pairs"] == len(rep["pairs"]) > 0
    assert rep["pairs"].dtype == np.int64
    assert _counts(flagged).sum() == rep["num_pairs"]
    assert _counts(rep["other"]).sum() == rep["num_pairs"]
    seg = rep["segments"]
    assert len(seg.cells[0].data) == rep["num_segments"] > 0
    # On the circle x = 0.5, r = sqrt(0.75).
    np.testing.assert_allclose(seg.points[:, 0], 0.5, atol=0.05)
    np.testing.assert_allclose(
        np.hypot(seg.points[:, 1], seg.points[:, 2]), np.sqrt(0.75), atol=0.05
    )
    # The merged surface finds the same pairs on its own.
    merged = mio.merge([a, b])
    _, self_rep = mio.find_intersections(merged, return_report=True)
    assert self_rep["num_pairs"] == rep["num_pairs"]
    na = len(a.cells[0].data)
    np.testing.assert_array_equal(self_rep["pairs"][:, 0], rep["pairs"][:, 0])
    np.testing.assert_array_equal(self_rep["pairs"][:, 1], rep["pairs"][:, 1] + na)


def test_cap_options_and_errors():
    a, b = _sphere(), _sphere(1.0)
    out, rep = mio.find_intersections(
        a, b, max_pairs=2, segments=False, attach_flags=False, return_report=True
    )
    assert rep["truncated"] and len(rep["pairs"]) == 2
    assert rep["segments"] is None
    assert "intersections:count" not in out.cell_data
    assert rep["warnings"]
    with pytest.raises(ValueError):
        mio.find_intersections(a, max_pairs=-1)
    with pytest.raises(ValueError):
        mio.find_intersections(mio.box())
    with pytest.raises(ValueError):
        mio.find_intersections(a, region="nope")


def test_disjoint_parts_and_a_fold():
    assert mio.find_intersections(_sphere(), _sphere(3.0), return_report=True)[1][
        "clean"
    ]
    fold = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0.5, 0.25, 0]], dtype=float),
        [("triangle", np.array([[0, 1, 2], [0, 1, 3]]))],
    )
    _, rep = mio.find_intersections(fold, return_report=True)
    assert rep["num_pairs"] == 1 and rep["num_coplanar"] == 1
