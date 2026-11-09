"""thickness (roadmap 3.1.4): wall thickness and clearance by a ray per facet."""

import math

import numpy as np
import pytest

import meshioplusplus as mio

_core = pytest.importorskip("meshioplusplus._core")
if not hasattr(_core, "thickness"):
    pytest.skip("core built without thickness", allow_module_level=True)


def _dist(mesh):
    return np.concatenate(mesh.cell_data["thickness:distance"])


def test_a_sphere_is_as_thick_as_its_diameter():
    out, rep = mio.thickness(mio.icosa_sphere(1.0, 3), return_report=True)
    assert rep["direction"] == "inward"
    assert rep["num_missed"] == 0
    assert rep["min"] == pytest.approx(2.0, abs=0.02)
    d = _dist(out)
    assert np.all(d <= 2.0) and np.all(d > 1.95)
    hit = np.concatenate(out.cell_data["thickness:hit_cell"])
    assert hit.dtype == np.int64 and np.all(hit >= 0)


def test_clearance_to_a_second_part_and_threshold():
    a = mio.icosa_sphere(1.0, 3)
    b = mio.icosa_sphere(1.0, 3, center=(3.0, 0.0, 0.0))
    out, rep = mio.thickness(a, b, threshold=1.2, return_report=True)
    assert rep["direction"] == "outward"
    assert 1.0 <= rep["min"] < 1.05
    assert rep["min_point"][0] > 0.9
    assert 0 < rep["num_below"] < rep["num_hit"]
    assert rep["num_missed"] > 0
    assert np.isnan(_dist(out)).sum() == rep["num_missed"]
    # Bounded rays miss.
    _, short = mio.thickness(a, b, max_distance=0.5, return_report=True)
    assert short["num_hit"] == 0


def test_open_sheet_and_errors():
    out, rep = mio.thickness(mio.rectangle((2, 2)), return_report=True)
    assert rep["num_hit"] == 0 and rep["num_missed"] == 4
    assert all(math.isnan(v) for v in _dist(out))
    with pytest.raises(ValueError):
        mio.thickness(mio.rectangle(), direction="sideways")
    with pytest.raises(ValueError):
        mio.thickness(mio.rectangle(), max_distance=-1.0)
    with pytest.raises(ValueError):
        mio.thickness(mio.box())
    out = mio.thickness(mio.icosa_sphere(1.0, 1), attach=False)
    assert "thickness:distance" not in out.cell_data


def _wedge(L=4.0, W=8.0, angle=np.radians(30.0), k=40):
    """A closed wedge: a fine bottom (z = 0) under a top z = x tan(angle)."""
    pts, tris = [], []
    for i in range(k + 1):
        for j in range(k + 1):
            pts.append((L * i / k, W * j / k, 0.0))
    for i in range(k):
        for j in range(k):
            a, b = i * (k + 1) + j, (i + 1) * (k + 1) + j
            tris += [(a, b + 1, b), (a, a + 1, b + 1)]
    h = L * np.tan(angle)
    base = len(pts)
    pts += [(0, 0, 0), (L, 0, 0), (L, 0, h), (0, W, 0), (L, W, 0), (L, W, h)]
    a0, b0, c0, a1, b1, c1 = range(base, base + 6)
    tris += [
        (a0, c1, a1),
        (a0, c0, c1),
        (b0, b1, c1),
        (b0, c1, c0),
        (a0, b0, c0),
        (a1, c1, b1),
    ]
    return mio.Mesh(np.array(pts, float), [("triangle", np.array(tris))])


def test_the_sphere_and_cone_methods_on_a_wedge():
    wedge = _wedge()
    angle = np.radians(30.0)
    cells = wedge.cells[0].data[: 2 * 40 * 40]
    centroids = wedge.points[cells].mean(axis=1)
    f = int(np.argmin(np.hypot(centroids[:, 0] - 2.0, centroids[:, 1] - 4.0)))
    x = centroids[f, 0]
    ray = _dist(mio.thickness(wedge))[f]
    assert ray == pytest.approx(x * np.tan(angle), rel=1e-9)
    out, rep = mio.thickness(wedge, method="sphere", return_report=True)
    assert rep["method"] == "sphere"
    assert _dist(out)[f] == pytest.approx(x * np.sin(angle), rel=0.02)
    assert rep["p2"] <= rep["p95"]
    cone = _dist(mio.thickness(wedge, method="cone", cone_angle=35.0, cone_rays=128))[f]
    assert x * np.sin(angle) - 1e-9 <= cone < 1.02 * x * np.sin(angle)
    assert np.isnan(_dist(mio.thickness(wedge, max_normal_deviation=20.0))[f])
    with pytest.raises(ValueError, match="method"):
        mio.thickness(wedge, method="laser")
    two = mio.thickness(wedge, two_sided=True, return_report=True)[1]
    assert two["max_distance"] == pytest.approx(
        0.2 * np.linalg.norm([4.0, 8.0, 4 * np.tan(angle)])
    )
