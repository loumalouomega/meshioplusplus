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
