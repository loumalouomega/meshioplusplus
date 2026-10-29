"""extrude and revolve (roadmap 3.1.2, 3.1.3) through the Python API.

The C++ core does the work (tests/cpp/test_extrude.cpp pins the geometry);
these check the Python surface: arguments, the report, regions and data
carried, and the refusals by name.
"""

import math

import numpy as np
import pytest

import meshioplusplus as mio

_core = pytest.importorskip("meshioplusplus._core")
if not hasattr(_core, "extrude"):  # pragma: no cover - an older core
    pytest.skip("core predates extrude", allow_module_level=True)


def _sides(mesh, name):
    for r in mesh.regions:
        if r.name == name and r.kind == "side":
            return {tuple(e) for e in np.asarray(r.entries).tolist()}
    raise KeyError(name)


def test_extruding_a_rectangle_gives_the_box():
    rect = mio.rectangle(dims=(3, 2))
    out, report = mio.extrude(rect, (0.0, 0.0, 1.0), num_layers=4, return_report=True)
    box = mio.box(dims=(3, 2, 4))
    np.testing.assert_array_equal(out.points, box.points)
    np.testing.assert_array_equal(out.cells[0].data, box.cells[0].data)
    assert report["num_hexahedra"] == 24 and report["num_flipped"] == 0
    assert report["min_scaled_jacobian"] == 1.0
    assert "num_on_axis" not in report and "full" not in report
    assert _sides(out, "bottom") == _sides(box, "zmin")
    assert _sides(out, "top") == _sides(box, "zmax")
    for name in ("xmin", "xmax", "ymin", "ymax"):
        assert _sides(out, name) == _sides(box, name)
    assert _sides(out, "side") == set().union(
        *(_sides(box, n) for n in ("xmin", "xmax", "ymin", "ymax"))
    )


def test_grading_offsets_and_layer_data():
    tri = mio.rectangle(dims=(2, 2), cell_type="triangle")
    tri.cell_data["id"] = [np.arange(8, dtype=np.int32)]
    out, report = mio.extrude(
        tri, (0, 0, 2), num_layers=3, ratio=2.0, record_layers=True, return_report=True
    )
    assert report["num_wedges"] == 24
    np.testing.assert_allclose(
        np.unique(out.points[:, 2]), [0.0, 2 / 7, 6 / 7, 2.0], atol=1e-15
    )
    np.testing.assert_array_equal(out.cell_data["id"][0], np.tile(np.arange(8), 3))
    np.testing.assert_array_equal(
        out.cell_data["extrude:layer"][0], np.repeat([1, 2, 3], 8)
    )
    by_offsets = mio.extrude(tri, (0, 0, 2), offsets=[1 / 7, 3 / 7, 1.0])
    np.testing.assert_allclose(by_offsets.points, out.points, atol=1e-15)


def test_a_curve_extrudes_to_a_surface():
    t = np.linspace(0.0, 2 * math.pi, 13)[:-1]
    pts = np.column_stack([np.cos(t), np.sin(t), np.zeros_like(t)])
    lines = np.array([[i, (i + 1) % 12] for i in range(12)])
    out, report = mio.extrude(
        mio.Mesh(pts, [("line", lines)]), num_layers=2, return_report=True
    )
    assert report["num_quads"] == 24 and report["min_scaled_jacobian"] == 0.0
    assert [c.type for c in out.cells] == ["quad"]
    assert len(_sides(out, "bottom")) == 12 and len(_sides(out, "top")) == 12


def test_revolve_a_rectangle_about_its_edge():
    out, report = mio.revolve(
        mio.rectangle(dims=(4, 3), size=(2.0, 1.5)), num_layers=12, return_report=True
    )
    assert report["full"] and report["num_on_axis"] == 4
    assert report["num_wedges"] == 36 and report["num_hexahedra"] == 108
    assert mio.check_quality(out)["num_inverted"] == 0
    _, conf = mio.check_conformity(out, return_report=True)
    assert conf["conforming"] and conf["num_boundary_facets"] == len(
        _sides(out, "side")
    )
    assert not any(r.name in ("start", "end") for r in out.regions)
    partial, rp = mio.revolve(
        mio.rectangle(dims=(2, 2)), angle=-90.0, num_layers=3, return_report=True
    )
    assert not rp["full"] and len(_sides(partial, "start")) == 4


def test_refusals_are_by_name():
    rect = mio.rectangle(dims=(2, 2))
    with pytest.raises(ValueError, match="parallel"):
        mio.extrude(rect, (0, 1, 0))
    with pytest.raises(ValueError, match="max_cells"):
        mio.extrude(rect, num_layers=10, max_cells=39)
    with pytest.raises(ValueError, match="last offset"):
        mio.extrude(rect, offsets=[0.5, 0.9])
    with pytest.raises(ValueError, match="direction"):
        mio.extrude(rect, (0, 1))
    with pytest.raises(NotImplementedError, match="volume"):
        mio.extrude(mio.box())
    with pytest.raises(ValueError, match="crosses the axis"):
        mio.revolve(mio.rectangle(dims=(1, 1), size=(2.0, 1.0), origin=(-1.0, 0, 0)))
    with pytest.raises(ValueError, match="both sides"):
        mio.revolve(mio.rectangle(dims=(2, 1), size=(2.0, 1.0), origin=(-1.0, 0, 0)))
    with pytest.raises(ValueError, match="3 steps"):
        mio.revolve(rect, num_layers=2)
