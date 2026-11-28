"""Tests for Roadmap Tier 9 Wave 3 operations (units, CRS, booleans, mixed_green)."""

from __future__ import annotations

import numpy as np
import pytest

import meshioplusplus as mio


def _quad_grid_3x3():
    pts = np.array(
        [[float(i), float(j), 0.0] for j in range(4) for i in range(4)],
        dtype=np.float64,
    )
    quads = np.array(
        [
            [j * 4 + i, j * 4 + i + 1, (j + 1) * 4 + i + 1, (j + 1) * 4 + i]
            for j in range(3)
            for i in range(3)
        ],
        dtype=np.int64,
    )
    return mio.Mesh(pts, [("quad", quads)])


def _closed_tet_surface():
    pts = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    tris = np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]], dtype=np.int64)
    return mio.Mesh(pts, [("triangle", tris)])


def test_convert_units_roundtrip_and_strict():
    pts = np.array(
        [[0.0, 0.0, 0.0], [1000.0, 0.0, 0.0], [0.0, 1000.0, 0.0], [0.0, 0.0, 1000.0]],
        dtype=np.float64,
    )
    m = mio.Mesh(
        pts,
        [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
        point_data={"stress": np.full(4, 250.0, dtype=np.float64)},
    )
    mio._units.set_length_unit(m, 1e-3)
    mio.set_field_unit(m, "stress", "MPa")
    assert mio.field_unit(m, "stress")[0] == pytest.approx(1e6)

    si, report = mio.convert_units(m, "SI", strict=True, return_report=True)
    assert report["num_arrays_converted"] == 1
    assert si.points[1, 0] == pytest.approx(1.0)
    assert si.point_data["stress"][0] == pytest.approx(250e6)

    back = mio.convert_units(si, "mm-t-s", strict=True)
    assert back.points[1, 0] == pytest.approx(1000.0)
    assert back.point_data["stress"][0] == pytest.approx(250.0)

    m.point_data["undeclared"] = np.ones(4, dtype=np.float64)
    with pytest.raises(ValueError):
        mio.convert_units(m, "SI", strict=True)
    relaxed, rep2 = mio.convert_units(m, "SI", strict=False, return_report=True)
    assert rep2["warnings"]
    assert relaxed.point_data["undeclared"][0] == pytest.approx(1.0)


def test_mesh_crs_and_convert_coordinates_crs():
    m = _closed_tet_surface()
    assert mio.mesh_crs(m) is None
    mio.set_mesh_crs(m, "EPSG:4326")
    assert mio.mesh_crs(m) == "EPSG:4326"

    same = mio.convert_coordinates(m, "crs", target_crs="EPSG:4326")
    assert mio.mesh_crs(same) == "EPSG:4326"
    np.testing.assert_allclose(same.points, m.points)


def test_mesh_boolean_without_provider():
    a = _closed_tet_surface()
    b = _closed_tet_surface()
    if not mio.has_boolean_provider():
        with pytest.raises((RuntimeError, NotImplementedError)):
            mio.mesh_boolean(a, b, op="union")


def test_refine_mixed_green_quad():
    m = _quad_grid_3x3()
    out = mio.refine(m, cells=[4], closure="mixed_green")
    by_type = {cb.type: len(cb.data) for cb in out.cells}
    assert by_type["quad"] == 8
    assert by_type["triangle"] == 20
    assert len(out.points) == 25
