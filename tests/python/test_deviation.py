"""Tests for ``deviation`` (roadmap 7.4.1, 7.4.2): the Python face over the
compiled core, the reference forms, the CLI verb. The C++ behaviour is pinned
in ``tests/cpp/test_deviation.cpp``."""

import json

import numpy as np
import pytest

import meshioplusplus as mio

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

pytestmark = pytest.mark.skipif(
    _core is None or not hasattr(_core, "deviation"), reason="needs the compiled core"
)


def test_an_offset_sphere_reads_its_offset():
    ref = mio.icosa_sphere(1.0, 3)
    scan = mio.icosa_sphere(1.01, 3)
    out, rep = mio.deviation(scan, ref, tolerance=0.02, return_report=True)
    d = out.point_data["deviation:distance"]
    np.testing.assert_allclose(d, 0.01, atol=1e-9)
    assert rep["mean"] == pytest.approx(0.01, abs=1e-9)
    assert rep["share_within"] == 1.0
    assert rep["num_measured"] == len(scan.points)
    assert set(rep["percentiles"]) == {"p1", "p5", "p50", "p95", "p99"}
    assert len(rep["histogram_edges"]) == 21 and len(rep["histogram_counts"]) == 20
    assert (
        sum(rep["histogram_counts"]) + rep["num_outside_range"] == rep["num_measured"]
    )
    assert rep["suggested_range"] == pytest.approx(0.01)
    # the scan comes back otherwise unchanged
    np.testing.assert_array_equal(out.points, scan.points)


def test_a_primitive_reference_in_every_form():
    scan = mio.icosa_sphere(1.01, 3)
    params = [0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0]
    for ref in (("sphere", params), {"type": 2, "params": params}):
        _, rep = mio.deviation(scan, ref, return_report=True)
        assert rep["mean"] == pytest.approx(0.01, abs=1e-12)
    fitted = mio.fit_primitives(mio.icosa_sphere(1.0, 4), patch_data="")
    _, rep = mio.deviation(scan, fitted, patch=0, return_report=True)
    assert rep["mean"] == pytest.approx(0.01, abs=1e-4)
    with pytest.raises(ValueError, match="no fitted patch"):
        mio.deviation(scan, fitted, patch=7)
    with pytest.raises(ValueError, match="unknown primitive"):
        mio.deviation(scan, ("blob", params))


def test_the_histogram_range_and_bad_options():
    ref = mio.icosa_sphere(1.0, 3)
    scan = mio.icosa_sphere(1.01, 3)
    _, rep = mio.deviation(scan, ref, bins=4, range=(0.0, 0.005), return_report=True)
    assert rep["num_outside_range"] == rep["num_measured"]
    assert rep["histogram_edges"] == pytest.approx(
        [0.0, 0.00125, 0.0025, 0.00375, 0.005]
    )
    with pytest.raises(ValueError, match="range"):
        mio.deviation(scan, ref, range=(1.0, 0.0))
    with pytest.raises(ValueError, match="bins"):
        mio.deviation(scan, ref, bins=0)
    with pytest.raises(ValueError, match="max_normal_angle"):
        mio.deviation(scan, ref, max_normal_angle=270.0)


def test_cli_verb(tmp_path, capsys):
    from meshioplusplus._cli import main

    ref, scan, out = tmp_path / "ref.vtu", tmp_path / "scan.vtu", tmp_path / "out.vtu"
    mio.write(ref, mio.icosa_sphere(1.0, 2))
    mio.write(scan, mio.icosa_sphere(1.02, 2))
    assert main([str_ for str_ in ("deviation", str(scan), str(ref), str(out))]) == 0
    assert "deviation:distance" in mio.read(out).point_data
    capsys.readouterr()
    code = main(
        [
            "deviation",
            str(scan),
            str(ref),
            "--tolerance",
            "0.05",
            "--min-share",
            "1",
            "--json",
        ]
    )
    report = json.loads(capsys.readouterr().out)
    assert code == 0 and report["passed"] is True
    code = main(
        ["deviation", str(scan), str(ref), "--tolerance", "1e-4", "--min-share", "1"]
    )
    assert code == 1


def test_deviation_regions_find_a_bump_and_a_dent():
    pts = np.array([[i, j, 0.0] for j in range(11) for i in range(11)], dtype=float)
    quads = np.array(
        [
            [j * 11 + i, j * 11 + i + 1, (j + 1) * 11 + i + 1, (j + 1) * 11 + i]
            for j in range(10)
            for i in range(10)
        ]
    )
    d = np.zeros(len(pts))
    i, j = pts[:, 0], pts[:, 1]
    d[(i >= 1) & (i <= 4) & (j >= 1) & (j <= 4)] = 2.0
    d[(i >= 6) & (i <= 9) & (j >= 6) & (j <= 9)] = -1.0
    plate = mio.Mesh(pts, [("quad", quads)], point_data={"deviation:distance": d})
    out, rep = mio.deviation_regions(plate, 0.5, return_report=True)
    assert [g["name"] for g in rep["regions"]] == [
        "deviation:excess_0",
        "deviation:missing_0",
    ]
    assert rep["regions"][0]["num_cells"] == 21
    assert rep["regions"][0]["volume"] == pytest.approx(30.0)
    assert rep["regions"][1]["volume"] == pytest.approx(-9.0)
    assert rep["regions"][0]["centroid"] == pytest.approx([2.5, 2.5, 0.0])
    names = {r.name for r in out.regions}
    assert {"deviation:excess_0", "deviation:missing_0"} <= names
    assert int(np.max(out.cell_data["deviation:region"][0])) == 1
    with pytest.raises(ValueError, match="tolerance"):
        mio.deviation_regions(plate, 0.0)
    with pytest.raises(ValueError, match="run deviation first"):
        mio.deviation_regions(mio.icosa_sphere(1.0, 1), 0.1)


def test_exaggerate_grows_an_offset_sphere():
    ref = mio.icosa_sphere(1.0, 3)
    measured = mio.deviation(mio.icosa_sphere(1.01, 3), ref)
    out, rep = mio.exaggerate(measured, scale=20.0, return_report=True)
    r = np.linalg.norm(np.asarray(out.points), axis=1)
    np.testing.assert_allclose(r, 1.01 + 20 * 0.01, atol=2e-3)
    assert rep["clamp"] == pytest.approx(0.01)
    assert rep["num_moved"] == len(out.points)
    np.testing.assert_allclose(
        np.linalg.norm(measured.points, axis=1), 1.01, atol=1e-12
    )
    with pytest.raises(ValueError, match="run deviation first"):
        mio.exaggerate(ref)


def test_colorize_bakes_banded_colours(tmp_path):
    measured = mio.deviation(mio.icosa_sphere(1.01, 3), mio.icosa_sphere(1.0, 3))
    measured.point_data["z"] = np.asarray(measured.points)[:, 2].copy()
    out, rep = mio.colorize(measured, "z", cmap="coolwarm", steps=5, return_report=True)
    rgb = np.stack([out.point_data[k] for k in ("red", "green", "blue")], axis=1)
    assert rgb.dtype == np.uint8
    assert len({tuple(c) for c in rgb}) == 5
    assert rep["location"] == "point"
    assert rep["vmin"] == pytest.approx(float(out.points[:, 2].min()))
    path = tmp_path / "c.ply"
    mio.write(path, out)
    back = mio.read(path)
    np.testing.assert_array_equal(back.point_data["red"], out.point_data["red"])
    with pytest.raises(ValueError, match="vmin"):
        mio.colorize(measured, "z", vmin=1.0)
    with pytest.raises(ValueError, match="nope"):
        mio.colorize(measured, "nope")
