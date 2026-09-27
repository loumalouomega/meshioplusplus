"""Tests for ``check_quality`` -- the quality gate (v16.24.0). The C++ side is
pinned in ``tests/cpp/test_quality_gate.cpp``; here the Python face, its numpy
twin against the core, the pipeline step and the specification grammar."""

import math

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _quality_gate as qg

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")


def two_triangles():
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [3, 0, 0], [5, 0, 0], [4, 0.1, 0]], float
    )
    return mio.Mesh(pts, [("triangle", np.array([[0, 1, 2], [3, 4, 5]]))])


def test_spec_grammar():
    t = qg.parse_quality_thresholds(
        "scaled_jacobian >= 0.2; quality:aspect_ratio<=5 @ 1%\n# comment\nskewness <= 0.5 @ 0.25"
    )
    assert [x["metric"] for x in t] == [
        "quality:scaled_jacobian",
        "quality:aspect_ratio",
        "quality:skewness",
    ]
    assert t[0]["min"] == 0.2 and t[0]["max"] is None
    assert t[1]["max_fraction"] == pytest.approx(0.01)
    for bad in ("bogus >= 1", "min_angle > 1", "min_angle >= x"):
        with pytest.raises(ValueError):
            qg.parse_quality_thresholds(bad)


@needs_core
@pytest.mark.parametrize(
    "spec",
    [
        "",
        "min_angle >= 30",
        "a;min_angle >= 30 @ 50%",
        "aspect_ratio <= 2; max_angle <= 100",
    ],
)
def test_parser_matches_the_core(spec):
    if spec.startswith("a;"):
        spec = spec[2:]
    py = qg.parse_quality_thresholds(spec)
    core = _core.parse_quality_thresholds(spec)
    assert py == [dict(d) for d in core]


def test_gate_counts_violations():
    r = mio.check_quality(two_triangles(), "min_angle >= 30")
    c = r["checks"][0]
    assert (c["evaluated"], c["violations"], c["worst_cell"]) == (2, 1, 1)
    assert not c["passed"] and not r["passed"]
    assert mio.check_quality(two_triangles(), "min_angle >= 30 @ 50%")["checks"][0][
        "passed"
    ]
    assert [c["metric"] for c in r["checks"]][1:] == ["inverted", "degenerate"]
    assert (
        len(
            mio.check_quality(two_triangles(), max_inverted=-1, max_degenerate=-1)[
                "checks"
            ]
        )
        == 0
    )


@needs_core
@pytest.mark.parametrize(
    "spec", ["min_angle >= 30", "aspect_ratio <= 3 @ 10%", "max_angle <= 170"]
)
def test_twin_matches_the_core(spec):
    mesh = two_triangles()
    core = mio.check_quality(mesh, spec)
    twin = qg._check_quality_py(mesh, qg.parse_quality_thresholds(spec), 0, 0)
    assert core["passed"] == twin["passed"]
    for a, b in zip(core["checks"], twin["checks"]):
        for k in ("name", "metric", "evaluated", "violations", "worst_cell", "passed"):
            assert a[k] == b[k], k
        for k in ("worst", "fraction", "min", "max"):
            assert (math.isnan(a[k]) and math.isnan(b[k])) or a[k] == pytest.approx(
                b[k]
            ), k


def test_pipeline_quality_gate(tmp_path):
    src = tmp_path / "in.vtu"
    mio.write(str(src), two_triangles(), compression=None)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "QualityGate", "Require": ["min_angle >= 1"]}],
        "Output": {"Path": str(tmp_path / "out.vtu")},
    }
    report = mio.run_pipeline(settings)
    assert report["steps"][0]["NumFailed"] == 0
    settings["Operations"][0]["Require"] = ["min_angle >= 30"]
    with pytest.raises(RuntimeError, match="QualityGate failed"):
        mio.run_pipeline(settings)
