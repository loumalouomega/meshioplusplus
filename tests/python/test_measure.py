"""measure, align_datums and deviation against an element (roadmap 7.2, 7.3.4,
7.4.2)."""

import json
import math

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "measure"), reason="needs the C++ core"
)


def block(h=20.0, tilt_deg=0.0, turn=0.0, shift=(0.0, 0.0, 0.0)):
    """A block of quads, one cell region per face; the xmax face leans."""
    d = h * math.tan(math.radians(tilt_deg))
    pts = np.array(
        [
            [0, 0, 0],
            [2, 0, 0],
            [2, 3, 0],
            [0, 3, 0],
            [0, 0, h],
            [2 + d, 0, h],
            [2 + d, 3, h],
            [0, 3, h],
        ],
        dtype=float,
    )
    c, s = math.cos(turn), math.sin(turn)
    pts[:, :2] = pts[:, :2] @ np.array([[c, s], [-s, c]])
    pts += shift
    quads = np.array(
        [
            [0, 3, 2, 1],
            [4, 5, 6, 7],
            [0, 4, 7, 3],
            [1, 2, 6, 5],
            [0, 1, 5, 4],
            [3, 7, 6, 2],
        ]
    )
    mesh = mp.Mesh(pts, [("quad", quads)])
    mesh.regions = [
        mp.Region(name, "cell", np.array([k]))
        for k, name in enumerate(["bottom", "top", "xmin", "xmax", "ymin", "ymax"])
    ]
    return mesh


SPEC = {
    "elements": [
        {"name": "A", "fit": "plane", "region": "bottom"},
        {"name": "B", "fit": "plane", "region": "xmax"},
        {"name": "C", "fit": "plane", "region": "ymin"},
        {"name": "T", "fit": "plane", "region": "top"},
        {"name": "edge", "construct": "intersection", "of": ["A", "B"]},
    ],
    "dimensions": [
        {
            "name": "height",
            "type": "distance",
            "between": ["A", "T"],
            "nominal": 20,
            "tolerance": 0.01,
        }
    ],
    "tolerances": [
        {
            "name": "perp",
            "type": "perpendicularity",
            "feature": "B",
            "datum": "A",
            "limit": 0.01,
        },
        {"name": "flat", "type": "flatness", "feature": "A"},
    ],
}


def test_the_tilted_side_fails_perpendicularity():
    out, rep = mp.measure(block(tilt_deg=0.05), SPEC)
    perp = rep["tolerances"][0]
    assert perp["value"] == pytest.approx(20.0 * math.tan(math.radians(0.05)), abs=1e-9)
    assert perp["pass"] is False
    assert rep["dimensions"][0]["pass"] is True
    assert rep["num_fail"] == 1
    table = mp.measure_table(out)
    assert [e["name"] for e in table] == ["A", "B", "C", "T", "edge"]
    assert table[4]["kind"] == "line"
    assert "measure:B" in [r.name for r in out.regions]


def test_spec_as_json_text_and_file(tmp_path):
    path = tmp_path / "spec.json"
    path.write_text(json.dumps(SPEC))
    _, a = mp.measure(block(), str(path))
    _, b = mp.measure(block(), json.dumps(SPEC))
    assert a == b
    with pytest.raises(ValueError, match="no element named"):
        mp.measure(block(), {"dimensions": [{"type": "angle", "between": ["X", "Y"]}]})


def test_align_datums_puts_the_faces_on_the_coordinate_planes():
    measured, _ = mp.measure(block(4.0, 0.0, 0.7, (5.0, -3.0, 2.0)), SPEC)
    aligned, matrix = mp.align_datums(measured, "A", "B", "C")
    p = aligned.points
    np.testing.assert_allclose(p[[0, 1, 2, 3], 2], 0.0, atol=1e-12)
    np.testing.assert_allclose(p[[1, 2, 5, 6], 0], 0.0, atol=1e-12)
    np.testing.assert_allclose(p[[0, 1, 4, 5], 1], 0.0, atol=1e-12)
    assert matrix.shape == (4, 4)
    with pytest.raises(ValueError, match="parallel"):
        mp.align_datums(measured, "A", "T")


def test_deviation_against_an_element():
    measured, _ = mp.measure(block(), SPEC)
    out, rep = mp.deviation(measured, measured, element="T", return_report=True)
    d = out.point_data["deviation:distance"]
    np.testing.assert_allclose(d[[4, 5, 6, 7]], 0.0, atol=1e-12)
    assert rep["max"] == pytest.approx(0.0, abs=1e-12)
    with pytest.raises(ValueError, match="edge"):
        mp.deviation(measured, measured, element="edge")


def test_both_pipeline_engines(tmp_path):
    src = tmp_path / "block.vtu"
    mp.write(src, block(tilt_deg=0.05), compression=None)
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps(SPEC))
    ops = [
        {"Op": "Measure", "Spec": str(spec)},
        {"Op": "Deviation", "Element": "T"},
        {"Op": "AlignDatums", "Primary": "A", "Secondary": "B", "Origin": "C"},
    ]
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": ops,
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep = mp.run_pipeline(settings)
    assert rep["steps"][0]["NumFail"] == 1.0
    settings["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    cpp = _core.run_pipeline_json(json.dumps(settings))
    assert cpp["steps"][0]["NumFail"] == 1.0
    a, b = mp.read(str(tmp_path / "py.vtu")), mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(a.points, b.points)
    settings["Operations"] = [{"Op": "Measure", "Spec": str(spec), "Require": True}]
    with pytest.raises(Exception, match="Measure failed"):
        mp.run_pipeline(settings)


def test_the_mcp_tools(tmp_path):
    from meshioplusplus.mcp import _tools

    src = tmp_path / "block.vtu"
    mp.write(src, block(), compression=None)
    rep = _tools.tool_measure(str(src), SPEC, output_path=str(tmp_path / "m.vtu"))
    assert rep["num_fail"] == 0
    al = _tools.tool_align_datums(
        str(tmp_path / "m.vtu"), str(tmp_path / "a.vtu"), "A", "B", "C"
    )
    assert len(al["matrix"]) == 4
