"""fit_curves (roadmap 7.6.2): lines, arcs and circles from chained sections."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "fit_curves"), reason="needs the C++ core"
)


def _prism(outline):
    n = len(outline)
    pts = [(x, y, 0.0) for x, y in outline] + [(x, y, 1.0) for x, y in outline]
    tris = []
    for i in range(n):
        j = (i + 1) % n
        tris += [(i, j, j + n), (i, j + n, i + n)]
    return mp.Mesh(np.array(pts, float), [("triangle", np.array(tris))])


def _slot(a=2.0, r=0.75, k=24):
    out = []
    for t in np.linspace(-np.pi / 2, np.pi / 2, k + 1):
        out.append((a + r * np.cos(t), r * np.sin(t)))
    for t in np.linspace(np.pi / 2, 3 * np.pi / 2, k + 1):
        out.append((-a + r * np.cos(t), r * np.sin(t)))
    return out


def _section(mesh):
    return mp.slice(mesh, (0, 0, 0.5), (0, 0, 1), chain=True)


def test_a_slot_is_two_lines_and_two_arcs():
    out, rep = mp.fit_curves(
        _section(_prism(_slot())), tolerance=1e-6, return_report=True
    )
    assert (rep["num_lines"], rep["num_arcs"], rep["num_circles"]) == (2, 2, 0)
    arcs = [p for p in rep["pieces"] if p["type"] == "arc"]
    for p in arcs:
        assert p["radius"] == pytest.approx(0.75, abs=1e-6)
        assert abs(p["sweep"]) == pytest.approx(180.0, abs=1e-4)
    assert np.all(out.cell_data["sketch:piece"][0] >= 0)
    assert out.field_data["sketch:params"].shape == (4, 10)
    assert len(rep["frames"]) == 1
    assert json.loads(mp._fit_curves.sketch_json(rep))["pieces"][0]["type"] in (
        "line",
        "arc",
    )


def test_refusals():
    with pytest.raises(ValueError, match="line mesh"):
        mp.fit_curves(mp.icosa_sphere(1.0, 1))
    with pytest.raises(ValueError, match="section:chain"):
        mp.fit_curves(mp.Mesh(np.eye(3), [("line", np.array([[0, 1], [1, 2]]))]))


def test_cli_pipeline_and_mcp(tmp_path):
    from meshioplusplus.mcp import _tools

    sec = str(tmp_path / "sec.vtu")
    mp.write(sec, _section(_prism(_slot())), compression=None)
    assert mp._cli.main(["fit-curves", sec, "--tolerance", "1e-6"]) in (0, None)
    res = _tools.tool_fit_curves(sec, tolerance=1e-6)
    assert res["num_arcs"] == 2
    src = str(tmp_path / "slot.vtu")
    mp.write(src, _prism(_slot()), compression=None)
    doc = {
        "Version": 1,
        "Input": {"Path": src},
        "Operations": [
            {"Op": "Slice", "Point": [0, 0, 0.5], "Normal": [0, 0, 1], "Chain": True},
            {"Op": "FitCurves", "Tolerance": 1e-6},
        ],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep_py = mp.run_pipeline(doc)
    assert rep_py["steps"][1]["NumArcs"] == 2.0
    if getattr(_core, "__has_json__", False):
        doc["Output"]["Path"] = str(tmp_path / "c.vtu")
        assert _core.run_pipeline_json(json.dumps(doc))["steps"] == rep_py["steps"]
