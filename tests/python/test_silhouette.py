"""silhouette (roadmap 7.6.3): the outline of a part along a direction."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "silhouette"), reason="needs the C++ core"
)


def test_a_sphere_is_one_circle_of_its_radius():
    n = np.array([1.0, 2.0, 3.0])
    out, rep = mp.silhouette(mp.icosa_sphere(2.0, 3), n, return_report=True)
    assert (rep["num_loops"], rep["num_outer"], rep["num_holes"]) == (1, 1, 0)
    n /= np.linalg.norm(n)
    h = out.points @ n
    r = np.linalg.norm(out.points - np.outer(h, n), axis=1)
    assert np.all(r <= 2.0 + 1e-12) and np.all(r >= 2.0 * 0.98)
    np.testing.assert_allclose(h, rep["min_height"])
    assert out.cell_data["section:chain"][0].dtype == np.int64
    np.testing.assert_allclose(out.field_data["silhouette:normal"], n)


def test_holes_offsets_and_the_prism():
    plate = mp.annulus(0.5, 1.0, 64, 2, cell_type="triangle")
    out, rep = mp.silhouette(plate, return_report=True)
    assert (rep["num_outer"], rep["num_holes"]) == (1, 1)
    assert sorted(set(out.cell_data["silhouette:hole"][0].tolist())) == [0, 1]
    sketch, fr = mp.fit_curves(out, tolerance=1e-3, return_report=True)
    assert fr["num_circles"] == 2
    cube = mp.box((2, 2, 2))
    _, mitre = mp.silhouette(cube, clearance=0.1, round_angle=80, return_report=True)
    assert mitre["area"] == pytest.approx(1.44, abs=1e-12)
    _, press = mp.silhouette(cube, clearance=-0.1, return_report=True)
    assert press["area"] == pytest.approx(0.64, abs=1e-12)
    prism, prep = mp.silhouette(
        mp.icosa_sphere(1.0, 2), extrude=True, return_report=True
    )
    tris = prism.cells[0].data
    edges = {}
    for t in tris:
        for k in range(3):
            e = (int(t[k]), int(t[(k + 1) % 3]))
            edges[e] = edges.get(e, 0) + 1
    assert all(n == 1 and (b, a) in edges for (a, b), n in edges.items())
    assert (prep["min_height"], prep["max_height"]) == (-1.0, 1.0)


def test_both_pipeline_engines_agree(tmp_path):
    src = tmp_path / "in.vtu"
    mp.write(src, mp.icosa_sphere(1.0, 2))
    step = {"Op": "Silhouette", "Normal": [0, 1, 1], "Clearance": 0.05}
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [step],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep = mp.run_pipeline(settings)
    assert rep["steps"][0]["NumLoops"] == 1
    settings["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    cpp = _core.run_pipeline_json(json.dumps(settings))
    assert cpp["steps"][0]["Area"] == rep["steps"][0]["Area"]
    assert mp.meshes_equal(
        mp.read(str(tmp_path / "py.vtu")), mp.read(str(tmp_path / "cpp.vtu"))
    )


def test_an_outline_crosses_to_dxf(tmp_path):
    sketch = mp.fit_curves(mp.silhouette(mp.cylinder(1.0, 2.0, 8, 2)), tolerance=1e-3)
    path = tmp_path / "outline.dxf"
    mp.write(path, sketch)
    back = mp.read(path)
    assert len(back.cells[0].data) > 0


def test_refusals():
    with pytest.raises(ValueError, match="normal"):
        mp.silhouette(mp.icosa_sphere(), (0, 0, 0))
    with pytest.raises(ValueError, match="arc_segments"):
        mp.silhouette(mp.icosa_sphere(), arc_segments=2)
    lines = mp.Mesh(np.eye(3), [("line", np.array([[0, 1], [1, 2]]))])
    with pytest.raises(ValueError, match="no surface cells"):
        mp.silhouette(lines)


def test_the_mcp_tool(tmp_path):
    from meshioplusplus.mcp import _tools

    src = tmp_path / "s.vtu"
    mp.write(src, mp.icosa_sphere(1.0, 2))
    rep = _tools.tool_silhouette(
        str(src), output_path=str(tmp_path / "o.vtu"), clearance=0.1
    )
    assert rep["num_loops"] == 1 and "output_path" in rep
