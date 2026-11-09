"""profile_solid through a CAD provider (roadmap 7.7.4)."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "profile_solid"), reason="needs the C++ core"
)
needs_provider = pytest.mark.skipif(
    "occt" not in mp.profile_solid_providers(), reason="needs OpenCASCADE"
)


def sketch(corners):
    corners = np.asarray(corners, dtype=float)
    pts = []
    for i in range(len(corners)):
        for s in range(4):
            pts.append(
                corners[i] + (corners[(i + 1) % len(corners)] - corners[i]) * s / 4
            )
    pts = np.array(pts)
    n = len(pts)
    m = mp.Mesh(pts, [("line", np.array([[i, (i + 1) % n] for i in range(n)]))])
    m.cell_data["section:chain"] = [np.zeros(n, dtype=np.int64)]
    return mp.fit_curves(m, tolerance=1e-6)


def rectangle():
    return sketch([[0, 0, 0], [2, 0, 0], [2, 1, 0], [0, 1, 0]])


def test_the_text_form():
    f = mp.profile_features_from_text("extrude:2:symmetric;revolve:0,0,0:0,0,1:90")
    assert f[0] == {"kind": "extrude", "distance": 2.0, "symmetric": True}
    assert f[1]["angle"] == 90.0 and f[1]["axis_direction"] == [0.0, 0.0, 1.0]
    with pytest.raises(ValueError, match="feature"):
        mp.profile_features_from_text("loft:1")


def test_refusals():
    with pytest.raises(ValueError, match="sketch"):
        mp.profile_solid(mp.icosa_sphere(1.0, 1), "extrude:1")
    with pytest.raises(ValueError, match="first feature"):
        mp.profile_solid(rectangle(), [{"kind": "fillet", "radius": 0.1}])
    with pytest.raises(NotImplementedError, match="provider"):
        mp.profile_solid(rectangle(), "extrude:1", provider="nowhere")


@needs_provider
def test_extrude_fillet_and_shell(tmp_path):
    box, rep = mp.profile_solid(
        rectangle(), "extrude:3", tmp_path / "box.step", return_report=True
    )
    assert rep["valid"] and rep["num_faces"] == 6
    assert rep["volume"] == pytest.approx(6.0)
    assert (tmp_path / "box.step").read_bytes().startswith(b"ISO-10303-21")
    _, rep = mp.profile_solid(rectangle(), "extrude:3;fillet:0.1", return_report=True)
    assert rep["valid"] and rep["volume"] < 6.0
    _, rep = mp.profile_solid(rectangle(), "extrude:3;shell:0.1", return_report=True)
    assert rep["volume"] == pytest.approx(6.0 - 1.8 * 0.8 * 2.9, rel=1e-6)


@needs_provider
def test_the_pipeline_and_mcp(tmp_path):
    from meshioplusplus.mcp import _tools

    src = tmp_path / "sketch.vtu"
    mp.write(src, rectangle())
    rep = _tools.tool_profile_solid(str(src), str(tmp_path / "solid.vtu"), "extrude:2")
    assert rep["volume"] == pytest.approx(4.0)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "ProfileSolid", "Features": "extrude:2"}],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    assert mp.run_pipeline(settings)["steps"][0]["Volume"] == pytest.approx(4.0)
