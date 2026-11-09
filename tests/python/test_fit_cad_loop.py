"""face_support, fit_cad_iterate and constrain_primitives (roadmap 7.7.1-7.7.3)."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "fit_cad_iterate"), reason="needs the C++ core"
)


def open_box():
    box = mp.box_surface((4, 4, 4))
    quads = box.cells[0].data
    keep = ~np.all(box.points[quads][:, :, 2] == 0.0, axis=1)
    holed = mp.Mesh(box.points, [(box.cells[0].type, quads[keep])])
    return mp.repair(holed, max_hole_edges=0, record_provenance=True)


def bent_box():
    box = mp.box_surface((8, 8, 8))
    p = box.points.copy()
    p[:, 2] = p[:, 2] * (1 + 0.06 * np.maximum(0.0, p[:, 0] - 0.5))
    return mp.Mesh(p, box.cells)


def fitted(mesh):
    return mp.fit_primitives(mp.segment(mesh))


def test_the_bottom_a_scan_missed_is_inferred():
    scan = open_box()
    solid = mp.fit_cad(fitted(scan), backend="native")
    marked, rep = mp.face_support(solid, scan, return_report=True)
    assert (rep["num_measured"], rep["num_partial"], rep["num_inferred"]) == (5, 0, 1)
    assert rep["share_within"] == pytest.approx(1.0)
    assert "fit_cad:face_support" in marked.field_data
    assert "fit_cad:support" in marked.cell_data
    _, rep2 = mp.face_support(solid, scan, hole_data="none", return_report=True)
    assert rep2["num_measured"] == 6
    with pytest.raises(ValueError, match="fit_cad:face"):
        mp.face_support(scan, scan)


def test_iterate_splits_a_merged_patch():
    out, rep = mp.fit_cad_iterate(
        bent_box(),
        segment_tolerance=0.02,
        fit_tolerance=1e-3,
        target_explained=0.999,
        refine_factor=0.1,
        backend="native",
        return_report=True,
    )
    rounds = rep["rounds"]
    assert len(rounds) == 2
    assert rounds[0]["explained"] < 0.9 < rounds[1]["explained"]
    assert rounds[0]["resegmented"]
    assert rep["stop_reason"] == "target" and rep["target_met"]
    assert rep["best_round"] == 1
    _, again = mp.fit_cad_iterate(
        bent_box(),
        segment_tolerance=0.02,
        fit_tolerance=1e-3,
        target_explained=0.999,
        refine_factor=0.1,
        backend="native",
        return_report=True,
    )
    assert [r["explained"] for r in again["rounds"]] == [r["explained"] for r in rounds]
    with pytest.raises(ValueError, match="part_class"):
        mp.fit_cad_iterate(bent_box(), part_class="plastic")


def test_constrain_primitives_forbid_and_report():
    rng = np.random.default_rng(3)
    box = mp.box_surface((4, 4, 4))
    jittered = mp.Mesh(
        box.points + rng.uniform(-2e-5, 2e-5, box.points.shape), box.cells
    )
    skin = fitted(jittered)
    out, rep = mp.constrain_primitives(skin, return_report=True)
    assert rep["num_constrained"] > 0
    assert rep["relations"] and not rep["unmet"]
    target = rep["relations"][-1]["patch_b"]
    ids = [int(i) for i in skin.field_data["fit:patch_id"]]
    forbid = [
        (i, target, kind)
        for i in ids
        if i != target
        for kind in ("parallel", "perpendicular")
    ]
    out2, rep2 = mp.constrain_primitives(skin, forbid=forbid, return_report=True)
    row = ids.index(target)
    np.testing.assert_array_equal(
        out2.field_data["fit:patch_params"][row],
        skin.field_data["fit:patch_params"][row],
    )
    assert all(r["patch_b"] != target for r in rep2["relations"])
    _, rep3 = mp.constrain_primitives(
        skin, force=[(ids[0], ids[1], "coaxial")], return_report=True
    )
    assert rep3["unmet"] == [(ids[0], ids[1], "coaxial")]
    with pytest.raises(ValueError, match="relation"):
        mp.constrain_primitives(skin, force=[(ids[0], ids[1], "tangent")])


def test_both_pipeline_engines(tmp_path):
    src = tmp_path / "bent.vtu"
    mp.write(src, bent_box(), compression=None)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {
                "Op": "FitCadIterate",
                "SegmentTolerance": 0.02,
                "FitTolerance": 0.001,
                "TargetExplained": 0.999,
                "RefineFactor": 0.1,
                "Backend": "native",
            },
            {"Op": "FaceSupport", "Scan": str(src)},
        ],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep = mp.run_pipeline(settings)
    assert rep["steps"][0]["NumRounds"] == 2.0
    assert rep["steps"][1]["NumMeasured"] == 7.0
    settings["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    cpp = _core.run_pipeline_json(json.dumps(settings))
    assert cpp["steps"][0]["NumRounds"] == 2.0
    a, b = mp.read(str(tmp_path / "py.vtu")), mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(a.points, b.points)


def test_the_mcp_tools(tmp_path):
    from meshioplusplus.mcp import _tools

    src = tmp_path / "bent.vtu"
    mp.write(src, bent_box(), compression=None)
    rep = _tools.tool_fit_cad_iterate(
        str(src), str(tmp_path / "solid.vtu"), backend="native"
    )
    assert rep["stop_reason"] == "target"
    rep = _tools.tool_face_support(str(tmp_path / "solid.vtu"), str(src))
    assert rep["num_measured"] >= 6
    mp.write(tmp_path / "skin.vtu", fitted(bent_box()))
    rep = _tools.tool_constrain_primitives(
        str(tmp_path / "skin.vtu"), str(tmp_path / "c.vtu")
    )
    assert "relations" in rep
