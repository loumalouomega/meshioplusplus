"""Tests for ``segment`` (roadmap §6.3.1): the Python face, its numpy twin
pinned against the compiled core, the CLI verb and the pipeline step. The C++
behaviour is pinned in ``tests/cpp/test_segment.cpp``."""

import json

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _segment as sg

from .analytic_surfaces import box_skin, cylinder, sphere_band, torus, turned_part

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")


def _opts(**kw):
    opts = {
        "feature_angle": 30.0,
        "normal_tolerance": 0.0,
        "curvature_tolerance": 0.0,
        "min_patch_cells": 1,
        "region": "",
        "prefix": "patch",
    }
    opts.update(kw)
    return opts


def capped_cylinder():
    side = cylinder(m=4)
    n = 32
    pts = np.vstack([side.points, [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]]])
    bottom, top, last = len(side.points), len(side.points) + 1, 4 * n
    tris = []
    for k in range(n):
        k1 = (k + 1) % n
        tris.append([bottom, k1, k])
        tris.append([top, last + k, last + k1])
    return mio.Mesh(pts, [("quad", side.cells[0].data), ("triangle", np.array(tris))])


def test_box_gives_six_faces_and_twelve_adjacencies():
    out, report = mio.segment(box_skin(3), return_report=True)
    assert report == {
        "num_patches": 6,
        "num_merged": 0,
        "num_feature_edges": 36,
        "num_boundary_edges": 0,
        "num_non_manifold_edges": 0,
    }
    adj = out.field_data["segment:adjacency"]
    assert adj.shape == (12, 3)
    assert np.all(adj[:, 0] < adj[:, 1]) and np.all(adj[:, 2] == 3)
    assert np.bincount(out.cell_data["segment:patch"][0]).tolist() == [9] * 6
    names = sorted(r.name for r in out.regions)
    assert names == [f"patch:{k}" for k in range(6)]


def test_capped_cylinder_gives_three_patches_and_points_are_untouched():
    mesh = capped_cylinder()
    out = mio.segment(mesh)
    side, caps = out.cell_data["segment:patch"]
    assert len(set(side.tolist())) == 1 and len(set(caps.tolist())) == 2
    np.testing.assert_array_equal(out.points, mesh.points)
    for a, b in zip(out.cells, mesh.cells):
        np.testing.assert_array_equal(a.data, b.data)


def test_tolerances_and_merging():
    assert mio.segment(cylinder(), return_report=True)[1]["num_patches"] == 1
    _, rep = mio.segment(cylinder(), normal_tolerance=40.0, return_report=True)
    assert rep["num_patches"] > 1
    _, before = mio.segment(torus(), normal_tolerance=40.0, return_report=True)
    _, after = mio.segment(
        torus(), normal_tolerance=40.0, min_patch_cells=40, return_report=True
    )
    assert after["num_merged"] > 0
    assert after["num_patches"] + after["num_merged"] == before["num_patches"]


def test_bad_options_and_volumes_are_refused():
    with pytest.raises(ValueError, match="feature angle"):
        mio.segment(box_skin(1), feature_angle=200.0)
    with pytest.raises(ValueError, match="min_patch_cells"):
        mio.segment(box_skin(1), min_patch_cells=0)
    with pytest.raises(ValueError, match="no cell region"):
        mio.segment(box_skin(1), region="nope")
    tet = mio.Mesh(np.eye(4, 3), [("tetra", np.array([[0, 1, 2, 3]]))])
    with pytest.raises(ValueError, match="volume block"):
        mio.segment(tet)
    with pytest.raises(ValueError, match="volume block"):
        sg._segment_py(tet, _opts())


@needs_core
@pytest.mark.parametrize(
    "mesh",
    [box_skin(4), capped_cylinder(), torus(), sphere_band(), turned_part(noise=1e-3)],
    ids=["box", "capped_cylinder", "torus", "sphere", "noisy_part"],
)
@pytest.mark.parametrize(
    "kw",
    [
        {},
        {"curvature_tolerance": 0.3},
        {"normal_tolerance": 35.0, "min_patch_cells": 6},
        {"feature_angle": 10.0},
    ],
    ids=["default", "curvature", "normal+merge", "tight"],
)
def test_twin_matches_the_core(mesh, kw):
    core, core_report = mio.segment(mesh, return_report=True, **kw)
    twin, twin_report = sg._segment_py(mesh, _opts(**kw))
    assert core_report == twin_report
    for a, b in zip(core.cell_data["segment:patch"], twin.cell_data["segment:patch"]):
        np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(
        core.field_data["segment:adjacency"], twin.field_data["segment:adjacency"]
    )
    ca = {r.name: np.asarray(r.entries) for r in core.regions}
    tb = {r.name: np.asarray(r.entries) for r in twin.regions}
    assert ca.keys() == tb.keys()
    for name in ca:
        np.testing.assert_array_equal(ca[name], tb[name])


def test_rerun_replaces_the_patch_regions():
    many = mio.segment(cylinder(), normal_tolerance=40.0)
    one = mio.segment(many)
    assert [r.name for r in one.regions if r.name.startswith("patch:")] == ["patch:0"]


def test_cli_verb(tmp_path, capsys):
    src = str(tmp_path / "box.vtu")
    mio.write(src, box_skin(2))
    out = str(tmp_path / "seg.vtu")
    assert mio._cli.main(["segment", src, out, "--json"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["num_patches"] == 6
    back = mio.read(out)
    assert sorted(r.name for r in back.regions if r.name.startswith("patch:")) == [
        f"patch:{k}" for k in range(6)
    ]


def test_pipeline_steps_in_both_engines(tmp_path):
    # Uncompressed input: the C++ leg reads it through _core, and the Windows CI
    # build of _core has no zlib to inflate the default compressed VTU. These
    # tests are about the steps, not the codec.
    src = tmp_path / "box.vtu"
    mio.write(src, box_skin(2), compression=None)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {"Op": "Segment", "FeatureAngle": 30.0},
            {"Op": "FitPrimitives", "Types": ["plane", "cylinder"]},
        ],
        "Output": {"Path": str(tmp_path / "out.vtu")},
    }
    reports = [mio.run_pipeline(settings)]
    if _core is not None and hasattr(_core, "run_pipeline_json"):
        settings["Output"]["Path"] = str(tmp_path / "out_cpp.vtu")
        reports.append(_core.run_pipeline_json(json.dumps(settings)))
        assert mio.meshes_equal(
            mio.read(str(tmp_path / "out.vtu")), mio.read(settings["Output"]["Path"])
        )
    for report in reports:
        seg, fit = report["steps"]
        assert seg["NumPatches"] == 6
        assert fit["NumPlane"] == 6 and fit["NumUnfitted"] == 0
    bad = dict(settings, Operations=[{"Op": "FitPrimitives", "Types": ["blob"]}])
    with pytest.raises(Exception, match="unknown primitive type 'blob'"):
        mio.run_pipeline(bad)
    if _core is not None and hasattr(_core, "run_pipeline_json"):
        with pytest.raises(Exception, match="unknown primitive type 'blob'"):
            _core.run_pipeline_json(json.dumps(bad))


def test_a_noisy_turned_part_splits_at_creases_and_at_curvature():
    out, report = mio.segment(
        turned_part(noise=1e-3),
        curvature_tolerance=1.0,
        min_patch_cells=20,
        return_report=True,
    )
    assert report["num_patches"] == 6
    sizes = np.bincount(np.concatenate(out.cell_data["segment:patch"]))
    # base 256, dome 384, wall 512, fillet 512, shoulder 192, cone 384 cells
    assert sorted(sizes.tolist()) == [192, 256, 384, 384, 512, 512]
