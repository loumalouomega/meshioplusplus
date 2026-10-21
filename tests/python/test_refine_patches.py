"""refine_patches: variational shape approximation of a segmentation.

The operation lives in the C++ core alone (no numpy twin): the argument checks
and the refusal without the core run everywhere, the rest against a core that
carries ``refine_patches``.
"""

import json
import sys

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus._cli._main import main as cli_main

from .analytic_surfaces import box_skin, filleted_plate, turned_part

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python install
    _core = None

needs_core = pytest.mark.skipif(
    not hasattr(_core, "refine_patches"), reason="core predates refine_patches"
)

REPORT_KEYS = {
    "num_patches",
    "num_iterations",
    "num_splits",
    "num_merges",
    "num_absorbed",
    "num_primitive",
    "num_freeform",
    "num_facets",
    "tolerance",
    "explained_area_fraction",
    "initial_explained_area_fraction",
    "hausdorff_max",
    "hausdorff_rms",
    "initial_energy",
    "energy",
    "stop",
}

# The arrays refinement writes, compared byte for byte between runs.
OUTPUT_ARRAYS = (
    "refine:patch_id",
    "refine:patch_proxy",
    "refine:patch_iterations",
    "refine:patch_rms",
    "refine:patch_max",
    "refine:patch_cells",
    "refine:patch_demotion",
    "refine:explained_area_fraction",
    "refine:hausdorff",
    "fit:patch_id",
    "fit:patch_type",
    "fit:patch_params",
)


def _labels(mesh):
    return np.concatenate([np.asarray(a) for a in mesh.cell_data["segment:patch"]])


@pytest.fixture(scope="module")
def segmented_box():
    # Every edge of the box is smooth at 100 degrees: one patch no primitive
    # explains, which the RANSAC split takes apart.
    return mio.segment(box_skin(3), feature_angle=100.0)


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"tolerance": 0.0}, "tolerance"),
        ({"tolerance": float("nan")}, "tolerance"),
        ({"normal_weight": -1.0}, "normal_weight"),
        ({"epsilon": -1.0}, "epsilon"),
        ({"epsilon": float("inf")}, "epsilon"),
        ({"feature_angle": 200.0}, "feature angle"),
        ({"max_iterations": -1}, "max_iterations"),
        ({"min_patch_cells": 0}, "min_patch_cells"),
        ({"ransac_seed": -1}, "ransac_seed"),
        ({"ransac_seed": 1 << 64}, "ransac_seed"),
        ({"ransac_seed": 1.5}, "ransac_seed"),
        ({"ransac_seed": True}, "ransac_seed"),
        ({"types": ["blob"]}, "unknown primitive type 'blob'"),
        ({"prefix": ""}, "prefix"),
    ],
)
def test_bad_arguments_are_refused_before_the_core(kwargs, match):
    seg = mio.segment(box_skin(2))
    with pytest.raises(ValueError, match=match) as info:
        mio.refine_patches(seg, **kwargs)
    assert str(info.value).startswith("meshio++: refine_patches: ")


def test_an_unsegmented_mesh_is_refused():
    with pytest.raises(ValueError, match="run segment first"):
        mio.refine_patches(box_skin(2))


def test_without_core_the_refusal_is_named(monkeypatch):
    """The operation has no numpy twin: with the core absent it says so."""
    seg = mio.segment(box_skin(2))
    monkeypatch.setitem(sys.modules, "meshioplusplus._core", None)
    monkeypatch.delattr(mio, "_core", raising=False)
    with pytest.raises(NotImplementedError, match="compiled C\\+\\+ core"):
        mio.refine_patches(seg)


def test_fit_cad_refine_options_are_checked_before_the_core():
    for kwargs, match in (
        ({"refine_tolerance": 0.0}, "refine_tolerance"),
        ({"normal_weight": -1.0}, "normal_weight"),
        ({"refine_epsilon": -1.0}, "refine_epsilon"),
        ({"max_refine_iterations": -1}, "max_refine_iterations"),
        ({"ransac_seed": -1}, "ransac_seed"),
    ):
        with pytest.raises(ValueError, match=match):
            mio.fit_cad(box_skin(2), refine_segmentation=True, **kwargs)


@needs_core
def test_an_unchanged_partition_comes_back_as_it_went_in():
    seg = mio.segment(box_skin(3))
    out, report = mio.refine_patches(seg, return_report=True)
    assert set(report) == REPORT_KEYS
    assert report["stop"] == "unchanged"
    assert report["num_patches"] == report["num_primitive"] == 6
    assert report["explained_area_fraction"] == pytest.approx(1.0, abs=1e-12)
    assert report["hausdorff_max"] <= report["tolerance"]
    np.testing.assert_array_equal(_labels(out), _labels(seg))
    # The geometry is untouched.
    np.testing.assert_array_equal(out.points, seg.points)
    for a, b in zip(out.cells, seg.cells):
        assert a.type == b.type
        np.testing.assert_array_equal(a.data, b.data)


@needs_core
def test_box_as_one_patch_is_split_into_its_six_sides(segmented_box):
    assert len(np.unique(_labels(segmented_box))) == 1
    out, report = mio.refine_patches(
        segmented_box, feature_angle=100.0, max_iterations=20, return_report=True
    )
    assert report["num_splits"] > 0
    assert report["initial_explained_area_fraction"] == 0.0
    assert report["num_patches"] == report["num_primitive"] == 6
    assert report["explained_area_fraction"] == pytest.approx(1.0, abs=1e-12)
    assert len(np.unique(_labels(out))) == 6
    np.testing.assert_array_equal(out.field_data["refine:patch_proxy"], [1] * 6)
    assert set(np.concatenate(out.cell_data["fit:type"]).tolist()) == {1}
    for name in OUTPUT_ARRAYS:
        assert name in out.field_data, name
    # Without splitting nothing can be found.
    _, off = mio.refine_patches(
        segmented_box,
        feature_angle=100.0,
        max_iterations=20,
        allow_split=False,
        return_report=True,
    )
    assert off["num_primitive"] == 0


@needs_core
def test_tangent_junctions_of_an_angle_segmentation_are_explained():
    seg = mio.segment(turned_part())
    _, report = mio.refine_patches(seg, return_report=True)
    assert report["explained_area_fraction"] >= 0.99
    assert report["explained_area_fraction"] > report["initial_explained_area_fraction"]
    assert report["num_splits"] > 0
    assert report["hausdorff_max"] <= report["tolerance"] * (1.0 + 1e-12)


@needs_core
def test_a_seed_gives_identical_output(segmented_box):
    kwargs = {"feature_angle": 100.0, "max_iterations": 20, "ransac_seed": 12345}
    a, ra = mio.refine_patches(segmented_box, return_report=True, **kwargs)
    b, rb = mio.refine_patches(segmented_box, return_report=True, **kwargs)
    assert ra == rb
    assert _labels(a).tobytes() == _labels(b).tobytes()
    for name in OUTPUT_ARRAYS:
        x = np.asarray(a.field_data[name])
        y = np.asarray(b.field_data[name])
        assert x.dtype == y.dtype and x.tobytes() == y.tobytes(), name


@needs_core
def test_types_and_prefix_are_honoured(segmented_box):
    out, report = mio.refine_patches(
        segmented_box,
        feature_angle=100.0,
        max_iterations=20,
        types=["plane"],
        prefix="face",
        return_report=True,
    )
    assert report["num_primitive"] == 6
    assert set(out.field_data["refine:patch_proxy"].tolist()) == {1}
    assert any(r.name.startswith("face:") for r in out.regions)


@needs_core
def test_cli_round_trip(tmp_path, capsys):
    src = tmp_path / "box.vtu"
    dst = tmp_path / "refined.vtu"
    mio.write(src, box_skin(3))
    # The input has no segment:patch, so the verb segments it first, at the
    # same feature angle.
    args = ["refine-patches", str(src), str(dst), "--feature-angle", "100"]
    assert cli_main(args + ["--max-iterations", "20", "--json"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert set(report) == REPORT_KEYS
    assert report["num_patches"] == report["num_primitive"] == 6
    out = mio.read(dst)
    assert len(np.unique(_labels(out))) == 6
    assert "refine:patch_proxy" in out.field_data
    # A hexadecimal seed, the flags and the human summary.
    flags = ["--ransac-seed", "0x2a", "--no-merge", "--no-split", "--absolute"]
    assert cli_main(args + flags + ["--tolerance", "0.01"]) == 0
    assert "refine_patches:" in capsys.readouterr().out


@needs_core
def test_to_cad_refine_flags(tmp_path, capsys):
    src = tmp_path / "plate.vtu"
    dst = tmp_path / "plate.step"
    mio.write(src, mio.segment(filleted_plate()))
    base = ["to-cad", str(src), str(dst), "--backend", "native", "--json"]
    assert cli_main(base) in (0, 2)
    plain = json.loads(capsys.readouterr().out)
    assert "refined" not in plain and "explained_area_fraction" not in plain
    flags = ["--refine", "--ransac-seed", "0x2a", "--max-refine-iterations", "5"]
    assert cli_main(base + flags) in (0, 2)
    refined = json.loads(capsys.readouterr().out)
    assert refined["refined"] is True
    assert refined["refinement_kept"] is True
    assert (
        refined["explained_area_fraction"] >= refined["initial_explained_area_fraction"]
    )
    for key in ("refine_iterations", "num_refine_splits", "num_refine_merges"):
        assert key in refined


@needs_core
def test_pipeline_step_in_both_engines(tmp_path):
    src = tmp_path / "part.vtu"
    mio.write(src, turned_part(), compression=None)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {"Op": "Segment"},
            {"Op": "RefinePatches", "RansacSeed": 42, "MaxIterations": 8},
        ],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    reports = [mio.run_pipeline(settings)]
    outputs = [mio.read(tmp_path / "py.vtu")]
    if hasattr(_core, "run_pipeline_json"):
        settings["Output"]["Path"] = str(tmp_path / "cpp.vtu")
        reports.append(_core.run_pipeline_json(json.dumps(settings)))
        outputs.append(mio.read(tmp_path / "cpp.vtu"))
    for report in reports:
        entry = report["steps"][-1]
        assert entry["op"] == "RefinePatches"
        assert entry["ExplainedAreaFraction"] >= 0.99
        assert entry["Stop"] in (0.0, 1.0, 2.0, 3.0)
        for key in (
            "NumPatches",
            "NumIterations",
            "NumSplits",
            "NumMerges",
            "NumPrimitive",
            "NumFreeform",
            "HausdorffMax",
        ):
            assert key in entry, key
    if len(reports) == 2:
        py, cpp = (r["steps"][-1] for r in reports)
        assert {k: py[k] for k in py if k[:1].isupper()} == {
            k: cpp[k] for k in cpp if k[:1].isupper()
        }
        # The two engines run the same core: the arrays agree exactly (the
        # files themselves may differ in the writer's own bytes).
        a, b = outputs
        np.testing.assert_array_equal(_labels(a), _labels(b))
        for name in OUTPUT_ARRAYS:
            np.testing.assert_array_equal(
                np.asarray(a.field_data[name]), np.asarray(b.field_data[name])
            )
    # The seed must be exact in a JSON number: above 2^53 it is refused.
    bad = dict(settings, Operations=[{"Op": "RefinePatches", "RansacSeed": 2**53 + 1}])
    with pytest.raises(Exception, match="RansacSeed"):
        mio.run_pipeline(bad)
    bad = dict(settings, Operations=[{"Op": "RefinePatches", "Types": ["blob"]}])
    with pytest.raises(Exception, match="unknown primitive type 'blob'"):
        mio.run_pipeline(bad)
