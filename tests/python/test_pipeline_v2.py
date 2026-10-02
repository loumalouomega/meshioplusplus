"""Spatial multi-input / fan-out parity, separate from transient sequences."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core
from meshioplusplus._pipeline import _V2_OP_TABLE

from . import helpers


@pytest.fixture(params=["python", "native"])
def runner(request):
    if request.param == "python":
        return mp.run_pipeline
    if not hasattr(_core, "pipeline_v2_op_table") or not getattr(
        _core, "__has_json__", False
    ):
        pytest.skip("requires rebuilt JSON-enabled pipeline v2")
    return lambda settings: _core.run_pipeline_json(json.dumps(settings))


def test_v2_vocabulary():
    if not hasattr(_core, "pipeline_v2_op_table"):
        pytest.skip("requires rebuilt pipeline v2")
    assert _core.pipeline_v2_op_table() == {k: list(v) for k, v in _V2_OP_TABLE.items()}


def _write(path, mesh):
    mp.write(path, mesh, compression=None)
    return str(path)


def test_merge_uses_current_mesh_then_ordered_inputs(tmp_path, runner):
    a = _write(tmp_path / "a.vtu", helpers.tet_mesh)
    second = mp.transform(helpers.tet_mesh, translate=[3, 0, 0])
    b = _write(tmp_path / "b.vtu", second)
    out = tmp_path / "merged.vtu"
    report = runner({
        "Version": 2, "Input": {"Path": a}, "Output": {"Path": str(out), "Codec": "none"},
        "Operations": [{"Op": "Merge", "Inputs": [b]}, {"Op": "Quality"}],
    })
    mesh = mp.read(out)
    np.testing.assert_array_equal(mesh.points, np.vstack([helpers.tet_mesh.points, second.points]))
    count = sum(len(c) for c in helpers.tet_mesh.cells)
    np.testing.assert_array_equal(np.concatenate(mesh.cell_data["source_mesh_id"]), np.repeat([0, 1], count))
    assert report["steps"][0] == {"op": "Merge", "NumInputs": 2.0}


def test_interpolate_source_input_to_current_target(tmp_path, runner):
    target = _write(tmp_path / "target.vtu", helpers.tet_mesh)
    source_mesh = helpers.tet_mesh.copy()
    source_mesh.point_data["value"] = np.arange(len(source_mesh.points), dtype=float)
    source = _write(tmp_path / "source.vtu", source_mesh)
    out = tmp_path / "sampled.vtu"
    report = runner({
        "Version": 2, "Input": {"Path": target}, "Output": {"Path": str(out), "Codec": "none"},
        "Operations": [{"Op": "Interpolate", "Inputs": [source], "Arrays": ["value"]}],
    })
    np.testing.assert_array_equal(mp.read(out).point_data["value"].ravel(), source_mesh.point_data["value"])
    assert report["steps"] == [{"op": "Interpolate"}]


def test_undo_green_reads_coarse_input(tmp_path, runner):
    coarse = helpers.tet_mesh.copy()
    fine = mp.refine(coarse, record_levels=True, record_hierarchy=True)
    a = _write(tmp_path / "fine.vtu", fine)
    b = _write(tmp_path / "coarse.vtu", coarse)
    out = tmp_path / "undone.vtu"
    report = runner({
        "Version": 2, "Input": {"Path": a}, "Output": {"Path": str(out), "Codec": "none"},
        "Operations": [{"Op": "UndoGreen", "Inputs": [b]}],
    })
    expected, counts = mp.undo_green(coarse, fine, return_report=True)
    assert sum(len(c) for c in mp.read(out).cells) == sum(len(c) for c in expected.cells)
    assert report["steps"] == [{"op": "UndoGreen", "NumGroupsUndone": float(counts["num_groups_undone"]), "NumCellsRemoved": float(counts["num_cells_removed"])}]


@pytest.mark.parametrize("op,token", [("Split", "key"), ("Partition", "part")])
def test_patterned_fanout(tmp_path, runner, op, token):
    source = _write(tmp_path / "input.vtu", helpers.tri_quad_mesh)
    step = {"Op": op}
    if op == "Partition":
        step.update(Nparts=2, Method="sfc", RecordIds=True)
    report = runner({
        "Version": 2, "Input": {"Path": source},
        "Output": {"Pattern": str(tmp_path / ("piece_{" + token + "}.vtu")), "Codec": "none"},
        "Operations": [step],
    })
    outputs = sorted(tmp_path.glob("piece_*.vtu"))
    assert len(outputs) == 2
    assert report["steps"] == [{"op": op, "NumPieces": 2.0}]
    meshes = [mp.read(p) for p in outputs]
    assert sum(sum(len(c) for c in m.cells) for m in meshes) == sum(len(c) for c in helpers.tri_quad_mesh.cells)
    if op == "Partition":
        assert all("partition:original_point_id" in m.point_data for m in meshes)


def test_sanitized_key_collisions_fail_before_any_write(tmp_path, runner):
    mesh = helpers.tet_mesh.copy()
    mesh.regions = [mp.Region("a/b", "cell", [0]), mp.Region("a_b", "cell", [0])]
    source = _write(tmp_path / "input.vtu", mesh)
    with pytest.raises(Exception, match="collides"):
        runner({
            "Version": 2, "Input": {"Path": source},
            "Output": {"Pattern": str(tmp_path / "piece_{key}.vtu"), "Codec": "none"},
            "Operations": [{"Op": "Split", "By": "regions"}],
        })
    assert not list(tmp_path.glob("piece_*"))


@pytest.mark.parametrize("operations,output,error", [
    ([{"Op": "Merge", "Inputs": []}], {"Path": "out.vtu"}, "Inputs"),
    ([{"Op": "Interpolate", "Inputs": ["a", "b"]}], {"Path": "out.vtu"}, "Inputs"),
    ([{"Op": "Split"}, {"Op": "Quality"}], {"Pattern": "out_{key}.vtu"}, "terminal"),
    ([{"Op": "Split"}], {"Pattern": "out_{part}.vtu"}, "needs"),
    ([{"Op": "Split"}], {"Pattern": "out.vtu"}, "requires"),
    ([], {"Path": "out.vtu", "Pattern": "out_{key}.vtu"}, "mutually exclusive"),
])
def test_v2_schema_errors_precede_io(runner, operations, output, error):
    with pytest.raises(Exception, match=error):
        runner({"Version": 2, "Input": {"Path": "missing.vtu"}, "Output": output, "Operations": operations})
