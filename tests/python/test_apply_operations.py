"""`apply_operations`: the settings pipeline applied to a mesh in memory
(roadmap 7.2.1) -- same vocabulary, validation and report as `run_pipeline`,
without the Input read and the Output write."""

import copy
import json
import pathlib

import numpy as np
import pytest

import meshioplusplus

from . import helpers

CAD = pathlib.Path(__file__).parent / "meshes" / "cad_json"


def _tets():
    mesh = copy.deepcopy(helpers.tet_mesh)
    mesh.point_data["temperature"] = np.linspace(0.0, 1.0, len(mesh.points))
    return mesh


def _cells(mesh):
    return sum(len(block.data) for block in mesh.cells)


def _snapshot(mesh):
    return (
        mesh.points.copy(),
        [(b.type, np.array(b.data, copy=True)) for b in mesh.cells],
        {k: np.array(v, copy=True) for k, v in mesh.point_data.items()},
        {k: [np.array(a, copy=True) for a in v] for k, v in mesh.cell_data.items()},
    )


def _assert_same_snapshot(a, b):
    np.testing.assert_array_equal(a[0], b[0])
    assert [t for t, _ in a[1]] == [t for t, _ in b[1]]
    for (_, x), (_, y) in zip(a[1], b[1]):
        np.testing.assert_array_equal(x, y)
    assert a[2].keys() == b[2].keys()
    for k in a[2]:
        np.testing.assert_array_equal(a[2][k], b[2][k])
    assert a[3].keys() == b[3].keys()
    for k in a[3]:
        for x, y in zip(a[3][k], b[3][k]):
            np.testing.assert_array_equal(x, y)


CHAIN = [
    {"Op": "Transform", "RotateAxis": [0, 0, 1], "RotateDegrees": 45},
    {"Op": "Gradient", "Array": "temperature"},
    {"Op": "Quality"},
    {"Op": "DataCalc", "Expr": "t2 = temperature * 2"},
    {"Op": "DataRename", "Point": ["temperature:temp"]},
    {"Op": "Clean"},
]


def test_refine_multiplies_tets_by_eight():
    """The probe of the ParaView Pipeline filter (roadmap 6.8)."""
    mesh = _tets()
    out, report = meshioplusplus.apply_operations(mesh, [{"Op": "Refine"}])
    assert _cells(out) == 8 * _cells(mesh)
    assert [s["op"] for s in report["steps"]] == ["Refine"]
    assert report["warnings"] == []


def test_matches_run_pipeline(tmp_path):
    mesh = _tets()
    in_path, out_path = tmp_path / "in.vtu", tmp_path / "out.vtu"
    meshioplusplus.write(str(in_path), mesh, compression=None)
    file_report = meshioplusplus.run_pipeline(
        {
            "Version": 1,
            "Input": {"Path": str(in_path)},
            "Operations": CHAIN,
            "Output": {"Path": str(out_path)},
        }
    )
    from_file = meshioplusplus.read(str(out_path))

    out, report = meshioplusplus.apply_operations(mesh, CHAIN)
    assert report == file_report
    np.testing.assert_allclose(out.points, from_file.points, atol=1e-12)
    assert [b.type for b in out.cells] == [b.type for b in from_file.cells]
    for a, b in zip(out.cells, from_file.cells):
        np.testing.assert_array_equal(a.data, b.data)
    assert set(out.point_data) == set(from_file.point_data)
    assert set(out.cell_data) == set(from_file.cell_data)
    assert "temp" in out.point_data and "t2" in out.point_data
    assert "temperature:gradient" in out.cell_data


def test_input_mesh_is_not_modified():
    mesh = _tets()
    before = _snapshot(mesh)
    meshioplusplus.apply_operations(
        mesh,
        CHAIN
        + [
            {"Op": "DataDrop", "Point": ["temp"]},
            {"Op": "Refine"},
            {"Op": "Smooth", "Iterations": 2},
        ],
    )
    _assert_same_snapshot(before, _snapshot(mesh))


def test_operations_as_json_text_and_settings_document():
    mesh = _tets()
    expected, _ = meshioplusplus.apply_operations(mesh, [{"Op": "Refine"}])
    from_text, _ = meshioplusplus.apply_operations(mesh, json.dumps([{"Op": "Refine"}]))
    from_doc, report = meshioplusplus.apply_operations(
        mesh,
        {
            "Version": 1,
            "Input": {"Path": "ignored.vtu"},
            "Operations": [{"Op": "Refine"}],
            "Output": {"Path": "ignored.vtu"},
        },
    )
    assert _cells(from_text) == _cells(expected) == _cells(from_doc)
    assert [s["op"] for s in report["steps"]] == ["Refine"]


def test_empty_operations_return_the_mesh_unchanged():
    mesh = _tets()
    out, report = meshioplusplus.apply_operations(mesh, [])
    assert report == {"steps": [], "warnings": []}
    np.testing.assert_array_equal(out.points, mesh.points)


def test_schema_errors_raise_before_any_step_runs(monkeypatch):
    from meshioplusplus import _pipeline

    calls = []
    real = _pipeline._apply_step

    def spy(*args, **kwargs):
        calls.append(args[1]["Op"])
        return real(*args, **kwargs)

    monkeypatch.setattr(_pipeline, "_apply_step", spy)
    mesh = _tets()
    with pytest.raises(ValueError, match="unknown operation 'Frobnicate'"):
        meshioplusplus.apply_operations(mesh, [{"Op": "Refine"}, {"Op": "Frobnicate"}])
    with pytest.raises(ValueError, match="unknown parameter 'Bogus'"):
        meshioplusplus.apply_operations(
            mesh, [{"Op": "Quality"}, {"Op": "Clean", "Bogus": 1}]
        )
    with pytest.raises(ValueError, match="merge"):
        meshioplusplus.apply_operations(mesh, [{"Op": "Merge"}])
    assert calls == []


def test_bad_arguments():
    mesh = _tets()
    with pytest.raises(ValueError, match="Version"):
        meshioplusplus.apply_operations(mesh, [], version=3)
    with pytest.raises(ValueError, match="unknown key 'Bogus'|Bogus"):
        meshioplusplus.apply_operations(mesh, {"Operations": [], "Bogus": 1})
    with pytest.raises(ValueError, match="list of steps"):
        meshioplusplus.apply_operations(mesh, 42)
    with pytest.raises(ValueError, match="Generate"):
        meshioplusplus.apply_operations(None, [{"Op": "Quality"}])
    with pytest.raises(TypeError, match="Mesh"):
        meshioplusplus.apply_operations("mesh.vtu", [{"Op": "Quality"}])


def test_fan_out_steps_are_refused_by_name():
    mesh = _tets()
    with pytest.raises(ValueError, match="'Split'.*run_pipeline"):
        meshioplusplus.apply_operations(
            mesh, [{"Op": "Split", "By": "block"}], version=2
        )
    with pytest.raises(ValueError, match="'Partition'.*run_pipeline"):
        meshioplusplus.apply_operations(
            mesh,
            [{"Op": "Partition", "Nparts": 2, "GhostLayers": 1}],
            version=2,
        )


def test_generate_starts_from_none():
    mesh, report = meshioplusplus.apply_operations(
        None,
        [
            {"Op": "Generate", "Shape": "box", "Dims": [2, 2, 2], "CellType": "tetra"},
            {"Op": "Quality"},
        ],
    )
    assert _cells(mesh) > 0
    assert "quality:scaled_jacobian" in mesh.cell_data
    assert [s["op"] for s in report["steps"]] == ["Generate", "Quality"]


def test_model_step_without_a_model_raises_by_name():
    with pytest.raises(ValueError, match="MeshCadSurface.*mesh.info"):
        meshioplusplus.apply_operations(
            _tets(), [{"Op": "MeshCadSurface", "TargetSize": 0.3}]
        )


def test_model_step_uses_the_model_in_mesh_info():
    """A CAD read carries its B-rep in ``mesh.info``; the in-memory run meshes
    it exactly as the file run does."""
    from meshioplusplus import _core

    if not hasattr(_core, "mesh_cad_surface"):
        pytest.skip("core predates mesh_cad_surface")
    path = CAD / "sphere.cad.json"
    model = meshioplusplus.read(str(path))
    ops = [{"Op": "MeshCadSurface", "TargetSize": 0.3}, {"Op": "Quality"}]
    out, report = meshioplusplus.apply_operations(model, ops)
    assert [s["op"] for s in report["steps"]] == ["MeshCadSurface", "Quality"]
    assert _cells(out) > 0
    # A mesh-changing step drops the model, as in a file run.
    with pytest.raises(ValueError, match="drops the model"):
        meshioplusplus.apply_operations(
            model, [{"Op": "Quality"}, {"Op": "MeshCadSurface", "TargetSize": 0.3}]
        )


def test_version_two_merge_reads_its_inputs(tmp_path):
    extra = tmp_path / "extra.vtu"
    other = _tets()
    other.points = other.points + 10.0
    meshioplusplus.write(str(extra), other, compression=None)
    mesh = _tets()
    out, report = meshioplusplus.apply_operations(
        mesh, [{"Op": "Merge", "Inputs": [str(extra)]}], version=2
    )
    assert len(out.points) == 2 * len(mesh.points)
    assert report["steps"][0]["NumInputs"] == 2.0


def test_exported_and_exempt_from_mcp():
    assert "apply_operations" in meshioplusplus.__all__
    assert meshioplusplus.apply_operations is not None
