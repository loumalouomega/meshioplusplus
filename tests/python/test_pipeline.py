"""The settings.json pipeline: the pure-Python runner (`run_pipeline`), its
parity with the C++ engine (`_core.run_pipeline_file`), the transcribed op/key
vocabulary table, both CLI verbs' flags, and the strict schema errors.
"""

import copy
import json
import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._pipeline import _EXCLUDED_OPS, _OP_TABLE

from . import helpers


@pytest.fixture
def settings_env(tmp_path):
    """A tet mesh with a point_data field on disk, plus in/out paths."""
    mesh = copy.deepcopy(helpers.tet_mesh)
    mesh.point_data["temperature"] = np.linspace(0.0, 1.0, len(mesh.points))
    in_path = tmp_path / "in.vtu"
    meshioplusplus.write(str(in_path), mesh)
    return {
        "in": str(in_path),
        "out": str(tmp_path / "out.vtu"),
        "tmp": tmp_path,
    }


def make_settings(env, operations):
    return {
        "Version": 1,
        "Input": {"Path": env["in"]},
        "Operations": operations,
        "Output": {"Path": env["out"]},
    }


# --------------------------------------------------------------------------- #
# The vocabulary table is a transcription -- pin it against the C++ owner.    #
# --------------------------------------------------------------------------- #
def test_op_table_matches_core():
    if not hasattr(_core, "pipeline_op_table"):
        pytest.skip("_core predates the pipeline")
    core_table = _core.pipeline_op_table()
    assert {k: list(v) for k, v in _OP_TABLE.items()} == core_table


def test_excluded_ops_error_by_name():
    for op, message in _EXCLUDED_OPS.items():
        with pytest.raises(ValueError, match="CLI verb"):
            meshioplusplus.run_pipeline(
                {
                    "Input": {"Path": "a"},
                    "Output": {"Path": "b"},
                    "Operations": [{"Op": op}],
                }
            )
        assert "CLI verb" in message


# --------------------------------------------------------------------------- #
# The runner itself                                                           #
# --------------------------------------------------------------------------- #
def test_run_pipeline_chains_and_reports(settings_env):
    report = meshioplusplus.run_pipeline(
        make_settings(
            settings_env,
            [
                {"Op": "Transform", "RotateAxis": [0, 0, 1], "RotateDegrees": 45},
                {"Op": "Gradient", "Array": "temperature"},
                {"Op": "Quality"},
                {"Op": "Clean"},
            ],
        )
    )
    assert [s["op"] for s in report["steps"]] == [
        "Transform",
        "Gradient",
        "Quality",
        "Clean",
    ]
    assert report["steps"][3]["PointsWelded"] == 0
    out = meshioplusplus.read(settings_env["out"])
    assert "temperature:gradient" in out.cell_data
    assert "quality:scaled_jacobian" in out.cell_data


def test_refine_record_hierarchy_step(settings_env):
    """RecordHierarchy round-trips through a pipeline document -- both the
    pure-Python runner and (when available) the C++ engine, since the op is
    dispatched generically off the same _OP_TABLE/pipeline_op_table pair
    test_op_table_matches_core pins."""
    settings = make_settings(settings_env, [{"Op": "Refine", "RecordHierarchy": True}])
    meshioplusplus.run_pipeline(settings)
    out = meshioplusplus.read(settings_env["out"])
    assert "refine:cell_id" in out.cell_data
    assert "refine:parent_id" in out.cell_data
    assert "refine:entity" in out.point_data  # the multigrid-stencil fix

    if hasattr(_core, "run_pipeline_json"):
        settings["Output"]["Path"] = str(settings_env["tmp"] / "out_cpp.vtu")
        _core.run_pipeline_json(json.dumps(settings))
        out_cpp = meshioplusplus.read(settings["Output"]["Path"])
        assert "refine:cell_id" in out_cpp.cell_data
        assert "refine:parent_id" in out_cpp.cell_data


def test_subdivide_pipeline_step(settings_env):
    """Subdivide dispatches generically off the same _OP_TABLE/
    pipeline_op_table pair test_op_table_matches_core pins. Since subdivide
    has no numpy fallback, both the Python runner and the C++ engine call the
    same underlying `_core.subdivide` -- so, unlike EstimateError, the two
    engines' output is bit-for-bit comparable, not merely within tolerance."""
    settings = make_settings(
        settings_env, [{"Op": "Subdivide", "RecordParentIds": True}]
    )
    report = meshioplusplus.run_pipeline(settings)
    assert report["steps"][0]["op"] == "Subdivide"
    out = meshioplusplus.read(settings_env["out"])
    # helpers.tet_mesh has two tetra cells, 4 faces each -> 8 children.
    assert len(out.cells[0].data) == 8
    assert "subdivide:parent_cell" in out.cell_data

    if hasattr(_core, "run_pipeline_json"):
        settings["Output"]["Path"] = str(settings_env["tmp"] / "out_cpp.vtu")
        _core.run_pipeline_json(json.dumps(settings))
        out_cpp = meshioplusplus.read(settings["Output"]["Path"])
        assert meshioplusplus.meshes_equal(out, out_cpp)


def test_agglomerate_pipeline_step(settings_env):
    """Agglomerate dispatches generically off the same _OP_TABLE/
    pipeline_op_table pair, and like Subdivide has no numpy fallback, so both
    engines call the same underlying `_core.agglomerate`."""
    settings = make_settings(
        settings_env, [{"Op": "Agglomerate", "TargetGroupSize": 2}]
    )
    report = meshioplusplus.run_pipeline(settings)
    assert report["steps"][0]["op"] == "Agglomerate"
    out = meshioplusplus.read(settings_env["out"])
    # helpers.tet_mesh has two face-adjacent tetra cells -> one merged cell.
    assert len(out.cells[0].data) == 1

    if hasattr(_core, "run_pipeline_json"):
        settings["Output"]["Path"] = str(settings_env["tmp"] / "out_cpp.vtu")
        _core.run_pipeline_json(json.dumps(settings))
        out_cpp = meshioplusplus.read(settings["Output"]["Path"])
        assert meshioplusplus.meshes_equal(out, out_cpp)


def test_decimate_volume_pipeline_step(settings_env):
    """DecimateVolume dispatches generically off the same _OP_TABLE/
    pipeline_op_table pair, and like Subdivide/Agglomerate has no numpy
    fallback, so both engines call the same underlying
    `_core.decimate_volume`."""
    settings = make_settings(
        settings_env,
        [{"Op": "DecimateVolume", "TargetCells": 1, "PreserveFeatures": False}],
    )
    report = meshioplusplus.run_pipeline(settings)
    assert report["steps"][0]["op"] == "DecimateVolume"
    out = meshioplusplus.read(settings_env["out"])
    # helpers.tet_mesh has two face-adjacent tetra cells -> collapses to one.
    assert len(out.cells[0].data) <= 1

    if hasattr(_core, "run_pipeline_json"):
        settings["Output"]["Path"] = str(settings_env["tmp"] / "out_cpp.vtu")
        _core.run_pipeline_json(json.dumps(settings))
        out_cpp = meshioplusplus.read(settings["Output"]["Path"])
        assert meshioplusplus.meshes_equal(out, out_cpp)


def test_estimate_error_pipeline_step(settings_env):
    """EstimateError dispatches generically off the same _OP_TABLE/
    pipeline_op_table pair test_op_table_matches_core pins. Unlike
    test_refine_record_hierarchy_step this does NOT assert cross-engine
    equality of the produced arrays -- estimate_error's own
    test_cpp_matches_python (test_error.py) already established that the two
    engines only agree to within a numeric tolerance (the same,
    already-accepted precedent data_average's own measure weighting has), not
    bit-for-bit, so this only checks that both engines run and attach the
    right arrays."""
    settings = make_settings(
        settings_env,
        [
            {
                "Op": "EstimateError",
                "Array": "temperature",
                "Marking": "absolute",
                "MarkingValue": 1e-9,
            }
        ],
    )
    report = meshioplusplus.run_pipeline(settings)
    assert report["steps"][0]["op"] == "EstimateError"
    assert "GlobalError" in report["steps"][0]
    out = meshioplusplus.read(settings_env["out"])
    assert "error:zz" in out.cell_data
    assert "error:marked" in out.cell_data

    if hasattr(_core, "run_pipeline_json"):
        settings["Output"]["Path"] = str(settings_env["tmp"] / "out_cpp.vtu")
        cpp_report = _core.run_pipeline_json(json.dumps(settings))
        assert cpp_report["steps"][0]["op"] == "EstimateError"
        out_cpp = meshioplusplus.read(settings["Output"]["Path"])
        assert "error:zz" in out_cpp.cell_data
        assert "error:marked" in out_cpp.cell_data


def test_hessian_pipeline_step(settings_env):
    """Hessian dispatches generically off the same _OP_TABLE/
    pipeline_op_table pair test_op_table_matches_core pins. Unlike
    EstimateError, its Python fallback (`_hessian_py`) composes the public
    `gradient()` function directly -- which itself prefers `_core.gradient`
    when available -- so with `_core` installed both pipeline engines
    literally call the same underlying gradient() twice, and the two
    engines' output is bit-for-bit comparable, not merely within
    tolerance (test_hessian.py::test_cpp_matches_python established this)."""
    settings = make_settings(settings_env, [{"Op": "Hessian", "Array": "temperature"}])
    report = meshioplusplus.run_pipeline(settings)
    assert report["steps"][0]["op"] == "Hessian"
    assert "NumSkipped" in report["steps"][0]
    out = meshioplusplus.read(settings_env["out"])
    assert "temperature:hessian" in out.cell_data

    if hasattr(_core, "run_pipeline_json"):
        settings["Output"]["Path"] = str(settings_env["tmp"] / "out_cpp.vtu")
        cpp_report = _core.run_pipeline_json(json.dumps(settings))
        assert cpp_report["steps"][0]["op"] == "Hessian"
        out_cpp = meshioplusplus.read(settings["Output"]["Path"])
        assert meshioplusplus.meshes_equal(out, out_cpp)


def test_run_pipeline_accepts_path_text_and_dict(settings_env, tmp_path):
    settings = make_settings(settings_env, [{"Op": "Quality"}])
    # dict
    meshioplusplus.run_pipeline(settings)
    # JSON text
    meshioplusplus.run_pipeline(json.dumps(settings))
    # file path (str and PathLike)
    path = tmp_path / "settings.json"
    path.write_text(json.dumps(settings))
    meshioplusplus.run_pipeline(str(path))
    meshioplusplus.run_pipeline(path)


def test_path_overrides(settings_env, tmp_path):
    settings = make_settings(settings_env, [])
    other_out = tmp_path / "other.vtu"
    meshioplusplus.run_pipeline(settings, output_path=str(other_out))
    assert other_out.exists()
    # The override wins over a bogus path in the document too.
    settings["Input"]["Path"] = "/no/such/input.vtu"
    meshioplusplus.run_pipeline(
        settings, input_path=settings_env["in"], output_path=str(other_out)
    )


def test_isosurface_and_slice_ops(settings_env):
    report = meshioplusplus.run_pipeline(
        make_settings(
            settings_env,
            [{"Op": "Isosurface", "Array": "temperature", "Isovalue": 0.5}],
        )
    )
    assert report["steps"][0]["ContourCells"] >= 1
    report = meshioplusplus.run_pipeline(
        make_settings(
            settings_env,
            [{"Op": "Slice", "Point": [0, 0, 100], "Normal": [0, 0, 1]}],
        )
    )
    assert report["steps"][0]["SectionFaces"] == 0
    assert any("section is empty" in w for w in report["warnings"])


def test_section_is_an_alias_of_slice(settings_env):
    report = meshioplusplus.run_pipeline(
        make_settings(
            settings_env,
            [{"Op": "Section", "Point": [0, 0, 0.25], "Normal": [0, 0, 1]}],
        )
    )
    assert report["steps"][0]["op"] == "Section"
    assert report["steps"][0]["SectionFaces"] > 0


def test_partition_attaches_labels(settings_env):
    meshioplusplus.run_pipeline(
        make_settings(settings_env, [{"Op": "Partition", "Nparts": 2}])
    )
    out = meshioplusplus.read(settings_env["out"])
    assert "partition:part" in out.cell_data


def test_data_ops(settings_env):
    meshioplusplus.run_pipeline(
        make_settings(
            settings_env,
            [
                {"Op": "DataCalc", "Expr": "t2 = temperature * 2"},
                {"Op": "DataRename", "Point": ["t2:doubled"]},
                {"Op": "DataCondition", "Mode": "normalize", "Names": ["doubled"]},
                {"Op": "ToCell", "Names": ["doubled"]},
            ],
        )
    )
    out = meshioplusplus.read(settings_env["out"])
    assert "doubled" in out.point_data
    assert "doubled" in out.cell_data


def test_output_encoding_and_codec(settings_env):
    settings = make_settings(settings_env, [])
    settings["Output"]["Encoding"] = "binary"
    settings["Output"]["Codec"] = "zlib"
    meshioplusplus.run_pipeline(settings)
    assert pathlib.Path(settings_env["out"]).exists()
    # A codec on a format with no block codec is an error, never ignored.
    settings["Output"]["Path"] = str(settings_env["tmp"] / "out.msh")
    settings["Output"]["Codec"] = "zstd"
    with pytest.raises(ValueError):
        meshioplusplus.run_pipeline(settings)


def test_transform_requires_exactly_one_source(settings_env):
    with pytest.raises(ValueError, match="exactly one"):
        meshioplusplus.run_pipeline(
            make_settings(
                settings_env,
                [{"Op": "Transform", "Translate": [1, 0, 0], "Scale": 2}],
            )
        )
    with pytest.raises(ValueError, match="exactly one"):
        meshioplusplus.run_pipeline(make_settings(settings_env, [{"Op": "Transform"}]))


# --------------------------------------------------------------------------- #
# Strict schema errors, always naming the offender                            #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "mutate, match",
    [
        (lambda s: s.update(Bogus=1), "unknown key 'Bogus'"),
        (lambda s: s.update(Version=2), "Version"),
        (lambda s: s.pop("Input"), "Input is required"),
        (lambda s: s["Output"].pop("Path"), "Output.Path is required"),
        (
            lambda s: s.update(Operations=[{"Op": "Nope"}]),
            "unknown operation 'Nope'",
        ),
        (
            lambda s: s.update(Operations=[{"Op": "Clean", "Foo": 1}]),
            "unknown parameter 'Foo'",
        ),
        (
            lambda s: s.update(Operations=[{"Op": "Clean", "Atol": "tight"}]),
            "'Atol' must be a number",
        ),
        (
            lambda s: s.update(Operations=[{"Op": "Clean", "Weld": 1}]),
            "'Weld' must be a boolean",
        ),
        (lambda s: s["Output"].update(Codec="brotli"), "Codec"),
        (lambda s: s["Input"].update(Options={"Mmap": "maybe"}), "Mmap"),
    ],
)
def test_schema_errors(settings_env, mutate, match):
    settings = make_settings(settings_env, [])
    mutate(settings)
    with pytest.raises(ValueError, match=match):
        meshioplusplus.run_pipeline(settings)


def test_lenient_is_a_named_warning_in_the_python_runner(settings_env):
    settings = make_settings(settings_env, [])
    settings["Input"]["Options"] = {"Lenient": True}
    report = meshioplusplus.run_pipeline(settings)
    assert any("Lenient" in w for w in report["warnings"])


# --------------------------------------------------------------------------- #
# C++ engine parity                                                           #
# --------------------------------------------------------------------------- #
needs_core_json = pytest.mark.skipif(
    not getattr(_core, "__has_json__", False),
    reason="_core built without MESHIOPLUSPLUS_WITH_JSON",
)
# The two tests below drive _core's own read/write, so the registry-default
# binary+zlib VTU writer must be compiled in too (Windows CI builds _core
# with native paths, including zlib, off) -- not just the JSON parser.
needs_core_json_and_zlib = pytest.mark.skipif(
    not (
        getattr(_core, "__has_json__", False) and getattr(_core, "__has_zlib__", False)
    ),
    reason="_core built without MESHIOPLUSPLUS_WITH_JSON or MESHIOPLUSPLUS_WITH_ZLIB",
)


@needs_core_json_and_zlib
def test_cpp_matches_python(settings_env, tmp_path):
    operations = [
        {"Op": "ConvertCells", "Mode": "simplexify"},
        {"Op": "Quality"},
        {"Op": "Gradient", "Array": "temperature"},
        {"Op": "Clean"},
    ]
    settings = make_settings(settings_env, operations)
    py_report = meshioplusplus.run_pipeline(settings)
    py_mesh = meshioplusplus.read(settings_env["out"])

    cpp_out = tmp_path / "cpp_out.vtu"
    settings["Output"]["Path"] = str(cpp_out)
    settings_path = tmp_path / "settings.json"
    settings_path.write_text(json.dumps(settings))
    cpp_report = _core.run_pipeline_file(str(settings_path))
    cpp_mesh = meshioplusplus.read(str(cpp_out))

    assert meshioplusplus.meshes_equal(py_mesh, cpp_mesh)
    assert py_report["steps"] == cpp_report["steps"]
    assert py_report["warnings"] == cpp_report["warnings"]


@needs_core_json_and_zlib
def test_core_run_pipeline_json(settings_env):
    settings = make_settings(settings_env, [{"Op": "Quality"}])
    report = _core.run_pipeline_json(json.dumps(settings))
    assert report["steps"] == [{"op": "Quality"}]
    assert (
        "quality:scaled_jacobian" in meshioplusplus.read(settings_env["out"]).cell_data
    )


def test_core_throws_by_name_when_compiled_out(settings_env):
    if getattr(_core, "__has_json__", True):
        pytest.skip("_core carries the JSON parser")
    with pytest.raises(Exception, match="MESHIOPLUSPLUS_WITH_JSON"):
        _core.run_pipeline_json("{}")


# --------------------------------------------------------------------------- #
# The CLI verb                                                                #
# --------------------------------------------------------------------------- #
def test_cli_pipeline_verb(settings_env, tmp_path, capsys):
    from meshioplusplus._cli import main

    settings_path = tmp_path / "settings.json"
    settings_path.write_text(
        json.dumps(make_settings(settings_env, [{"Op": "Quality"}]))
    )
    assert main(["pipeline", str(settings_path)]) == 0
    out = capsys.readouterr().out
    assert "step 1: Quality" in out

    other = tmp_path / "cli_out.vtu"
    assert main(["pipeline", str(settings_path), "--output", str(other), "--json"]) == 0
    assert other.exists()
    assert json.loads(capsys.readouterr().out)["steps"] == [{"op": "Quality"}]


def test_compute_sdf_step_replaces_the_geometry(tmp_path):
    """``ComputeSdf`` is ``Voxelize``'s sibling in the single-mesh chain.

    Both take a mesh in and hand a *different* mesh back rather than
    transforming the input's geometry, which is exactly the shape a chain wants.
    """
    surface = meshioplusplus.convert_cells(
        meshioplusplus.extract_surface(meshioplusplus.grid([2, 2, 2])),
        mode="simplexify",
    )
    src = tmp_path / "in.vtu"
    out = tmp_path / "out.vtu"
    meshioplusplus.write(src, surface)
    report = meshioplusplus.run_pipeline(
        {
            "Version": 1,
            "Input": {"Path": str(src)},
            "Operations": [{"Op": "ComputeSdf", "Resolution": [4, 4, 4]}],
            "Output": {"Path": str(out)},
        }
    )
    assert report["steps"][0]["op"] == "ComputeSdf"
    assert report["steps"][0]["MaxDepth"] == 0
    back = meshioplusplus.read(out)
    assert back.cells[0].type == "hexahedron"
    assert len(back.cells[0].data) == 64
    assert "sdf:distance" in back.point_data


def test_crop_step_takes_a_data_predicate(tmp_path):
    mesh = meshioplusplus.grid([4, 4, 4])
    mesh.cell_data["t"] = [np.arange(64, dtype=np.float64)]
    src = tmp_path / "in.vtu"
    out = tmp_path / "out.vtu"
    meshioplusplus.write(src, mesh)
    report = meshioplusplus.run_pipeline(
        {
            "Version": 1,
            "Input": {"Path": str(src)},
            "Operations": [{"Op": "Crop", "Where": "t", "Compare": "<", "Value": 10.0}],
            "Output": {"Path": str(out)},
        }
    )
    assert report["steps"][0]["CellsKept"] == 10.0

    # `Mode` means nothing for a per-cell predicate, so it is refused rather
    # than ignored.
    with pytest.raises(ValueError, match="Mode"):
        meshioplusplus.run_pipeline(
            {
                "Version": 1,
                "Input": {"Path": str(src)},
                "Operations": [
                    {"Op": "Crop", "Where": "t", "Value": 10.0, "Mode": "any"}
                ],
                "Output": {"Path": str(out)},
            }
        )
