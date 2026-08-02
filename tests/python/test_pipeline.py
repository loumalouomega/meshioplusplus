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
