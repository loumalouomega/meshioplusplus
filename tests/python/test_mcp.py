"""Tests for the MCP server (meshioplusplus.mcp).

Two halves, mirroring test_interop.py's split:

* The **pure half** exercises the payload layer (``meshioplusplus.mcp._tools``)
  which imports no third-party library beyond numpy — it runs in the default
  CI matrix with the ``mcp`` SDK absent. Every report is asserted to be
  strict JSON (``json.dumps(..., allow_nan=False)``).
* The **gated half** (``pytest.importorskip("mcp")``) builds the FastMCP
  server and exercises tool listing, tool calls, resources and the error
  payload shape.

The parity guard (`test_every_operation_has_a_tool`) is the enforcement
mechanism behind CLAUDE.md's "keep the MCP server in sync" rule: a new public
operation in ``meshioplusplus.__all__`` fails here until it is claimed by a
tool's ``wraps`` (or consciously exempted in ``_NOT_TOOLS``).
"""

import json
import os

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus.mcp import TOOL_REGISTRY, _tools

# --------------------------------------------------------------------------- #
# Fixtures                                                                    #
# --------------------------------------------------------------------------- #


@pytest.fixture(autouse=True)
def _reset_root():
    yield
    _tools.set_root(None)


def _mixed_mesh():
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1]], dtype=float
    )
    cells = [
        ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
        ("tetra", np.array([[0, 1, 2, 4]])),
    ]
    return meshioplusplus.Mesh(
        points,
        cells,
        point_data={"t": np.linspace(0.0, 1.0, 5)},
        cell_data={"c": [np.array([1.0, 2.0]), np.array([3.0])]},
    )


@pytest.fixture()
def mesh_file(tmp_path):
    path = str(tmp_path / "in.vtu")
    meshioplusplus.write(path, _mixed_mesh())
    return path


def _dump(report):
    """Assert the report is strict JSON and return it."""
    json.dumps(report, allow_nan=False)
    return report


# --------------------------------------------------------------------------- #
# Pure half: sanitizer and sandbox                                            #
# --------------------------------------------------------------------------- #


def test_json_safe_numpy_and_non_finite():
    out = _dump(
        _tools._json_safe(
            {
                "i": np.int64(3),
                "f": np.float32(1.5),
                "nan": float("nan"),
                "inf": np.float64("inf"),
                "arr": np.array([1.0, float("-inf")]),
                "nested": {"t": (np.int32(1), [np.float64(2.0)])},
            }
        )
    )
    assert out["i"] == 3 and out["f"] == 1.5
    assert out["nan"] is None and out["inf"] is None
    assert out["arr"] == [1.0, None]
    assert out["nested"]["t"] == [1, [2.0]]
    assert out["non_finite_replaced"] == 3


def test_json_safe_truncates_large_arrays():
    out = _dump(_tools._json_safe({"big": np.arange(5000).reshape(100, 50)}))
    wrapper = out["big"]
    assert wrapper["truncated"] is True
    assert wrapper["size"] == 5000
    assert wrapper["shape"] == [100, 50]
    assert wrapper["preview"] == list(range(50))


def test_root_confinement(tmp_path, mesh_file):
    _tools.set_root(str(tmp_path))
    # relative paths resolve against the root
    assert _tools.tool_info("in.vtu")["num_points"] == 5
    # absolute paths inside the root are fine
    assert _tools.tool_info(mesh_file)["num_points"] == 5
    with pytest.raises(ValueError, match="outside the configured root"):
        _tools.tool_info("/etc/passwd")
    with pytest.raises(ValueError, match="outside the configured root"):
        _tools.tool_info("../escape.vtu")


def test_missing_input_is_a_clean_error(tmp_path):
    with pytest.raises(ValueError, match="input file not found"):
        _tools.tool_info(str(tmp_path / "nope.vtu"))


# --------------------------------------------------------------------------- #
# Pure half: inspection tools                                                 #
# --------------------------------------------------------------------------- #


def test_formats_payload():
    out = _dump(_tools.tool_formats())
    assert "vtu" in out["readable"] and "vtu" in out["writable"]
    assert out["extensions"][".vtu"] == ["vtu"]


def test_sniff(mesh_file):
    out = _dump(_tools.tool_sniff(mesh_file))
    assert out["format"] == "vtu"
    assert out["from_extension"] == ["vtu"]


def test_info(mesh_file):
    out = _dump(_tools.tool_info(mesh_file))
    assert out["num_points"] == 5
    assert out["num_cells"] == 3
    assert out["point_data_names"] == ["t"]


def test_stats(mesh_file):
    out = _dump(_tools.tool_stats(mesh_file))
    assert out["num_points"] == 5
    assert out["cell_type_counts"] == {"triangle": 2, "tetra": 1}


def test_quality_summary_and_annotated_output(mesh_file, tmp_path):
    out = _dump(_tools.tool_quality(mesh_file))
    assert out["num_cells"] == 3
    assert any(m.startswith("quality:") for m in out["metrics"])
    assert isinstance(out["cell_arrays"], str)  # omitted, with a pointer
    annotated = str(tmp_path / "annotated.vtu")
    out = _dump(_tools.tool_quality(mesh_file, output_path=annotated))
    assert os.path.isfile(out["output_path"])
    reread = meshioplusplus.read(annotated)
    assert any(k.startswith("quality:") for k in reread.cell_data)


def test_data_info(mesh_file):
    out = _dump(_tools.tool_data_info(mesh_file))
    names = {(a["location"], a["name"]) for a in out["arrays"]}
    assert ("point_data", "t") in names and ("cell_data", "c") in names


def test_bandwidth(mesh_file):
    out = _dump(_tools.tool_bandwidth(mesh_file))
    assert out["bandwidth"] == 4


def test_data_preview(mesh_file):
    out = _dump(_tools.tool_data_preview(mesh_file, "t", location="point", limit=3))
    assert out["num_rows"] == 5
    assert out["values"] == [0.0, 0.25, 0.5]
    out = _dump(_tools.tool_data_preview(mesh_file, "c", location="cell"))
    assert out["values"] == [[1.0], [2.0], [3.0]]  # block-major concatenation
    with pytest.raises(ValueError, match="no point data array named 'zz'"):
        _tools.tool_data_preview(mesh_file, "zz")
    with pytest.raises(ValueError, match="unknown location"):
        _tools.tool_data_preview(mesh_file, "t", location="vertex")


def test_diff(mesh_file, tmp_path):
    same = str(tmp_path / "same.vtu")
    meshioplusplus.write(same, _mixed_mesh())
    out = _dump(_tools.tool_diff(mesh_file, same))
    assert out["equal"] is True
    moved = _mixed_mesh()
    moved.points = moved.points + 1.0
    other = str(tmp_path / "other.vtu")
    meshioplusplus.write(other, moved)
    out = _dump(_tools.tool_diff(mesh_file, other))
    assert out["equal"] is False and out["verdict"] == "different"


# --------------------------------------------------------------------------- #
# Pure half: conversion                                                       #
# --------------------------------------------------------------------------- #


def test_convert_roundtrip(mesh_file, tmp_path):
    out_path = str(tmp_path / "out.vtk")
    out = _dump(_tools.tool_convert(mesh_file, out_path))
    assert out["output_format"] == "vtk"
    assert os.path.isfile(out["output_path"])
    assert out["num_points"] == 5 and out["num_cells"] == 3


def test_convert_ascii_variant(mesh_file, tmp_path):
    out = _dump(_tools.tool_convert(mesh_file, str(tmp_path / "a.vtu"), mode="ascii"))
    with open(out["output_path"], "rb") as f:
        assert b"ascii" in f.read()


def test_convert_variant_errors(mesh_file, tmp_path):
    with pytest.raises(ValueError, match="has no ascii variant"):
        _tools.tool_convert(mesh_file, str(tmp_path / "a.obj"), mode="ascii")
    with pytest.raises(ValueError, match="only vti/vtu/vtp"):
        _tools.tool_convert(mesh_file, str(tmp_path / "a.vtk"), compression="zstd")
    with pytest.raises(ValueError, match="unknown mode"):
        _tools.tool_convert(mesh_file, str(tmp_path / "a.vtu"), mode="fast")


def test_pipeline_tool(mesh_file, tmp_path):
    out_path = str(tmp_path / "piped.vtu")
    settings_path = str(tmp_path / "settings.json")
    with open(settings_path, "w") as f:
        json.dump(
            {
                "Input": {"Path": mesh_file},
                "Operations": [{"Op": "Quality"}, {"Op": "Clean"}],
                "Output": {"Path": out_path},
            },
            f,
        )
    out = _dump(_tools.tool_pipeline(settings_path))
    assert os.path.isfile(out["output_path"])
    assert [s["op"] for s in out["steps"]] == ["Quality", "Clean"]
    # overrides win over the document's paths
    other = str(tmp_path / "other.vtu")
    _dump(_tools.tool_pipeline(settings_path, output_path=other))
    assert os.path.isfile(other)


def test_pipeline_tool_sandboxes_the_inner_paths(mesh_file, tmp_path, monkeypatch):
    # A settings document naming a path outside the root must fail exactly the
    # way a path argument would -- the sandbox covers the document's insides.
    settings_path = str(tmp_path / "settings.json")
    with open(settings_path, "w") as f:
        json.dump(
            {
                "Input": {"Path": "/etc/passwd"},
                "Operations": [],
                "Output": {"Path": str(tmp_path / "o.vtu")},
            },
            f,
        )
    monkeypatch.setattr(_tools, "_ROOT", str(tmp_path))
    with pytest.raises(ValueError, match="outside the configured root"):
        _tools.tool_pipeline(settings_path)


# --------------------------------------------------------------------------- #
# Pure half: mesh operations                                                  #
# --------------------------------------------------------------------------- #


def test_clean_report(mesh_file, tmp_path):
    out = _dump(_tools.tool_clean(mesh_file, str(tmp_path / "clean.vtu")))
    for key in (
        "points_welded",
        "points_removed_orphan",
        "cells_dropped_degenerate",
        "cells_dropped_duplicate",
    ):
        assert key in out


def test_reorder_reports_bandwidth(mesh_file, tmp_path):
    out = _dump(_tools.tool_reorder(mesh_file, str(tmp_path / "ro.vtu")))
    assert out["bandwidth_after"] <= out["bandwidth_before"]


def test_slice(mesh_file, tmp_path):
    out = _dump(
        _tools.tool_slice(
            mesh_file, str(tmp_path / "s.vtu"), origin=[0, 0, 0.5], normal=[0, 0, 1]
        )
    )
    assert out["num_cells"] >= 1
    assert all(b["type"] in ("triangle", "quad") for b in out["cell_blocks"])


def test_split_writes_one_file_per_piece(mesh_file, tmp_path):
    out_dir = str(tmp_path / "pieces")
    out = _dump(_tools.tool_split(mesh_file, by="type", output_dir=out_dir))
    assert sorted(out["pieces"]) == ["tetra", "triangle"]
    for key, path in out["pieces"].items():
        assert os.path.isfile(path)
        assert os.path.basename(path) == f"in_{key}.vtu"
        assert out["summaries"][key]["num_cells"] >= 1


def test_merge(mesh_file, tmp_path):
    out = _dump(_tools.tool_merge([mesh_file, mesh_file], str(tmp_path / "m.vtu")))
    assert out["num_points"] == 10 and out["num_cells"] == 6
    with pytest.raises(ValueError, match="at least two"):
        _tools.tool_merge([mesh_file], str(tmp_path / "m2.vtu"))


def test_partition_writes_nparts_files(mesh_file, tmp_path):
    out_dir = str(tmp_path / "parts")
    out = _dump(_tools.tool_partition(mesh_file, 2, method="sfc", output_dir=out_dir))
    assert len(out["parts"]) == 2
    assert all(os.path.isfile(p) for p in out["parts"])
    assert sum(s["num_cells"] for s in out["summaries"]) == 3


def test_crop(mesh_file, tmp_path):
    out = _dump(
        _tools.tool_crop(
            mesh_file, str(tmp_path / "c.vtu"), bbox=[0, 0, 0, 1, 1, 0.5], mode="any"
        )
    )
    assert out["num_cells"] == 3
    with pytest.raises(ValueError, match="both plane_origin and plane_normal"):
        _tools.tool_crop(mesh_file, str(tmp_path / "c2.vtu"), plane_origin=[0, 0, 0])


def test_transform(mesh_file, tmp_path):
    out_path = str(tmp_path / "t.vtu")
    _dump(_tools.tool_transform(mesh_file, out_path, translate=[10.0, 0.0, 0.0]))
    moved = meshioplusplus.read(out_path)
    assert moved.points[:, 0].min() >= 10.0
    with pytest.raises(ValueError, match="both rotate_axis and rotate_degrees"):
        _tools.tool_transform(mesh_file, out_path, rotate_axis=[0, 0, 1])


# --------------------------------------------------------------------------- #
# Pure half: data operations                                                  #
# --------------------------------------------------------------------------- #


def test_data_calc_cli_spelling(mesh_file, tmp_path):
    out_path = str(tmp_path / "calc.vtu")
    out = _dump(_tools.tool_data_calc(mesh_file, out_path, "t2 = t * 2"))
    assert "t2" in out["point_data"]
    reread = meshioplusplus.read(out_path)
    assert np.allclose(reread.point_data["t2"], 2 * reread.point_data["t"])


def test_data_manage_rename(mesh_file, tmp_path):
    out = _dump(
        _tools.tool_data_manage(
            mesh_file, str(tmp_path / "dm.vtu"), rename=[["point", "t", "temp"]]
        )
    )
    assert out["point_data"] == ["temp"]
    assert out["renamed"] == [["point_data:t", "point_data:temp"]]
    with pytest.raises(ValueError, match="at least one of keep/drop/rename"):
        _tools.tool_data_manage(mesh_file, str(tmp_path / "dm2.vtu"))


def test_data_convert(mesh_file, tmp_path):
    out = _dump(
        _tools.tool_data_convert(mesh_file, str(tmp_path / "dc.vtu"), "point_to_cell")
    )
    assert "t" in out["cell_data"]
    with pytest.raises(ValueError, match="unknown direction"):
        _tools.tool_data_convert(mesh_file, str(tmp_path / "dc2.vtu"), "sideways")


def test_gradient(mesh_file, tmp_path):
    out = _dump(_tools.tool_gradient(mesh_file, str(tmp_path / "g.vtu"), "t"))
    assert "t:gradient" in out["cell_data"]
    # The mixed fixture's triangle block is below the tet block's dimension, so
    # the skip counter must be reported rather than silently swallowed.
    assert out["num_skipped"] == 2
    assert out["num_fallback"] == 0

    div = _dump(
        _tools.tool_gradient(
            mesh_file,
            str(tmp_path / "d.vtu"),
            "t",
            operator="gradient",
            location="point",
        )
    )
    assert "t:gradient" in div["point_data"]

    with pytest.raises(ValueError):
        _tools.tool_gradient(mesh_file, str(tmp_path / "x.vtu"), "nope")


def test_data_condition(mesh_file, tmp_path):
    out_path = str(tmp_path / "cond.vtu")
    _dump(
        _tools.tool_data_condition(
            output_path=out_path,
            input_path=mesh_file,
            operation="normalize",
            arrays=["t"],
            lo=0.0,
            hi=10.0,
        )
    )
    reread = meshioplusplus.read(out_path)
    assert np.isclose(reread.point_data["t"].max(), 10.0)


# --------------------------------------------------------------------------- #
# Parity guard: every public operation must be claimed by a tool              #
# --------------------------------------------------------------------------- #

# Public API names that deliberately have no path-based MCP tool. Add here
# only when the API is genuinely not expressible as file-in/file-out (the
# in-memory interop/GPU handoffs, the interactive viewer, classes, plumbing).
_NOT_TOOLS = {
    "_cli",
    "write_points_cells",  # constructor variant of write(); convert covers I/O
    "register_format",
    "deregister_format",
    "partition_labels",  # in-memory variant of partition
    # A generator of (time, Mesh) pairs: the meshes never touch the filesystem,
    # so there is no path-in/path-out form. The `sequence` tool covers every
    # file-level use (fan-in, fan-out, N->N).
    "read_sequence",
    "view",  # interactive; screenshot is the headless tool
    "has_viewer",
    "to_pyvista",
    "from_pyvista",
    "to_trimesh",
    "from_trimesh",
    # An in-memory handoff to a live Blender datablock: there is no
    # path-in/path-out form, and `convert` already covers the file case.
    "to_blender",
    "from_blender",
    "has_blender",
    "to_arrow",
    "from_arrow",
    "read_parquet",  # tabular import, no Mesh output to write
    "to_pandas",
    "to_polars",
    "has_pyvista",
    "has_trimesh",
    "has_arrow",
    "has_pandas",
    "has_polars",
    "has_open3d",
    "has_dolfinx",
    "edge_index",  # in-memory (2, E) array; no path-based form
    "feature_matrix",  # in-memory matrix + schema; no path-based form
    "FeatureMatrix",
    "has_zarr",
    "to_dlpack",
    "to_cupy",
    "from_cupy",
    "has_cupy",
    "has_cuda_device",
    "to_torch",
    "to_jax",
    "has_torch",
    "has_jax",
    "Mesh",
    "CellBlock",
    "Region",
    "TimeSeries",  # in-memory random-access class; the `sequence` tool covers
    # the path-in/path-out cases (fan-in/fan-out/per-step chain)
    "DatasetEntry",  # class; the dataset_* tools cover the path-based cases
    "ReadError",
    "WriteError",
    "topological_dimension",
    "__version__",
}


def test_every_operation_has_a_tool():
    everything = meshioplusplus.__all__
    formats = set(everything[: everything.index("_cli")])
    ops = set(everything) - formats - _NOT_TOOLS
    covered = set()
    for spec in TOOL_REGISTRY.values():
        covered.update(spec["wraps"])
    missing = sorted(ops - covered)
    assert not missing, (
        f"public operations without an MCP tool: {missing} — add a tool to "
        "meshioplusplus/mcp/_tools.py (and register it in _server.py) or, if "
        "the API is genuinely not path-expressible, add a conscious exemption "
        "to _NOT_TOOLS above. See CLAUDE.md 'Keep the MCP server in sync'."
    )


def test_wraps_names_are_real_public_api():
    public = set(meshioplusplus.__all__)
    for name, spec in TOOL_REGISTRY.items():
        stale = [w for w in spec["wraps"] if w not in public]
        assert not stale, f"tool '{name}' claims unknown API names: {stale}"


def test_registry_entries_are_wellformed():
    for name, spec in TOOL_REGISTRY.items():
        assert callable(spec["fn"]), name
        assert spec["gated"] in (None, "arrow", "viewer"), name


def test_every_tool_function_is_callable():
    """Import-level smoke: call every non-gated tool's `fn` attribute lookup.

    `_tools.py` imports the public API by name, so a tool added without its
    import is an `F821` that only fires when that tool actually *runs* — which a
    per-tool test suite catches only for the tools it happens to exercise. This
    walks the whole registry instead.
    """
    import inspect

    for name, spec in TOOL_REGISTRY.items():
        fn = spec["fn"]
        # Every name the function body resolves from module globals must exist.
        for var in fn.__code__.co_names:
            if var in fn.__globals__ or hasattr(_tools, var):
                continue
            # Attribute names and locals also land in co_names; only flag a bare
            # global that is neither defined nor a builtin.
            assert not (
                var.islower() and var in getattr(_tools, "__annotations__", {})
            ), f"{name}: unresolved global {var}"
        assert inspect.signature(fn) is not None


def test_compute_sdf_tool(tmp_path):
    # A *surface*: distance is measured to a skin, not to a volume.
    skin = tmp_path / "skin.vtu"
    meshioplusplus.write(
        skin,
        meshioplusplus.convert_cells(
            meshioplusplus.extract_surface(meshioplusplus.grid([2, 2, 2])),
            mode="simplexify",
        ),
    )
    mesh_file = str(skin)
    out = tmp_path / "field.vtu"
    report = _tools.tool_compute_sdf(mesh_file, str(out), resolution=[4, 4, 4])
    assert out.exists()
    assert report["dims"] == [4, 4, 4]
    assert report["max_depth"] == 0
    assert report["num_cells"] == 64

    tree = tmp_path / "tree.vtu"
    report = _tools.tool_compute_sdf(
        mesh_file, str(tree), structure="octree", root_resolution=4, max_depth=2
    )
    assert report["max_depth"] == 2
    assert 64 < report["num_cells"] < 4096

    # resolution/cell_size size a voxel grid; an octree's finest cell is already
    # determined by root_resolution and max_depth.
    with pytest.raises(ValueError, match="octree"):
        _tools.tool_compute_sdf(
            mesh_file, str(tree), structure="octree", resolution=[4, 4, 4]
        )


def test_crop_tool_takes_a_data_predicate(mesh_file, tmp_path):
    tagged = tmp_path / "tagged.vtu"
    mesh = meshioplusplus.read(mesh_file)
    mesh.cell_data["t"] = [
        np.arange(len(cb.data), dtype=np.float64) for cb in mesh.cells
    ]
    meshioplusplus.write(tagged, mesh)

    out = tmp_path / "kept.vtu"
    report = _tools.tool_crop(
        str(tagged), str(out), where_array="t", where_compare="<", where_value=1.0
    )
    assert report["num_cells"] == len(mesh.cells)  # one cell per block matches


# --------------------------------------------------------------------------- #
# Gated half: the FastMCP server (needs the [mcp] extra)                      #
# --------------------------------------------------------------------------- #


def _server():
    pytest.importorskip("mcp")
    import meshioplusplus.mcp as mmcp

    return mmcp.create_server()


def _run(coro):
    import asyncio

    return asyncio.run(coro)


def _tool_json(result):
    # FastMCP.call_tool returns a content list (older SDKs) or a
    # (content, structured) tuple (newer); normalize to the parsed JSON.
    content = result[0] if isinstance(result, tuple) else result
    text = next(c.text for c in content if getattr(c, "type", "") == "text")
    return json.loads(text)


def test_server_lists_every_registered_tool():
    server = _server()
    tools = _run(server.list_tools())
    assert sorted(t.name for t in tools) == sorted(TOOL_REGISTRY)
    for t in tools:
        assert (t.description or "").strip(), f"tool '{t.name}' has no description"
        assert t.inputSchema, f"tool '{t.name}' has no input schema"


def test_server_call_tool_info_and_convert(mesh_file, tmp_path):
    server = _server()
    report = _tool_json(_run(server.call_tool("info", {"input_path": mesh_file})))
    assert report["num_points"] == 5
    out_path = str(tmp_path / "out.vtk")
    report = _tool_json(
        _run(
            server.call_tool(
                "convert", {"input_path": mesh_file, "output_path": out_path}
            )
        )
    )
    assert os.path.isfile(report["output_path"])


def test_server_error_payload_shape(tmp_path):
    server = _server()
    report = _tool_json(
        _run(server.call_tool("info", {"input_path": str(tmp_path / "nope.vtu")}))
    )
    assert report["error_type"] == "ValueError"
    assert "input file not found" in report["error"]


def test_server_gated_tool_names_the_extra(mesh_file, tmp_path):
    if meshioplusplus.has_viewer():
        pytest.skip("polyscope installed; the gated error path is not reachable")
    server = _server()
    report = _tool_json(
        _run(
            server.call_tool(
                "screenshot",
                {"input_path": mesh_file, "output_path": str(tmp_path / "s.png")},
            )
        )
    )
    assert report["error_type"] == "ImportError"
    assert "polyscope" in report["error"]


def test_server_resources():
    server = _server()
    resources = _run(server.list_resources())
    uris = sorted(str(r.uri) for r in resources)
    assert uris == ["meshioplusplus://formats", "meshioplusplus://version"]
    content = _run(server.read_resource("meshioplusplus://formats"))
    payload = json.loads(list(content)[0].content)
    assert payload == _tools.formats_payload()


def test_has_mcp_and_named_error_without_it():
    import meshioplusplus.mcp as mmcp

    assert isinstance(mmcp.has_mcp(), bool)
    if not mmcp.has_mcp():
        with pytest.raises(ImportError, match=r"pip install meshioplusplus\[mcp\]"):
            mmcp.create_server()


# --------------------------------------------------------------------------- #
# Pure half: dataset-manifest tools                                           #
# --------------------------------------------------------------------------- #
def _case_files(tmp_path, n=3):
    cases = tmp_path / "cases"
    cases.mkdir(exist_ok=True)
    for i in range(n):
        meshioplusplus.write(str(cases / f"case_{i}.vtu"), _mixed_mesh())
    return cases


def test_dataset_add_list_update_round_trip(tmp_path, monkeypatch):
    _case_files(tmp_path)
    monkeypatch.chdir(tmp_path)
    manifest = str(tmp_path / "m.json")

    out = _dump(
        _tools.tool_dataset_add(
            manifest,
            input_pattern="cases/case_*.vtu",
            entry_id="sweep",
            tags=["raw"],
            metadata={"Re": 100},
        )
    )
    assert out["entry_id"] == "sweep" and out["num_steps"] == 3
    # the stored source is manifest-relative AND "/"-separated regardless of
    # platform, so the manifest is genuinely portable (doc/datasets.md) --
    # not os.path.join, which would be backslash-joined on Windows.
    doc = json.loads(open(manifest, encoding="utf-8").read())
    assert doc["Entries"][0]["Source"]["Pattern"] == "cases/case_*.vtu"

    _dump(
        _tools.tool_dataset_add(
            manifest,
            input_paths=["cases/case_0.vtu", "cases/case_1.vtu"],
            entry_id="pair",
            group="g/h",
        )
    )
    out = _dump(
        _tools.tool_dataset_update(
            manifest, entry_ids=["sweep"], split="train", add_tags=["v1"]
        )
    )
    assert out["splits"] == {"train": 1, "None": 1}
    out = _dump(_tools.tool_dataset_list(manifest, split="train", resolve=True))
    assert out["num_matching"] == 1
    entry = out["entries"][0]
    assert entry["Id"] == "sweep" and entry["Tags"] == ["raw", "v1"]
    assert len(entry["resolved"]) == 3
    assert all(item["time_source"] for item in entry["resolved"])


def test_dataset_update_assign_and_annotate(tmp_path, monkeypatch):
    _case_files(tmp_path, n=6)
    monkeypatch.chdir(tmp_path)
    manifest = str(tmp_path / "m.json")
    for i in range(6):
        _tools.tool_dataset_add(
            manifest, input_paths=[f"cases/case_{i}.vtu"], entry_id=f"c{i}"
        )
    out = _dump(
        _tools.tool_dataset_update(
            manifest, assign_splits={"train": 0.5, "test": 0.5}, seed=1
        )
    )
    assert out["splits"] == {"train": 3, "test": 3}
    _dump(
        _tools.tool_dataset_update(
            manifest, entry_ids=["c0"], notes="odd", metadata={"Ma": 0.3}
        )
    )
    m = meshioplusplus.DatasetManifest.load(manifest)
    assert m["c0"].notes == "odd" and m["c0"].metadata == {"Ma": 0.3}
    with pytest.raises(ValueError, match="exactly one entry_id"):
        _tools.tool_dataset_update(manifest, all_entries=True, notes="x")
    with pytest.raises(ValueError, match="entry_ids or all_entries"):
        _tools.tool_dataset_update(manifest, split="train")


def test_dataset_tools_validate_input_shape(tmp_path):
    manifest = str(tmp_path / "m.json")
    with pytest.raises(ValueError, match="exactly one of input_pattern"):
        _tools.tool_dataset_add(manifest)
    with pytest.raises(ValueError, match="exactly one of input_pattern"):
        _tools.tool_dataset_add(manifest, input_pattern="a*.vtu", input_paths=["a.vtu"])


def test_dataset_add_sandboxes_the_sources(tmp_path, monkeypatch):
    # A source outside the root must fail exactly the way a path argument
    # would — before anything is written to the manifest.
    root = tmp_path / "root"
    root.mkdir()
    outside = tmp_path / "outside.vtu"
    meshioplusplus.write(str(outside), _mixed_mesh())
    monkeypatch.setattr(_tools, "_ROOT", str(root))
    manifest = str(root / "m.json")
    with pytest.raises(ValueError, match="outside the configured root"):
        _tools.tool_dataset_add(manifest, input_paths=[str(outside)])
    assert not os.path.exists(manifest)


def test_dataset_list_sandboxes_hand_edited_sources(tmp_path, monkeypatch):
    # A hand-edited manifest naming files outside the root: listing is fine
    # (the document is data), but resolving its plan must refuse — the
    # pipeline-tool rule applied to manifests.
    root = tmp_path / "root"
    root.mkdir()
    outside = tmp_path / "outside.vtu"
    meshioplusplus.write(str(outside), _mixed_mesh())
    manifest = root / "m.json"
    manifest.write_text(
        json.dumps(
            {
                "Version": 1,
                "Entries": [{"Id": "x", "Source": {"Path": str(outside)}}],
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(_tools, "_ROOT", str(root))
    _dump(_tools.tool_dataset_list(str(manifest)))  # data-only: allowed
    with pytest.raises(ValueError, match="outside the configured root"):
        _tools.tool_dataset_list(str(manifest), resolve=True)
