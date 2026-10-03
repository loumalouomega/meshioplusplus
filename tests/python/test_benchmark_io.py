"""The I/O sweep must exercise its explicitly supported single-type formats."""

import copy
import csv
import importlib.util
import itertools
import math
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]


@pytest.fixture
def bench(monkeypatch):
    # bench.py also imports the sibling inputs and the conformance declaration.
    # Restore its sys.path changes when the test finishes.
    monkeypatch.setattr(sys, "path", sys.path.copy())
    spec = importlib.util.spec_from_file_location(
        "benchmark_io", ROOT / "benchmark/bench.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_mfm_single_type_selection_does_not_change_conformance(bench):
    import conformance_spec as cs

    before = copy.deepcopy(cs.SPEC["mfm"])
    formats = bench.all_formats()
    assert formats.count(("mfm", "volume")) == 1
    assert cs.SPEC["mfm"] == before
    assert before["error"] == "WriteError"
    assert set(before["cells"].values()) == {"?"}


def test_unknown_conformance_cells_are_not_implicitly_benchmarked(bench, monkeypatch):
    import conformance_spec as cs

    monkeypatch.setattr(
        cs,
        "SPEC",
        {
            "mfm": cs.SPEC["mfm"],
            "untested": {"cells": {"tetra": "?"}, "error": "WriteError"},
        },
    )
    assert bench.all_formats() == [("mfm", "volume")]


@pytest.mark.parametrize("formats", [None, ["mfm"]])
def test_mfm_sweep_produces_real_rows(bench, monkeypatch, formats):
    import conformance_spec as cs

    monkeypatch.setattr(cs, "SPEC", {"mfm": cs.SPEC["mfm"]})
    monkeypatch.setitem(bench.SIZES, "S", 3)
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    rows = bench.run_all(sizes=["S"], formats=formats, repeats=1)
    assert len(rows) == 1
    row = rows[0]
    assert (row["format"], row["size"], row["input"]) == ("mfm", "S", "volume")
    assert row["cells"] == 6 * (3 - 1) ** 3
    assert row["bytes"] > 0
    assert row["read_s"] > 0
    assert row["write_s"] > 0


def test_mfm_cli_writes_a_nonempty_csv(bench, monkeypatch, tmp_path):
    monkeypatch.setitem(bench.SIZES, "S", 3)
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    target = tmp_path / "mfm.csv"
    assert bench.main(["--formats", "mfm", "--repeats", "1", "--out", str(target)]) == 0
    with target.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 1
    assert rows[0]["format"] == "mfm"
    assert int(rows[0]["cells"]) == 48


def test_mfm_benchmark_input_roundtrips_natively(bench, monkeypatch, tmp_path):
    from meshioplusplus import _core
    from meshioplusplus.mfm import _py_read

    monkeypatch.setitem(bench.SIZES, "S", 3)
    mesh = bench._inputs("S")["volume"]
    target = tmp_path / "grid.mfm"
    _core.mfm_write(str(target), mesh, ".16e")
    native = _core.mfm_read(str(target))
    reference = _py_read(target)
    for result in (native, reference):
        np.testing.assert_array_equal(result.points, mesh.points)
        assert [block.type for block in result.cells] == ["tetra"]
        np.testing.assert_array_equal(result.cells[0].data, mesh.cells[0].data)
        np.testing.assert_array_equal(
            result.cell_data["mfm:ref"][0], np.ones(48, dtype=np.int64)
        )


@pytest.mark.parametrize("binary", [False, True])
def test_gmsh_single_type_benchmark_keeps_mixed_conformance(
    bench, monkeypatch, tmp_path, binary
):
    import conformance_spec as cs

    from meshioplusplus import _core, gmsh

    original = copy.deepcopy(cs.SPEC["gmsh"])
    assert bench.all_formats().count(("gmsh", "volume")) == 1
    monkeypatch.setitem(bench.SIZES, "S", 3)
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    source = bench._inputs("S")["volume"]
    path = tmp_path / "mesh.msh"
    gmsh.write(path, source, binary=binary)
    result = _core.gmsh_read(str(path))
    np.testing.assert_array_equal(result.points, source.points)
    assert [block.type for block in result.cells] == ["tetra"]
    np.testing.assert_array_equal(result.cells[0].data, source.cells[0].data)
    assert cs.SPEC["gmsh"] == original
    assert original["error"] == "WriteError"


def test_dex_sweep_supplies_required_field_without_mutating_inputs(bench, monkeypatch):
    import conformance_spec as cs

    original = copy.deepcopy(cs.SPEC["dex"])
    monkeypatch.setattr(cs, "SPEC", {"dex": cs.SPEC["dex"]})
    monkeypatch.setitem(bench.SIZES, "S", 3)
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    assert bench.all_formats() == [("dex", "points")]
    source = bench._inputs("S")["points"]
    prepared = bench._input_for_format(source, "dex")
    assert source.point_data == {}
    np.testing.assert_array_equal(prepared.point_data["benchmark:field"], np.arange(27))
    rows = bench.run_all(sizes=["S"], formats=["dex"], repeats=1)
    assert len(rows) == 1 and rows[0]["bytes"] > 0
    assert rows[0]["input"] == "points"
    assert cs.SPEC["dex"] == original


def test_mfm_evidence_has_a_complete_matching_digest_matrix():
    """Check evidence completeness and parity, never a noisy timing threshold."""
    with (ROOT / "benchmark/text_io_mfm.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    configs = [("seq", "1"), *itertools.product(("openmp", "tbb"), ("1", "4", "8"))]
    expected = {
        (stage, backend, threads, size)
        for stage, (backend, threads), size in itertools.product(
            ("before", "after"), configs, ("M", "L")
        )
    }
    actual, digests, modules = {}, {}, {}
    for row in rows:
        key = (row["stage"], row["backend"], row["threads"], row["size"])
        assert key not in actual
        actual[key] = row
        assert (
            row["source_commit"]
            == {
                "before": "3b0791205bd68cdd50d3127d1aa9d08047ba9900",
                "after": "2d39e4020d2cac0a8372e170eb9ce92296741614",
            }[row["stage"]]
        )
        assert row["mesh_backend"] == "MESHIO"
        assert row["timing_repeats"] == "7"
        assert row["timing_warmup"] == "1"
        assert row["allocation_repeats"] == "3"
        for name in ("file_sha256", "mesh_sha256", "module_sha256"):
            digest = row[name]
            assert len(digest) == 64 and all(c in "0123456789abcdef" for c in digest)
        pair = (row["file_sha256"], row["mesh_sha256"])
        assert digests.setdefault(row["size"], pair) == pair
        assert (
            modules.setdefault((row["stage"], row["backend"]), row["module_sha256"])
            == row["module_sha256"]
        )
        assert int(row["cells"]) == {"M": 257250, "L": 998250}[row["size"]]
        assert int(row["allocations"]) > 0
        assert int(row["requested_bytes"]) > 0
        for name in ("read_s", "write_s"):
            assert math.isfinite(float(row[name])) and float(row[name]) > 0
    assert actual.keys() == expected
    assert len(set(digests.values())) == 2
    for (backend, threads), size in itertools.product(configs, ("M", "L")):
        before = actual[("before", backend, threads, size)]
        after = actual[("after", backend, threads, size)]
        assert int(after["allocations"]) < int(before["allocations"])
        assert int(after["requested_bytes"]) < int(before["requested_bytes"])


def test_mfm_seq_confirmation_keeps_both_measurement_orders():
    with (ROOT / "benchmark/text_io_mfm_seq_confirmation.csv").open(
        newline=""
    ) as stream:
        rows = list(csv.DictReader(stream))
    expected = set(itertools.product(("1", "2"), ("before", "after"), ("M", "L")))
    actual = {(row["round"], row["stage"], row["size"]) for row in rows}
    assert len(rows) == len(actual) == 8
    assert actual == expected
    for row in rows:
        assert row["backend"] == "seq" and row["threads"] == "1"
        assert row["timing_repeats"] == "7"
        for name in ("read_s", "write_s"):
            assert math.isfinite(float(row[name])) and float(row[name]) > 0


@pytest.mark.parametrize(
    "batch,formats,decreases",
    [
        ("batch2", ("off", "medit", "xyz", "unv"), ("off", "medit", "xyz", "unv")),
        (
            "batch3",
            ("ensight", "vtk", "tecplot", "ansys", "patran"),
            ("ensight", "vtk", "ansys"),
        ),
        ("batch4", ("femap", "pcd", "radioss"), ("femap", "pcd")),
        ("batch5", ("dex", "ip", "flux", "permas"), ("dex", "ip", "flux", "permas")),
        (
            "batch6",
            ("abaqus", "ansysInp", "marc", "nastran", "netgen"),
            ("abaqus", "ansysInp", "marc", "nastran", "netgen"),
        ),
        ("batch7", ("lsdyna", "radioss", "gid"), ("lsdyna",)),
        (
            "batch8",
            ("lsdyna", "ansysInp", "marc", "radioss", "patran", "nastran", "unv"),
            ("lsdyna", "ansysInp", "marc", "radioss", "nastran"),
        ),
        (
            "batch9",
            ("gmsh", "vtk", "ensight", "unv", "mfem", "mdpa", "gid", "femap"),
            (),
        ),
        (
            "batch10",
            ("abaqus", "nastran", "netgen", "flac3d"),
            ("abaqus", "nastran", "netgen", "flac3d"),
        ),
        (
            "batch11",
            ("gmsh", "obj", "ply", "openfoam", "dex"),
            ("obj", "dex"),
        ),
        ("batch12_float64_int64", ("off", "ip", "flux", "permas"), ()),
    ],
)
def test_text_io_evidence_is_complete_and_bit_identical(batch, formats, decreases):
    with (ROOT / f"benchmark/text_io_{batch}.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    configs = [("seq", "1")] + list(
        itertools.product(("openmp", "tbb"), ("1", "4", "8"))
    )
    expected = {
        (stage, backend, threads, fmt, size)
        for stage, (backend, threads), fmt, size in itertools.product(
            ("before", "after"), configs, formats, ("M", "L")
        )
    }
    actual, digests, modules = {}, {}, {}
    for row in rows:
        key = tuple(row[k] for k in ("stage", "backend", "threads", "format", "size"))
        assert key not in actual
        actual[key] = row
        digest = row["file_sha256"], row["mesh_sha256"]
        assert digests.setdefault((row["format"], row["size"]), digest) == digest
        assert (
            modules.setdefault((row["stage"], row["backend"]), row["module_sha256"])
            == row["module_sha256"]
        )
        assert row["timing_repeats"] == "7" and row["allocation_repeats"] == "3"
        if batch == "batch9":
            assert row["reader"] == (
                "mdpa_read" if row["format"] == "mdpa" else "strict Python dispatch"
            )
        for name in ("read_s", "write_s"):
            assert math.isfinite(float(row[name])) and float(row[name]) > 0
        for name in ("file_sha256", "mesh_sha256", "module_sha256"):
            assert len(row[name]) == 64 and int(row[name], 16) >= 0
    assert actual.keys() == expected
    for (backend, threads), fmt, size in itertools.product(
        configs, formats, ("M", "L")
    ):
        before = actual[("before", backend, threads, fmt, size)]
        after = actual[("after", backend, threads, fmt, size)]
        for field in ("allocations", "requested_bytes"):
            byte_only = {"batch7": ("radioss", "gid"), "batch8": ("patran",)}
            if fmt in decreases or (
                field == "requested_bytes" and fmt in byte_only.get(batch, ())
            ):
                assert int(after[field]) < int(before[field])
            else:
                assert int(after[field]) == int(before[field])
    with (ROOT / f"benchmark/text_io_{batch}_seq_confirmation.csv").open(
        newline=""
    ) as stream:
        confirm = list(csv.DictReader(stream))
    keys = {(row["round"], row["stage"], row["format"], row["size"]) for row in confirm}
    assert len(confirm) == len(keys) == 8 * len(formats)
    assert keys == set(
        itertools.product(("1", "2"), ("before", "after"), formats, ("M", "L"))
    )
    for row in confirm:
        assert row["backend"] == "seq" and row["threads"] == "1"
        assert (row["file_sha256"], row["mesh_sha256"]) == digests[
            row["format"], row["size"]
        ]
