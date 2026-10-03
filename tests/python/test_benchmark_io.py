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


def test_mfm_before_evidence_has_a_complete_matching_digest_matrix():
    """Check evidence completeness and parity, never a noisy timing threshold."""
    with (ROOT / "benchmark/text_io_mfm.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    configs = [("seq", "1"), *itertools.product(("openmp", "tbb"), ("1", "4", "8"))]
    expected = {
        (backend, threads, size)
        for (backend, threads), size in itertools.product(configs, ("M", "L"))
    }
    actual, digests, modules = {}, {}, {}
    for row in rows:
        key = (row["backend"], row["threads"], row["size"])
        assert key not in actual
        actual[key] = row
        assert row["stage"] == "before"
        assert row["source_commit"] == "3b0791205bd68cdd50d3127d1aa9d08047ba9900"
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
            modules.setdefault(row["backend"], row["module_sha256"])
            == row["module_sha256"]
        )
        assert int(row["cells"]) == {"M": 257250, "L": 998250}[row["size"]]
        assert int(row["allocations"]) > 0
        assert int(row["requested_bytes"]) > 0
        for name in ("read_s", "write_s"):
            assert math.isfinite(float(row[name])) and float(row[name]) > 0
    assert actual.keys() == expected
    assert len(set(digests.values())) == 2
