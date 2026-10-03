"""The trend store must not silently overwrite history or invent timings."""

import copy
import csv
import hashlib
import importlib.util
import io
import json
import pathlib

import pytest

REPO = pathlib.Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "bench_trends", REPO / "tools/bench_trends.py"
)
trends = importlib.util.module_from_spec(spec)
spec.loader.exec_module(trends)


def record(kind="operations", run="123"):
    row = {
        "backend": "seq",
        "threads": "1",
        "op": "clean",
        "cells": "100",
        "median_s": "0.1",
        "runs": "5",
        "digest": "0123456789abcdef",
    }
    if kind == "formats":
        row = {
            "size": "M",
            "format": "vtu",
            "input": "volume",
            "cells": "100",
            "bytes": "1024",
            "write_s": "0.2",
            "read_s": "0.3",
        }
    rows = [row]
    data = {
        "schema": 1,
        "kind": kind,
        "run_id": run,
        "attempt": "1",
        "sha": "a" * 40,
        "timestamp": "2026-10-02T12:00:00+00:00",
        "environment": {"cpu": "test"},
        "parameters": {"runs": 5},
        "rows": rows,
    }
    data["raw_csv_sha256"] = hashlib.sha256(raw_csv_bytes(data)).hexdigest()
    return data


def raw_csv_bytes(data):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=sorted(trends.FIELDS[data["kind"]]))
    writer.writeheader()
    writer.writerows(data["rows"])
    return stream.getvalue().encode()


def write_record(directory, data):
    data = copy.deepcopy(data)
    raw = raw_csv_bytes(data)
    # Recompute the checksum when the caller mutated rows without updating it,
    # unless the caller explicitly removed it to test rejection.
    if "raw_csv_sha256" not in data:
        data["raw_csv_sha256"] = hashlib.sha256(raw).hexdigest()
    kind = data["kind"]
    (directory / f"{kind}.json").write_text(json.dumps(data))
    (directory / f"{kind}.csv").write_bytes(raw)


@pytest.mark.parametrize("kind", ["formats", "operations"])
def test_ingest_and_export_are_idempotent(tmp_path, kind):
    incoming, store, site = [tmp_path / name for name in ("in", "store", "site")]
    incoming.mkdir()
    data = record(kind)
    write_record(incoming, data)
    trends.ingest(incoming, store)
    trends.ingest(incoming, store)
    trends.export(store, site)
    index = json.loads((site / "index.json").read_text())
    assert len(index["runs"]) == 1
    assert index["runs"][0]["rows"] == data["rows"]
    assert index["runs"][0]["raw_csv_sha256"] == data["raw_csv_sha256"]
    assert (site / "runs/123-1" / f"{kind}.csv").exists()


def test_conflicting_history_is_not_overwritten(tmp_path):
    incoming, store = tmp_path / "in", tmp_path / "store"
    incoming.mkdir()
    data = record()
    write_record(incoming, data)
    trends.ingest(incoming, store)
    original = (store / "runs/123-1/operations.json").read_bytes()
    data["rows"][0]["median_s"] = "10"
    data.pop("raw_csv_sha256")
    data.pop("_raw_csv", None)
    # Recompute the sidecar so JSON and CSV agree, but differ from history.
    write_record(incoming, data)
    with pytest.raises(ValueError, match="immutable"):
        trends.ingest(incoming, store)
    assert (store / "runs/123-1/operations.json").read_bytes() == original


@pytest.mark.parametrize(
    "field,value",
    [
        ("run_id", "../escape"),
        ("attempt", "0"),
        ("sha", "not-a-commit"),
        ("schema", 2),
        ("timestamp", "2026-10-02T12:00:00"),
        ("rows", []),
        ("raw_csv_sha256", "not-a-hash"),
    ],
)
def test_invalid_record_is_rejected(field, value):
    data = record()
    data[field] = value
    with pytest.raises(ValueError):
        trends.validate(data)


@pytest.mark.parametrize("value", ["NaN", "Infinity", "-1"])
def test_invalid_measurement_is_rejected(value):
    data = record()
    data["rows"][0]["median_s"] = value
    with pytest.raises(ValueError, match="median_s"):
        trends.validate(data)


def test_duplicate_measurement_is_rejected():
    data = record()
    data["rows"].append(copy.deepcopy(data["rows"][0]))
    with pytest.raises(ValueError, match="duplicate"):
        trends.validate(data)


def test_validation_happens_before_any_write(tmp_path):
    incoming, store = tmp_path / "in", tmp_path / "store"
    incoming.mkdir()
    write_named(incoming, "a", record())
    write_named(incoming, "b", record(run="../invalid"))
    with pytest.raises(ValueError):
        trends.ingest(incoming, store)
    assert not store.exists()


def test_environments_and_attempts_are_separated(tmp_path):
    incoming, store, site = [tmp_path / name for name in ("in", "store", "site")]
    incoming.mkdir()
    a, b, c = record(), record(run="124"), record()
    b["environment"]["cpu"] = "other CPU"
    c["attempt"] = "2"
    c["sha"] = "b" * 40
    for i, data in enumerate((a, b, c)):
        # Recompute the checksum after mutating environment/attempt metadata
        # (CSV content is unchanged, so the hash stays valid).
        data["raw_csv_sha256"] = hashlib.sha256(raw_csv_bytes(data)).hexdigest()
        write_named(incoming, str(i), data)
    trends.ingest(incoming, store)
    trends.export(store, site)
    records = json.loads((site / "index.json").read_text())["runs"]
    assert records[0]["environment_id"] == records[1]["environment_id"]
    assert records[0]["environment_id"] != records[2]["environment_id"]
    assert (store / "runs/123-2/operations.json").exists()


def write_named(directory, stem, data):
    (directory / f"{stem}.json").write_text(json.dumps(data))
    (directory / f"{stem}.csv").write_bytes(raw_csv_bytes(data))


def test_empty_export_supports_first_deployment(tmp_path):
    trends.export(tmp_path / "absent", tmp_path / "site")
    assert json.loads((tmp_path / "site/index.json").read_text()) == {
        "schema": 1,
        "runs": [],
    }


def test_missing_csv_sidecar_is_rejected(tmp_path):
    incoming, store = tmp_path / "in", tmp_path / "store"
    incoming.mkdir()
    (incoming / "operations.json").write_text(json.dumps(record()))
    with pytest.raises(ValueError, match="sidecar"):
        trends.ingest(incoming, store)
    assert not store.exists()


def test_checksum_mismatch_is_rejected(tmp_path):
    incoming, store = tmp_path / "in", tmp_path / "store"
    incoming.mkdir()
    data = record()
    data["raw_csv_sha256"] = "0" * 64
    write_named(incoming, "operations", data)
    with pytest.raises(ValueError, match="checksum"):
        trends.ingest(incoming, store)
    assert not store.exists()


def test_csv_rows_must_match_metadata(tmp_path):
    incoming, store = tmp_path / "in", tmp_path / "store"
    incoming.mkdir()
    data = record()
    # CSV with different timing than the JSON rows, but with a matching
    # checksum so the row comparison is what fails.
    other = copy.deepcopy(data["rows"])
    other[0]["median_s"] = "9.9"
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=sorted(trends.FIELDS[data["kind"]]))
    writer.writeheader()
    writer.writerows(other)
    raw = stream.getvalue().encode()
    data["raw_csv_sha256"] = hashlib.sha256(raw).hexdigest()
    (incoming / "operations.json").write_text(json.dumps(data))
    (incoming / "operations.csv").write_bytes(raw)
    with pytest.raises(ValueError, match="disagree"):
        trends.ingest(incoming, store)
    assert not store.exists()


def test_unexpected_csv_header_is_rejected(tmp_path):
    incoming, store = tmp_path / "in", tmp_path / "store"
    incoming.mkdir()
    data = record()
    (incoming / "operations.json").write_text(json.dumps(data))
    (incoming / "operations.csv").write_bytes(b"backend,threads\nseq,1\n")
    with pytest.raises(ValueError, match="header|checksum"):
        trends.ingest(incoming, store)
    assert not store.exists()


def test_record_location_must_match_metadata(tmp_path):
    incoming, store, site = [tmp_path / name for name in ("in", "store", "site")]
    incoming.mkdir()
    write_named(incoming, "operations", record())
    trends.ingest(incoming, store)
    # Move the record away from runs/<run>-<attempt>/<kind>.json.
    src = store / "runs/123-1/operations.json"
    dst = store / "runs/999-9/operations.json"
    dst.parent.mkdir(parents=True)
    dst.write_bytes(src.read_bytes())
    (store / "runs/999-9/operations.csv").write_bytes(
        (store / "runs/123-1/operations.csv").read_bytes()
    )
    src.unlink()
    (store / "runs/123-1/operations.csv").unlink()
    with pytest.raises(ValueError, match="location"):
        trends.export(store, site)
