"""Capture, validate and publish immutable benchmark runs (stdlib only)."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.metadata
import io
import json
import math
import os
import pathlib
import platform
import re
import shlex
import shutil
import subprocess

FIELDS = {
    "formats": {"size", "format", "input", "cells", "bytes", "write_s", "read_s"},
    "operations": {"backend", "threads", "op", "cells", "median_s", "runs", "digest"},
}


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"


def command(*args):
    try:
        return subprocess.check_output(
            args, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def build_metadata(directory):
    """Record the configuration and actual flags of a measured native build."""
    directory = pathlib.Path(directory)
    cache = {}
    for line in (directory / "CMakeCache.txt").read_text().splitlines():
        if line and not line.startswith(("#", "//")) and "=" in line:
            key, value = line.split("=", 1)
            cache[key.split(":", 1)[0]] = value
    entries = json.loads((directory / "compile_commands.json").read_text())
    flags = next(
        entry.get("arguments", entry.get("command"))
        for entry in entries
        if entry["file"].endswith("/benchmark/bench_ops.cpp")
    )
    # Paths vary per runner/worktree; keep flags but not source/output paths.
    tokens = flags if isinstance(flags, list) else shlex.split(flags)
    return {
        "compiler": command(cache["CMAKE_CXX_COMPILER"], "--version"),
        "build_type": cache["CMAKE_BUILD_TYPE"],
        "options": {
            key: value
            for key, value in sorted(cache.items())
            if key.startswith("MESHIOPLUSPLUS_")
            and not key.endswith(("_DIR", "_ROOT", "_LIBRARY", "_INCLUDE_DIR"))
        },
        "flags": [
            token for token in tokens if token.startswith(("-D", "-O", "-f", "-std="))
        ],
    }


def environment(kind, builds=None):
    cpu = platform.processor()
    cpuinfo = pathlib.Path("/proc/cpuinfo")
    if cpuinfo.exists():
        cpu = next(
            (
                s.split(":", 1)[1].strip()
                for s in cpuinfo.read_text().splitlines()
                if s.startswith("model name")
            ),
            cpu,
        )
    libraries = {}
    for name in ("libhdf5-dev", "libnetcdf-dev", "zlib1g-dev", "libtbb-dev"):
        libraries[name] = command("dpkg-query", "-W", "-f=${Version}", name)
    result = {
        "os": platform.system(),
        "arch": platform.machine(),
        "cpu": cpu,
        "cpus": os.cpu_count(),
        "runner": os.environ.get("RUNNER_ENVIRONMENT", "local"),
        "image": os.environ.get("ImageOS", "unknown"),
        "compiler": command(os.environ.get("CXX", "c++"), "--version").splitlines()[0],
        "libraries": libraries,
    }
    if kind == "formats":
        from meshioplusplus import _core

        result["python"] = platform.python_version()
        result["python_libraries"] = {}
        for name in ("numpy", "h5py", "netCDF4"):
            try:
                result["python_libraries"][name] = importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError:
                result["python_libraries"][name] = "absent"
        result["backend"] = _core.__parallel_backend__
        result["mesh_backend"] = _core.__mesh_backend__
        result["capabilities"] = {
            name: bool(getattr(_core, name))
            for name in dir(_core)
            if name.startswith("__has_")
        }
        result["threads"] = os.environ.get("OMP_NUM_THREADS", "unset")
    else:
        result["builds"] = {
            backend: build_metadata(directory)
            for backend, directory in sorted((builds or {}).items())
        }
    return result


def validate(record):
    if record.get("schema") != 1 or record.get("kind") not in FIELDS:
        raise ValueError("unsupported benchmark schema/kind")
    for key in ("run_id", "attempt"):
        if not re.fullmatch(r"[1-9][0-9]*", str(record.get(key, ""))):
            raise ValueError(f"invalid {key}")
    if not re.fullmatch(r"[0-9a-f]{40}", record.get("sha", "")):
        raise ValueError("invalid commit SHA")
    if not re.fullmatch(r"[0-9a-f]{64}", record.get("raw_csv_sha256", "")):
        raise ValueError("missing/invalid raw CSV checksum")
    stamp = dt.datetime.fromisoformat(record["timestamp"].replace("Z", "+00:00"))
    if stamp.tzinfo is None:
        raise ValueError("timestamp must have a timezone")
    if not isinstance(record.get("environment"), dict) or not record["environment"]:
        raise ValueError("missing environment")
    if not isinstance(record.get("parameters"), dict) or not record["parameters"]:
        raise ValueError("missing parameters")
    if not isinstance(record.get("rows"), list) or not record["rows"]:
        raise ValueError("empty benchmark run")
    kind = record["kind"]
    seen = set()
    for row in record["rows"]:
        if set(row) != FIELDS[kind] or not all(
            isinstance(v, str) for v in row.values()
        ):
            raise ValueError("unexpected CSV columns/types")
        integers = (
            ("cells", "bytes") if kind == "formats" else ("cells", "threads", "runs")
        )
        for key in integers:
            if not row[key].isdigit() or (
                key in ("threads", "runs") and int(row[key]) == 0
            ):
                raise ValueError(f"invalid {key}")
        times = ("write_s", "read_s") if kind == "formats" else ("median_s",)
        for key in times:
            value = float(row[key])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid {key}")
        if kind == "formats":
            if row["size"] not in ("S", "M", "L"):
                raise ValueError("invalid size")
            identity = tuple(row[k] for k in ("size", "format", "input", "cells"))
        else:
            if row["backend"] not in ("seq", "openmp", "tbb", "stl", "kokkos"):
                raise ValueError("invalid backend")
            if not re.fullmatch(r"[0-9a-fA-F]{16}", row["digest"]):
                raise ValueError("missing/invalid determinism digest")
            identity = tuple(row[k] for k in ("backend", "threads", "op", "cells"))
        if identity in seen:
            raise ValueError("duplicate measurement")
        seen.add(identity)
    canonical(record)  # reject non-finite metadata too
    return record


def csv_rows(raw, kind):
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8"), newline=""))
    if (
        reader.fieldnames is None
        or set(reader.fieldnames) != FIELDS[kind]
        or len(reader.fieldnames) != len(FIELDS[kind])
    ):
        raise ValueError("unexpected CSV header")
    return list(reader)


def raw_csv(path, record):
    try:
        raw = pathlib.Path(path).read_bytes()
    except OSError:
        raise ValueError("missing raw CSV sidecar")
    if hashlib.sha256(raw).hexdigest() != record.get("raw_csv_sha256", ""):
        raise ValueError("raw CSV checksum mismatch")
    if csv_rows(raw, record["kind"]) != record["rows"]:
        raise ValueError("raw CSV rows disagree with metadata")
    return raw


def capture(source, kind, output, parameters, build_root=None):
    source, output = pathlib.Path(source), pathlib.Path(output)
    raw = source.read_bytes()
    rows = csv_rows(raw, kind)
    builds = None
    if kind == "operations":
        if build_root is None:
            raise ValueError("operations capture needs --build-root")
        builds = {
            row["backend"]: pathlib.Path(build_root) / f"bench-ops-{row['backend']}"
            for row in rows
        }
    record = validate(
        {
            "schema": 1,
            "kind": kind,
            "run_id": os.environ["GITHUB_RUN_ID"],
            "attempt": os.environ["GITHUB_RUN_ATTEMPT"],
            "sha": os.environ["GITHUB_SHA"],
            "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
            "environment": environment(kind, builds),
            "parameters": parameters,
            "rows": rows,
            "raw_csv_sha256": hashlib.sha256(raw).hexdigest(),
        }
    )
    output.mkdir(parents=True, exist_ok=True)
    (output / f"{kind}.json").write_text(canonical(record))
    (output / f"{kind}.csv").write_bytes(raw)


def ingest(source, store):
    """Validate all inputs before writing; repeated identical records are no-ops."""
    source, store = pathlib.Path(source), pathlib.Path(store)
    pending = {}
    for path in sorted(source.rglob("*.json")):
        record = validate(json.loads(path.read_text()))
        directory = f"{record['run_id']}-{record['attempt']}"
        target = store / "runs" / directory / f"{record['kind']}.json"
        values = {
            target: canonical(record).encode(),
            target.with_suffix(".csv"): raw_csv(path.with_suffix(".csv"), record),
        }
        for destination, value in values.items():
            if destination in pending and pending[destination] != value:
                raise ValueError("conflicting input records")
            if destination.exists() and destination.read_bytes() != value:
                raise ValueError(f"immutable run already exists: {destination}")
            pending[destination] = value
    if not pending:
        raise ValueError("no benchmark records")
    for path, value in pending.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(value)


def export(store, output):
    """Produce static site data. Missing measurements stay absent, never zero."""
    store, output = pathlib.Path(store), pathlib.Path(output)
    records = []
    for path in sorted(store.glob("runs/*/*.json")):
        record = validate(json.loads(path.read_text()))
        expected = f"{record['run_id']}-{record['attempt']}/{record['kind']}.json"
        if path.relative_to(store / "runs").as_posix() != expected:
            raise ValueError("record location disagrees with metadata")
        raw_csv(path.with_suffix(".csv"), record)
        records.append(record)
    records.sort(
        key=lambda r: (r["timestamp"], int(r["run_id"]), int(r["attempt"]), r["kind"])
    )
    output.mkdir(parents=True, exist_ok=True)
    for record in records:
        # Library versions/machine changes split series; the project's release
        # version and commit do not, since those are what the graph compares.
        comparison = {
            "environment": record["environment"],
            "parameters": {
                key: value
                for key, value in record["parameters"].items()
                if key not in ("sizes", "tiers")
            },
        }
        record["environment_id"] = hashlib.sha256(
            canonical(comparison).encode()
        ).hexdigest()[:16]
    (output / "index.json").write_text(canonical({"schema": 1, "runs": records}))
    if (store / "runs").exists():
        shutil.copytree(store / "runs", output / "runs", dirs_exist_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    cap = commands.add_parser("capture")
    cap.add_argument("kind", choices=FIELDS)
    cap.add_argument("source", type=pathlib.Path)
    cap.add_argument("output", type=pathlib.Path)
    cap.add_argument("--parameters", required=True, type=json.loads)
    cap.add_argument("--build-root", type=pathlib.Path)
    for name in ("ingest", "export"):
        sub = commands.add_parser(name)
        sub.add_argument("source", type=pathlib.Path)
        sub.add_argument("output", type=pathlib.Path)
    args = parser.parse_args()
    if args.command == "capture":
        capture(args.source, args.kind, args.output, args.parameters, args.build_root)
    else:
        {"ingest": ingest, "export": export}[args.command](args.source, args.output)


if __name__ == "__main__":
    main()
