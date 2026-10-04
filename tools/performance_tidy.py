"""Run the installed Clang-Tidy performance checks without modifying sources."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import pathlib
import re
import shlex
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
CHECKS = "-*,performance-*"
DIAGNOSTIC = re.compile(r"^(.*):(\d+):(\d+): (warning|error): (.*?) \[([^\]]+)\]$")


def tidy_config():
    return {
        "Checks": CHECKS,
        "WarningsAsErrors": "",
        "HeaderFilterRegex": re.escape(str(ROOT))
        + r"/(src/cpp/(include|src|cli)|bindings)/",
    }


def tidy_arguments():
    # Explicit configuration keeps unrelated .clang-tidy policy out of this audit.
    return [f"--config={json.dumps(tidy_config())}", f"--checks={CHECKS}"]


def first_party(path):
    """Exclude vendored code, generated headers and non-C/C++ translation units."""
    try:
        rel = path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return False
    return "third_party" not in path.parts and rel.startswith(
        ("src/cpp/src/", "src/cpp/include/", "src/cpp/cli/", "bindings/")
    )


def source_path(entry):
    path = pathlib.Path(entry["file"])
    return (pathlib.Path(entry["directory"]) / path).resolve()


def configure(build_dir, backend):
    """A separate Clang/SEQ tree avoids GCC-only module-scanning arguments."""
    subprocess.run(
        [
            "cmake",
            "-S",
            str(ROOT),
            "-B",
            str(build_dir),
            "-G",
            "Ninja",
            "-DCMAKE_BUILD_TYPE=Release",
            f"-DCMAKE_C_COMPILER={os.environ.get('CC', 'clang')}",
            f"-DCMAKE_CXX_COMPILER={os.environ.get('CXX', 'clang++')}",
            "-DCMAKE_CXX_SCAN_FOR_MODULES=OFF",
            "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
            "-DMESHIOPLUSPLUS_BUILD_PYTHON=OFF",
            "-DMESHIOPLUSPLUS_BUILD_C_API=ON",
            "-DMESHIOPLUSPLUS_BUILD_CLI=ON",
            "-DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ",
            f"-DMESHIOPLUSPLUS_MESH_BACKEND={backend}",
        ],
        check=True,
    )


def configuration(build_dir, entries):
    cache = build_dir / "CMakeCache.txt"
    values = {}
    if cache.exists():
        for line in cache.read_text().splitlines():
            match = re.match(r"([^:#]+):[^=]+=(.*)", line)
            if match and (
                match[1].startswith("MESHIOPLUSPLUS_")
                or match[1]
                in ("CMAKE_BUILD_TYPE", "CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER")
            ):
                values[match[1]] = match[2]
    defines = set()
    for entry in entries:
        args = entry.get("arguments")
        if args is None:
            args = shlex.split(entry["command"])
        defines.update(arg[2:] for arg in args if arg.startswith("-D"))
    # The defines, not just WITH_*=ON cache requests, show what was really found.
    return {"cmake_cache": values, "compile_defines": sorted(defines)}


def parse_diagnostics(output, directory):
    records = []
    for line in output.splitlines():
        match = DIAGNOSTIC.match(line)
        if not match:
            continue
        path = (directory / match[1]).resolve()
        if not first_party(path):
            continue
        for check in match[6].split(","):
            if check.startswith("performance-"):
                records.append(
                    {
                        "file": path.relative_to(ROOT).as_posix(),
                        "line": int(match[2]),
                        "column": int(match[3]),
                        "check": check,
                        "message": match[5],
                        "review": (
                            "installed-header"
                            if any(
                                path.is_relative_to(ROOT / area)
                                for area in ("src/cpp/include", "bindings/c/include")
                            )
                            else "implementation"
                        ),
                    }
                )
    return records


def scan(tool, build_dir, source, directory, log_dir):
    command = [
        tool,
        "-p",
        str(build_dir),
        *tidy_arguments(),
        "--quiet",
        "--use-color=false",
        str(source),
    ]
    result = subprocess.run(
        command,
        cwd=directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
    )
    rel = source.relative_to(ROOT).as_posix()
    log_path = log_dir / (rel.replace("/", "__") + ".log")
    log_path.write_text(result.stdout)
    # Keep failures visible even if a tool/version incorrectly returns zero.
    failed = result.returncode != 0 or re.search(
        r"(^|\n)(?:.*?:\d+:\d+: )?(?:fatal )?error:", result.stdout
    )
    return {
        "file": rel,
        "returncode": result.returncode,
        "failed": bool(failed),
        "log": str(log_path),
        "diagnostics": parse_diagnostics(result.stdout, directory),
    }


def positive_int(value):
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-dir", type=pathlib.Path, help="reuse a configured compilation database"
    )
    parser.add_argument(
        "--mesh-backend",
        choices=("MESHIO", "NATIVE", "KRATOS"),
        default="MESHIO",
        help="mesh backend for the default tree (ignored with --build-dir)",
    )
    parser.add_argument(
        "--out", type=pathlib.Path, help="JSON report (plus sibling logs)"
    )
    parser.add_argument("--jobs", type=positive_int, default=2)
    parser.add_argument(
        "--file",
        action="append",
        type=pathlib.Path,
        help="scan only this TU; repeatable",
    )
    args = parser.parse_args(argv)
    tool = os.environ.get("CLANG_TIDY", "clang-tidy")
    build_dir = (
        args.build_dir
        or ROOT / "build" / f"performance-tidy-{args.mesh_backend.lower()}"
    ).resolve()
    database = build_dir / "compile_commands.json"
    out = (args.out or build_dir / "performance-tidy.json").resolve()
    if out in (database, build_dir / "CMakeCache.txt") or any(
        out.is_relative_to(ROOT / area)
        for area in ("src", "bindings", "tests", "tools")
    ):
        parser.error("the report must not overwrite build configuration or source code")
    out.parent.mkdir(parents=True, exist_ok=True)
    # Fail closed even when configuration, database validation or tool startup fails.
    out.write_text(
        json.dumps({"schema_version": 1, "complete": False}, indent=2) + "\n"
    )
    if args.build_dir is None:
        configure(build_dir, args.mesh_backend)
    entries = json.loads(database.read_text())
    eligible = [
        entry
        for entry in entries
        if first_party(source_path(entry))
        and source_path(entry).suffix in (".c", ".cc", ".cpp", ".cxx")
    ]
    sources = {}
    for entry in eligible:
        source = source_path(entry)
        directory = pathlib.Path(entry["directory"]).resolve()
        if not source.is_file() or not directory.is_dir():
            parser.error(f"stale compile command: {source} in {directory}")
        sources[source] = directory
    if args.file:
        wanted = {path.resolve() for path in args.file}
        missing = wanted - sources.keys()
        if missing:
            parser.error("no first-party compile command for: " + str(sorted(missing)))
        sources = {path: sources[path] for path in wanted}
    if not sources:
        parser.error("compilation database contains no first-party C/C++ sources")

    version = subprocess.check_output([tool, "--version"], text=True).strip()
    check_output = subprocess.check_output(
        [tool, *tidy_arguments(), "--list-checks"], text=True, cwd=ROOT
    )
    checks = sorted(set(re.findall(r"^\s+(performance-\S+)$", check_output, re.M)))
    if not checks:
        parser.error("Clang-Tidy has no enabled performance checks")
    log_dir = out.with_name(out.stem + "-logs")
    log_dir.mkdir(exist_ok=True)
    report = {
        "schema_version": 1,
        "complete": False,
        "clang_tidy": version,
        "clang_tidy_config": tidy_config(),
        "checks": checks,
        "build_dir": str(build_dir),
        "compile_commands_sha256": hashlib.sha256(database.read_bytes()).hexdigest(),
        "configuration": configuration(
            build_dir, [entry for entry in eligible if source_path(entry) in sources]
        ),
        "scope": "selected" if args.file else "database",
        "database_compile_commands": len(eligible),
        "database_translation_units": len({source_path(entry) for entry in eligible}),
        "translation_units": [],
        "diagnostics": [],
    }
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"performance-tidy: {len(sources)} TUs, {len(checks)} checks, {args.jobs} jobs",
        flush=True,
    )
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [
            pool.submit(scan, tool, build_dir, source, directory, log_dir)
            for source, directory in sorted(sources.items())
        ]
        results = [future.result() for future in futures]

    diagnostics = {}
    for result in results:
        for diagnostic in result.pop("diagnostics"):
            key = tuple(
                diagnostic[field] for field in ("file", "line", "column", "check")
            )
            if key not in diagnostics:
                diagnostics[key] = {**diagnostic, "translation_units": []}
            origins = diagnostics[key]["translation_units"]
            if result["file"] not in origins:
                origins.append(result["file"])
    failures = [result["file"] for result in results if result["failed"]]
    report.update(
        complete=not failures,
        translation_units=results,
        diagnostics=[diagnostics[key] for key in sorted(diagnostics)],
    )
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(f"performance-tidy: {len(diagnostics)} unique diagnostics; report: {out}")
    if failures:
        print(
            f"performance-tidy: INCOMPLETE: {len(failures)} TUs failed; see {log_dir}",
            file=sys.stderr,
        )
        return 1
    print("performance-tidy: scan complete; findings are advisory, not automatic fixes")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as exc:
        print(f"performance-tidy: {exc}", file=sys.stderr)
        sys.exit(1)
