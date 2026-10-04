"""Performance audits are advisory, but an incomplete scan must never pass."""

import csv
import importlib.util
import itertools
import json
import pathlib
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]


@pytest.fixture
def audit(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "performance_tidy", ROOT / "tools/performance_tidy.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    root = tmp_path / "repo with spaces"
    build = root / "build"
    build.mkdir(parents=True)
    monkeypatch.setattr(module, "ROOT", root)
    sources = [
        root / "src/cpp/src/first.cpp",
        root / "bindings/python/second.cpp",
        root / "src/cpp/third_party/vendor.cpp",
        root / "src/single_include/generated.cpp",
    ]
    header = root / "src/cpp/include/meshioplusplus/header.hpp"
    for path in [*sources, header]:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("// unchanged\n")
    entries = [
        {
            "directory": str(root),
            "file": str(path.relative_to(root)),
            "arguments": ["clang++", "-DMESHIOPLUSPLUS_PARALLEL_SEQ", str(path)],
        }
        for path in sources
    ]
    entries.append(entries[0].copy())  # Multiple compile commands for the same TU.
    (build / "compile_commands.json").write_text(json.dumps(entries))
    (build / "CMakeCache.txt").write_text(
        "MESHIOPLUSPLUS_WITH_HDF5:BOOL=ON\nCMAKE_CXX_COMPILER:FILEPATH=clang++\n"
    )
    calls = []

    def check_output(command, **kwargs):
        if "--version" in command:
            return "LLVM version test\n"
        assert "--list-checks" in command
        return (
            "Enabled checks:\n"
            "    performance-unnecessary-copy-initialization\n"
            "    performance-enum-size\n"
        )

    def run(command, **kwargs):
        calls.append(command)
        output = (
            f"{header}:12:3: warning: expensive copy "
            "[performance-unnecessary-copy-initialization]\n"
            "    source excerpt\n"
            f"{root / 'src/cpp/third_party/vendor.hpp'}:1:1: warning: vendor "
            "[performance-enum-size]\n"
        )
        return subprocess.CompletedProcess(command, 0, stdout=output)

    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    monkeypatch.setattr(module.subprocess, "run", run)
    monkeypatch.setenv("CLANG_TIDY", "fixture-clang-tidy")
    out = root / "report.json"
    return module, root, build, sources, out, calls


def test_complete_advisory_report_deduplicates_headers(audit):
    module, root, build, sources, out, calls = audit
    assert module.main(["--build-dir", str(build), "--out", str(out)]) == 0
    report = json.loads(out.read_text())
    assert report["complete"]
    assert report["scope"] == "database"
    assert report["database_compile_commands"] == 3
    assert report["database_translation_units"] == 2
    assert report["clang_tidy_config"]["WarningsAsErrors"] == ""
    assert len(report["translation_units"]) == 2
    assert len(calls) == 2
    assert len(report["diagnostics"]) == 1
    diagnostic = report["diagnostics"][0]
    assert diagnostic["review"] == "installed-header"
    assert len(diagnostic["translation_units"]) == 2
    assert len(report["compile_commands_sha256"]) == 64
    config = report["configuration"]
    assert config["cmake_cache"]["MESHIOPLUSPLUS_WITH_HDF5"] == "ON"
    assert config["compile_defines"] == ["MESHIOPLUSPLUS_PARALLEL_SEQ"]
    for call in calls:
        assert call[0] == "fixture-clang-tidy"
        assert "--checks=-*,performance-*" in call
        config = json.loads(
            next(arg[9:] for arg in call if arg.startswith("--config="))
        )
        assert config == report["clang_tidy_config"]
        assert not any(arg.startswith("--fix") for arg in call)
        assert pathlib.Path(call[-1]) in sources[:2]
    for source in sources:
        assert source.read_text() == "// unchanged\n"
    for result in report["translation_units"]:
        assert pathlib.Path(result["log"]).is_file()


@pytest.mark.parametrize("returncode", [0, 1])
def test_parse_errors_fail_even_when_the_tool_returns_zero(
    audit, monkeypatch, returncode
):
    module, root, build, sources, out, calls = audit

    def run(command, **kwargs):
        return subprocess.CompletedProcess(
            command,
            returncode,
            stdout="error: unknown argument: '-fmodules-ts' [clang-diagnostic-error]\n",
        )

    monkeypatch.setattr(module.subprocess, "run", run)
    assert module.main(["--build-dir", str(build), "--out", str(out)]) == 1
    report = json.loads(out.read_text())
    assert not report["complete"]
    assert all(item["failed"] for item in report["translation_units"])
    assert (
        "unknown argument"
        in pathlib.Path(report["translation_units"][0]["log"]).read_text()
    )


def test_nonzero_exit_without_diagnostics_fails(audit, monkeypatch):
    module, root, build, sources, out, calls = audit
    monkeypatch.setattr(
        module.subprocess,
        "run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, -9, stdout=""),
    )
    assert module.main(["--build-dir", str(build), "--out", str(out)]) == 1
    assert not json.loads(out.read_text())["complete"]


def test_interrupted_run_invalidates_an_earlier_success(audit, monkeypatch):
    module, root, build, sources, out, calls = audit
    out.write_text('{"complete": true}')

    def run(command, **kwargs):
        raise OSError("fixture tool cannot start")

    monkeypatch.setattr(module.subprocess, "run", run)
    with pytest.raises(OSError, match="cannot start"):
        module.main(["--build-dir", str(build), "--out", str(out)])
    assert not json.loads(out.read_text())["complete"]


@pytest.mark.parametrize("stage", ["database", "version"])
def test_preflight_failure_invalidates_an_earlier_success(audit, monkeypatch, stage):
    module, root, build, sources, out, calls = audit
    out.write_text('{"complete": true}')
    if stage == "database":
        (build / "compile_commands.json").write_text("malformed JSON")
    else:

        def check_output(*args, **kwargs):
            raise OSError("fixture missing clang-tidy")

        monkeypatch.setattr(module.subprocess, "check_output", check_output)
    with pytest.raises((OSError, ValueError)):
        module.main(["--build-dir", str(build), "--out", str(out)])
    assert not json.loads(out.read_text())["complete"]
    assert not calls


@pytest.mark.parametrize("target", ["source", "database"])
def test_report_cannot_overwrite_code_or_compile_commands(audit, target):
    module, root, build, sources, out, calls = audit
    path = sources[0] if target == "source" else build / "compile_commands.json"
    previous = path.read_bytes()
    with pytest.raises(SystemExit, match="2"):
        module.main(["--build-dir", str(build), "--out", str(path)])
    assert path.read_bytes() == previous
    assert not calls


def test_selected_file_is_marked_as_partial_coverage(audit):
    module, root, build, sources, out, calls = audit
    assert (
        module.main(
            ["--build-dir", str(build), "--out", str(out), "--file", str(sources[0])]
        )
        == 0
    )
    report = json.loads(out.read_text())
    assert report["complete"]  # Complete within the explicitly selected scope.
    assert report["scope"] == "selected"
    assert report["database_translation_units"] == 2
    assert len(report["translation_units"]) == 1


@pytest.mark.parametrize("index", [2, 3])
def test_vendored_and_generated_sources_cannot_be_selected(audit, index):
    module, root, build, sources, out, calls = audit
    with pytest.raises(SystemExit, match="2"):
        module.main(["--build-dir", str(build), "--file", str(sources[index])])
    assert not calls


def test_stale_compile_command_is_not_silently_skipped(audit):
    module, root, build, sources, out, calls = audit
    sources[0].unlink()
    with pytest.raises(SystemExit, match="2"):
        module.main(["--build-dir", str(build)])
    assert not calls


def test_empty_database_fails(audit):
    module, root, build, sources, out, calls = audit
    (build / "compile_commands.json").write_text("[]")
    with pytest.raises(SystemExit, match="2"):
        module.main(["--build-dir", str(build)])
    assert not calls


def test_no_performance_checks_fails(audit, monkeypatch):
    module, root, build, sources, out, calls = audit
    monkeypatch.setattr(
        module.subprocess, "check_output", lambda *a, **kw: "no checks\n"
    )
    with pytest.raises(SystemExit, match="2"):
        module.main(["--build-dir", str(build)])
    assert not calls


def test_multiple_check_names_and_relative_diagnostic_paths(audit):
    module, root, build, sources, out, calls = audit
    records = module.parse_diagnostics(
        "src/cpp/src/first.cpp:42:5: warning: single character "
        "[performance-faster-string-find,performance-prefer-single-char-overloads]",
        root,
    )
    assert len(records) == 2
    assert {record["check"] for record in records} == {
        "performance-faster-string-find",
        "performance-prefer-single-char-overloads",
    }
    assert all(record["review"] == "implementation" for record in records)


def test_installed_c_header_is_not_classified_as_an_implementation(audit):
    module, root, build, sources, out, calls = audit
    records = module.parse_diagnostics(
        "bindings/c/include/meshioplusplus/meshioplusplus.h:42:5: warning: enum width "
        "[performance-enum-size]",
        root,
    )
    assert len(records) == 1
    assert records[0]["review"] == "installed-header"


def test_default_configuration_is_separate_and_clang_compatible(audit):
    module, root, build, sources, out, calls = audit
    module.configure(build, "NATIVE")
    command = calls[0]
    assert command[0] == "cmake"
    assert "-DCMAKE_CXX_SCAN_FOR_MODULES=OFF" in command
    assert "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON" in command
    assert "-DMESHIOPLUSPLUS_MESH_BACKEND=NATIVE" in command
    assert "-DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ" in command
    assert "-DMESHIOPLUSPLUS_BUILD_PYTHON=OFF" in command


def test_fluent_batch_evidence_has_a_complete_matching_digest_matrix():
    """Validate the recorded evidence, not a noisy timing threshold for CI."""
    with (ROOT / "benchmark/performance_tidy_fluent.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    configs = [("seq", "1"), *itertools.product(("openmp", "tbb"), ("1", "4", "8"))]
    cases = list(itertools.product(("M", "L"), ("uniform", "mixed"), ("False", "True")))
    expected = {
        (stage, backend, threads, *case)
        for stage in ("before", "after")
        for backend, threads in configs
        for case in cases
    }
    actual, digests = {}, {}
    for row in rows:
        case = tuple(row[field] for field in ("size", "kind", "binary"))
        key = (row["stage"], row["backend"], row["threads"], *case)
        assert key not in actual
        actual[key] = row
        digest = row["sha256"]
        assert len(digest) == 64 and all(c in "0123456789abcdef" for c in digest)
        assert digests.setdefault(case, digest) == digest
    assert actual.keys() == expected
    assert len(set(digests.values())) == 8
    reductions = {
        ("M", "uniform"): (18, 1068148),
        ("L", "uniform"): (20, 4395604),
        ("M", "mixed"): (38, 21951300),
        ("L", "mixed"): (40, 30754276),
    }
    for backend, threads in configs:
        for (size, kind), saved in reductions.items():
            for binary in ("False", "True"):
                before = actual[("before", backend, threads, size, kind, binary)]
                after = actual[("after", backend, threads, size, kind, binary)]
                assert (
                    tuple(
                        int(before[field]) - int(after[field])
                        for field in ("allocations", "requested_bytes")
                    )
                    == saved
                )
