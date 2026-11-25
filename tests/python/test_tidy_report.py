"""Unit tests for the generalized Clang-Tidy runner tools/tidy_report.py (roadmap 3.5.1)."""

from __future__ import annotations

import importlib.util
import pathlib

REPO = pathlib.Path(__file__).resolve().parents[2]
_SPEC = importlib.util.spec_from_file_location(
    "tidy_report", REPO / "tools" / "tidy_report.py"
)
assert _SPEC and _SPEC.loader
tr = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(tr)


def test_parse_diagnostics_supports_arbitrary_check_families():
    src = REPO / "src" / "cpp" / "src" / "operations" / "smooth.cpp"
    out = (
        f"{src}:42:7: warning: narrowing conversion [bugprone-narrowing-conversions]\n"
        f"{src}:99:3: warning: unnecessary copy [performance-unnecessary-value-param]\n"
        f"/usr/include/vector:10:1: warning: system header [bugprone-foo]\n"
    )
    bugprone = tr.parse_diagnostics(
        out, REPO, enabled_checks={"bugprone-narrowing-conversions"}
    )
    assert len(bugprone) == 1
    assert bugprone[0]["check"] == "bugprone-narrowing-conversions"
    assert bugprone[0]["file"] == "src/cpp/src/operations/smooth.cpp"
    assert bugprone[0]["review"] == "implementation"

    both = tr.parse_diagnostics(
        out,
        REPO,
        enabled_checks={
            "bugprone-narrowing-conversions",
            "performance-unnecessary-value-param",
        },
    )
    assert len(both) == 2


def test_compare_baseline_ignores_line_numbers_and_flags_new_findings():
    base = [
        {
            "file": "src/cpp/src/a.cpp",
            "line": 10,
            "check": "bugprone-branch-clone",
            "message": "repeated branch",
        },
        {
            "file": "src/cpp/src/b.cpp",
            "line": 20,
            "check": "bugprone-use-after-move",
            "message": "used after move",
        },
    ]
    # First finding moved lines (still baselined); second resolved; third is new.
    cur = [
        {
            "file": "src/cpp/src/a.cpp",
            "line": 55,
            "check": "bugprone-branch-clone",
            "message": "repeated branch",
        },
        {
            "file": "src/cpp/src/c.cpp",
            "line": 12,
            "check": "bugprone-sizeof-expression",
            "message": "suspicious sizeof",
        },
    ]
    new_f, resolved = tr.compare_baseline(cur, base)
    assert [d["file"] for d in new_f] == ["src/cpp/src/c.cpp"]
    assert [d["file"] for d in resolved] == ["src/cpp/src/b.cpp"]
