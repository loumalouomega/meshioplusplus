"""Shared pytest configuration.

Reference fixtures under ``meshes/`` and ``input/`` live in Git LFS. A checkout
that skipped the fetch (CI uses ``lfs: false``) hands a reader a three-line
pointer file, which then fails as a parse error of the format under test. The
check below turns that into one message that names the command to run.
"""

import pathlib
import re

import pytest

_HERE = pathlib.Path(__file__).parent
_FIXTURE_ROOTS = ("meshes", "input")
_POINTER_MAGIC = b"version https://git-lfs.github.com/spec/v1"


def _is_pointer(path):
    try:
        with open(path, "rb") as f:
            return f.read(len(_POINTER_MAGIC)) == _POINTER_MAGIC
    except OSError:
        return False


def _pointer_directories():
    """Map ``<root>/<dir>`` (relative to this folder) to its LFS pointer files."""
    found = {}
    for root in _FIXTURE_ROOTS:
        base = _HERE / root
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file() and _is_pointer(path):
                top = path.relative_to(base).parts[0]
                found.setdefault(f"{root}/{top}", []).append(path)
    return found


def _modules_text(session):
    texts = []
    for path in {pathlib.Path(str(item.fspath)) for item in session.items}:
        try:
            texts.append(path.read_text(encoding="utf-8", errors="replace"))
        except OSError:
            continue
    return texts


def pytest_collection_finish(session):
    """Stop with one message when a collected module reads an unfetched fixture."""
    if not session.items:
        return
    pointers = _pointer_directories()
    if not pointers:
        return
    texts = _modules_text(session)
    hit = {}
    for key, files in pointers.items():
        root, name = key.split("/", 1)
        # Only a path built from the fixture root and the directory counts: a bare
        # "vtu" is a format name in most test modules, and must not stop the run.
        pattern = re.compile(
            r"""["']%(r)s["']\s*[/,]\s*["']%(n)s["']|\b%(r)s/%(n)s\b"""
            % {"r": re.escape(root), "n": re.escape(name)}
        )
        if any(pattern.search(t) for t in texts):
            hit[key] = files
    if not hit:
        return
    includes = ",".join(f"tests/python/{key}" for key in sorted(hit))
    count = sum(len(v) for v in hit.values())
    pytest.exit(
        f"{count} test fixture(s) are Git LFS pointers, not files "
        f"(in {', '.join(sorted(hit))}). Fetch them with:\n"
        f"    git lfs pull --include={includes}",
        returncode=4,
    )
