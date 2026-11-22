"""Write ``tests/python/meshes/MANIFEST.json``, the list of every reference fixture.

Each fixture under ``tests/python/meshes/`` and ``tests/python/input/`` gets a row
with its size, its SHA-256 and the README that documents its directory. A fixture
that is a Git LFS pointer (a checkout without ``git lfs pull``) is recorded with the
``oid`` and ``size`` the pointer carries, which are the hash and size of the real
file, so the manifest is the same with or without the LFS content::

    python tools/gen_corpus_manifest.py            # rewrite the manifest
    python tools/gen_corpus_manifest.py --check    # exit 1 if it is stale

The producer tool, the licence and the source URL of third-party fixtures are
kept in the per-directory ``README.md`` files (prose); ``test_corpus_manifest.py``
checks that every fixture is listed and every row exists, so a fixture cannot be
added without it appearing here and in the directory's README tables.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
ROOT = REPO / "tests" / "python"
SECTIONS = ("meshes", "input")
MANIFEST = ROOT / "meshes" / "MANIFEST.json"
_POINTER_RE = re.compile(
    rb"^version https://git-lfs\.github\.com/spec/v1\noid sha256:([0-9a-f]{64})\nsize (\d+)\n"
)
_SKIP_NAMES = {"MANIFEST.json", ".gitattributes", ".DS_Store"}


def fixture_paths():
    """Every fixture file as a sorted list of repository-relative POSIX paths."""
    found = []
    for section in SECTIONS:
        base = ROOT / section
        if not base.is_dir():
            continue
        for path in base.rglob("*"):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            if path.name in _SKIP_NAMES or path.name.upper().startswith("README"):
                continue
            found.append(path.relative_to(REPO).as_posix())
    return sorted(found)


def _digest(path):
    data = path.read_bytes()
    match = _POINTER_RE.match(data)
    if match:
        return match.group(1).decode(), int(match.group(2)), True
    return hashlib.sha256(data).hexdigest(), len(data), False


def _readme_for(path):
    """The README closest above the fixture, within its section, or ``None``."""
    section_root = ROOT / pathlib.PurePosixPath(path).parts[2]
    folder = (REPO / path).parent
    while True:
        for name in ("README.md", "README"):
            candidate = folder / name
            if candidate.is_file():
                return candidate.relative_to(REPO).as_posix()
        if folder == section_root or folder == REPO:
            return None
        folder = folder.parent


def build():
    rows = []
    for rel in fixture_paths():
        sha, size, _pointer = _digest(REPO / rel)
        rows.append(
            {"path": rel, "bytes": size, "sha256": sha, "readme": _readme_for(rel)}
        )
    return {"schema_version": 1, "fixtures": rows}


def render(doc):
    return json.dumps(doc, indent=1) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check", action="store_true", help="exit 1 if the manifest is stale"
    )
    args = parser.parse_args(argv)
    text = render(build())
    if args.check:
        current = MANIFEST.read_text(encoding="utf-8") if MANIFEST.is_file() else ""
        if current != text:
            print(
                f"{MANIFEST.relative_to(REPO)} is stale: run tools/gen_corpus_manifest.py"
            )
            return 1
        return 0
    MANIFEST.write_text(text, encoding="utf-8")
    print(
        f"wrote {MANIFEST.relative_to(REPO)} ({len(json.loads(text)['fixtures'])} fixtures)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
