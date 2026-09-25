"""Build per-format seed corpora for the reader fuzz target (doc/fuzzing.md).

Two sources, both small:

* every fixture under ``tests/python/meshes/<format>/`` whose directory is
  named after a native registry reader, capped at ``--max-bytes``;
* one file per writable format, written from a small mixed mesh by
  meshio++ itself, so formats without a fixture directory still start from a
  valid file instead of from nothing.

Usage::

    python tools/fuzz/seed_corpus.py OUT_DIR --replay BUILD/meshioplusplus_fuzz_replay \
        [--formats gmsh,vtk] [--max-bytes 262144]

writes ``OUT_DIR/<format>/<name>`` and prints one line per format.
"""

from __future__ import annotations

import argparse
import hashlib
import pathlib
import shutil
import subprocess
import sys
import tempfile

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
MESHES = REPO / "tests" / "python" / "meshes"

# Readers the campaign does not fuzz, with the reasons, in not_fuzzed.txt.
SKIP = {
    line.strip()
    for line in (pathlib.Path(__file__).parent / "not_fuzzed.txt")
    .read_text()
    .splitlines()
    if line.strip() and not line.startswith("#")
}


# Registry reader name -> the Python writer's name, where they differ.
_WRITER_NAME = {"dolfin": "dolfin-xml", "ansysinp": "ansysInp"}


def _seed_mesh(surface=False):
    import meshioplusplus

    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
        ]
    )
    if surface:
        return meshioplusplus.Mesh(
            points[:4, :2], [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))]
        )
    cells = [
        ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
        ("tetra", np.array([[0, 1, 3, 4]])),
    ]
    return meshioplusplus.Mesh(
        points,
        cells,
        point_data={"u": np.arange(6, dtype=float)},
        cell_data={"id": [np.array([1, 2]), np.array([3])]},
    )


def native_readers(replay: pathlib.Path) -> list[str]:
    """The readers compiled into the harness's build (``-list-formats``)."""
    out = subprocess.run(
        [str(replay), "-list-formats"], check=True, capture_output=True, text=True
    ).stdout
    return sorted(line.strip() for line in out.splitlines() if line.strip())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", type=pathlib.Path)
    ap.add_argument(
        "--replay",
        type=pathlib.Path,
        required=True,
        help="the meshioplusplus_fuzz_replay binary of the build to fuzz",
    )
    ap.add_argument("--formats", default="", help="comma-separated subset")
    ap.add_argument("--max-bytes", type=int, default=256 * 1024)
    args = ap.parse_args(argv)

    import meshioplusplus

    readers = native_readers(args.replay)
    wanted = [f for f in args.formats.split(",") if f] or readers
    written = meshioplusplus.formats()["writable"]
    for fmt in wanted:
        if fmt in SKIP or fmt not in readers:
            continue
        dst = args.out / fmt
        dst.mkdir(parents=True, exist_ok=True)
        n = 0
        ext = _extension_for(_WRITER_NAME.get(fmt, fmt))
        src = MESHES / fmt
        found = sorted(src.rglob("*")) if src.is_dir() else []
        # Files of the format filed under another directory (a .t19 under
        # marc/, a .vti under vtk/), found by the format's extension.
        if ext != ".dat":
            found += sorted(p for p in MESHES.rglob(f"*{ext}") if p.parent.name != fmt)[
                :20
            ]
        for f in found:
            if f.is_file() and 0 < f.stat().st_size <= args.max_bytes:
                digest = hashlib.sha1(f.read_bytes()).hexdigest()[:12]
                shutil.copyfile(f, dst / f"{digest}-{f.name}")
                n += 1
        writer = _WRITER_NAME.get(fmt, fmt)
        if writer in written:
            n += _generated(writer, ext, dst, args.max_bytes)
        print(f"{fmt}: {n} seed(s)")
    return 0


def _generated(writer, ext, dst, max_bytes) -> int:
    """Seeds meshio++ writes itself: the mixed mesh, else a 2-D triangle patch."""
    import meshioplusplus

    for surface in (False, True):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / f"seed{ext}"
            try:
                meshioplusplus.write(path, _seed_mesh(surface), file_format=writer)
            except Exception:  # a writer refusing this seed mesh: try the next
                continue
            n = 0
            for f in sorted(pathlib.Path(tmp).iterdir()):
                if f.is_file() and f.stat().st_size <= max_bytes:
                    shutil.copyfile(f, dst / f"generated-{f.name}")
                    n += 1
            return n
    print(f"{writer}: no generated seed (the writer refused both seed meshes)")
    return 0


def _extension_for(fmt: str) -> str:
    from meshioplusplus._helpers import extension_to_filetypes

    for ext, types in extension_to_filetypes.items():
        if fmt in types:
            return ext
    return ".dat"


if __name__ == "__main__":
    sys.exit(main())
