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


# Writers whose parsers have a region path the plain seed mesh never reaches:
# XDMF <Set>s, the VTU/VTP <FieldData> convention, MDPA sub model parts.
_REGION_SEED_WRITERS = ("mdpa", "xdmf", "vtu", "vtp")


def _bundle_pack(entries) -> bytes:
    """Pack ``[(name, bytes)]`` into the MIOB bundle tests/fuzz/fuzz_bundle.hpp reads."""
    import struct

    out = bytearray(b"MIOB\x01")
    out += struct.pack("<H", len(entries))
    for name, data in entries:
        name_b = name.encode()
        out += struct.pack("<H", len(name_b)) + name_b
        out += struct.pack("<I", len(data)) + data
    return bytes(out)


def _seed_mesh(surface=False, regions=False):
    import meshioplusplus
    from meshioplusplus._regions import Region

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
        mesh = meshioplusplus.Mesh(
            points[:4, :2], [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))]
        )
        if regions:
            mesh.regions = [
                Region("fixed", "point", np.array([0, 3])),
                Region("wall", "cell", np.array([1])),
                Region("edge", "side", np.array([[0, 0], [1, 2]])),
            ]
        return mesh
    cells = [
        ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
        ("tetra", np.array([[0, 1, 3, 4]])),
    ]
    mesh = meshioplusplus.Mesh(
        points,
        cells,
        point_data={"u": np.arange(6, dtype=float)},
        cell_data={"id": [np.array([1, 2]), np.array([3])]},
    )
    if regions:
        mesh.regions = [
            Region("fixed", "point", np.array([0, 1, 4])),
            Region("skin", "cell", np.array([0, 1])),
            Region("core", "cell", np.array([2])),
            Region("face", "side", np.array([[2, 0], [2, 3]])),
        ]
    return mesh


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
    ap.add_argument(
        "--library-readers",
        action="store_true",
        help="include library_formats.txt for an instrumented-library build",
    )
    args = ap.parse_args(argv)

    import meshioplusplus

    readers = native_readers(args.replay)
    wanted = [f for f in args.formats.split(",") if f] or readers
    written = meshioplusplus.formats()["writable"]
    skip = SKIP.copy()
    if args.library_readers:
        skip -= {
            line.strip()
            for line in (pathlib.Path(__file__).parent / "library_formats.txt")
            .read_text()
            .splitlines()
            if line.strip() and not line.startswith("#")
        }
    for fmt in wanted:
        if fmt in skip or fmt not in readers:
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
        if fmt == "vtx":
            # `.bp` directories, not lone files: pack each into one bundle input
            # so the harness reconstructs the directory the reader opens.
            n += _vtx_fixture_bundles(dst, args.max_bytes)
        else:
            for f in found:
                if (
                    f.is_file()
                    and 0 < f.stat().st_size <= args.max_bytes
                    and not f.read_bytes().startswith(
                        b"version https://git-lfs.github.com/spec/v1"
                    )
                ):
                    # A companion payload (an XDMF `.h5`/`.bin` next to its
                    # `.xdmf`) is meaningless alone: bundles below carry it.
                    if fmt == "xdmf" and f.suffix.lower() in (".h5", ".hdf5", ".bin"):
                        continue
                    digest = hashlib.sha1(f.read_bytes()).hexdigest()[:12]
                    shutil.copyfile(f, dst / f"{digest}-{f.name}")
                    n += 1
        writer = _WRITER_NAME.get(fmt, fmt)
        if writer in written and fmt != "xdmf":
            n += _generated(writer, ext, dst, args.max_bytes)
        if fmt == "xdmf":
            n += _xdmf_bundle_seeds(dst, args.max_bytes)
        print(f"{fmt}: {n} seed(s)")
    return 0


def _generated(writer, ext, dst, max_bytes) -> int:
    """Seeds meshio++ writes itself: the mixed mesh, else a 2-D triangle patch.

    A writer with a region path (``_REGION_SEED_WRITERS``) gets a second seed
    from a mesh that carries point, cell and side regions.
    """
    import meshioplusplus

    groups = [(False, "generated")]
    if writer in _REGION_SEED_WRITERS:
        groups.insert(0, (True, "generated-regions"))
    total = 0
    for regions, prefix in groups:
        for surface in (False, True):
            with tempfile.TemporaryDirectory() as tmp:
                path = pathlib.Path(tmp) / f"seed{ext}"
                try:
                    meshioplusplus.write(
                        path, _seed_mesh(surface, regions), file_format=writer
                    )
                except Exception:  # a writer refusing this seed mesh: try the next
                    continue
                for f in sorted(pathlib.Path(tmp).iterdir()):
                    if f.is_file() and f.stat().st_size <= max_bytes:
                        shutil.copyfile(f, dst / f"{prefix}-{f.name}")
                        total += 1
                break
        else:
            print(f"{writer}: no {prefix} seed (the writer refused both seed meshes)")
    return total


def _extension_for(fmt: str) -> str:
    from meshioplusplus._helpers import extension_to_filetypes

    for ext, types in extension_to_filetypes.items():
        if fmt in types:
            return ext
    return ".dat"


def _vtx_fixture_bundles(dst, max_bytes) -> int:
    """Pack each DOLFINx `.bp` fixture directory into one bundle seed."""
    total = 0
    for bp in sorted((MESHES / "vtx").glob("*.bp")):
        if not bp.is_dir():
            continue
        entries = []
        for f in sorted(bp.rglob("*")):
            if not f.is_file():
                continue
            rel = f.relative_to(bp)
            data = f.read_bytes()
            if data.startswith(b"version https://git-lfs.github.com/spec/v1"):
                continue
            entries.append((f"input.bp/{rel.as_posix()}", data))
        if not entries:
            continue
        blob = _bundle_pack(entries)
        if len(blob) <= max_bytes:
            digest = hashlib.sha1(blob).hexdigest()[:12]
            (dst / f"{digest}-{bp.name}.miob").write_bytes(blob)
            total += 1
    return total


def _xdmf_bundle_seeds(dst, max_bytes) -> int:
    """Positive XDMF bundles the production reader reconstructs: XML, Binary and HDF.

    The plain generated seeds above already cover single-file XML; bundles prove
    the companion path (heavy-data files, DataItem references) the roadmap calls
    out, with one bundle per data format plus a region-carrying HDF bundle.
    """
    import meshioplusplus

    total = 0
    for data_format in ("XML", "Binary", "HDF"):
        for regions, prefix in ((False, "bundle"), (True, "bundle-regions")):
            if regions and data_format != "HDF":
                continue
            with tempfile.TemporaryDirectory() as tmp:
                tmp_p = pathlib.Path(tmp)
                try:
                    meshioplusplus.write(
                        tmp_p / "seed.xdmf",
                        _seed_mesh(False, regions),
                        file_format="xdmf",
                        data_format=data_format,
                    )
                except Exception:
                    continue
                entries = []
                for f in sorted(tmp_p.rglob("*")):
                    if f.is_file() and f.stat().st_size <= max_bytes:
                        entries.append((f.name, f.read_bytes()))
                if not any(name.endswith(".xdmf") for name, _ in entries):
                    continue
                blob = _bundle_pack(entries)
                if len(blob) <= max_bytes:
                    (dst / f"{prefix}-{data_format.lower()}.miob").write_bytes(blob)
                    total += 1
    return total


if __name__ == "__main__":
    sys.exit(main())
