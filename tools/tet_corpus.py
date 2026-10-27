"""Tetrahedralize a seeded sample of Thingi10K and measure it (roadmap 6.7.8).

Thingi10K (Q. Zhou and A. Jacobson, 2016) is ten thousand triangle meshes
from a 3-D printing repository with every defect real surfaces have, and the
corpus tetrahedral meshers publish their success rates on. Each model carries
its own licence, so the manifest (``tools/tet_corpus/manifest.json``) names
the dataset's pinned revision, a seed and a sample size, and records ids and
SHA-256 hashes, never the files. Subcommands:

``fetch``
    download the dataset's metadata, draw the sample (``count`` ids, by the
    seed, from the closed, manifold, oriented STL models the metadata lists),
    and extract those models alone from the dataset's archive (streamed, or
    ``--archive`` for a local copy) into ``--dest``; a model whose hash the
    manifest has must match it, and the hashes it lacks are written to
    ``--record`` (commit them).
``run``
    tetrahedralize each fetched model in its own process under a time limit,
    so a crash or a hang is a result that names the model, and measure it:
    status, tetrahedra, whether the boundary was kept, the dihedral extremes
    and the time. The run exits non-zero when a model crashed or hung.
``summary``
    a Markdown summary of a run: the share meshed, the share with the
    boundary kept, the dihedral extremes and the time, per backend.

Only the standard library, NumPy and meshioplusplus are used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import pathlib
import random
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "tools" / "tet_corpus" / "manifest.json"


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"


def finite(value):
    """NaN and infinities as null, so a run is strict JSON."""
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {k: finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite(v) for v in value]
    return value


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_manifest(path):
    manifest = json.loads(pathlib.Path(path).read_text(encoding="utf-8"))
    if manifest.get("Version") != 1:
        raise SystemExit(f"{path}: 'Version' must be 1")
    for key in ("source", "sample"):
        if key not in manifest:
            raise SystemExit(f"{path}: '{key}' is required")
    manifest.setdefault("models", {})
    return manifest


def _open(url, timeout=600):
    return urllib.request.urlopen(url, timeout=timeout)


def _csv_rows(url):
    with _open(url) as r:
        return list(csv.DictReader(io.TextIOWrapper(r, encoding="utf-8")))


def _flag(row, key):
    try:
        return int(float(row.get(key) or 0))
    except ValueError:
        return 0


def eligible(geometry, summary):
    """The ids of the closed, manifold, oriented STL models, ascending."""
    ext = {}
    for row in summary:
        link = (row.get("Link") or "").strip()
        if row.get("ID") and link:
            ext[int(row["ID"])] = link.rsplit(".", 1)[-1].lower()
    ids = []
    for row in geometry:
        if not row.get("file_id"):
            continue
        fid = int(row["file_id"])
        if ext.get(fid) != "stl":
            continue
        if (
            _flag(row, "num_boundary_edges") == 0
            and _flag(row, "edge_manifold") == 1
            and _flag(row, "vertex_manifold") == 1
            and _flag(row, "oriented") == 1
        ):
            ids.append(fid)
    return sorted(ids)


def draw(ids, seed, count):
    """`count` of `ids` drawn by `seed`, ascending."""
    ids = sorted(ids)
    return sorted(random.Random(seed).sample(ids, min(count, len(ids))))


def fetch(args):
    manifest = load_manifest(args.manifest)
    src, sample = manifest["source"], manifest["sample"]
    base = args.base or src["base"]
    geometry = _csv_rows(f"{base}/{src['geometry']}")
    summary = _csv_rows(f"{base}/{src['summary']}")
    ids = draw(eligible(geometry, summary), sample["seed"], sample["count"])
    print(f"tet_corpus: {len(ids)} models drawn", file=sys.stderr)
    dest = pathlib.Path(args.dest)
    dest.mkdir(parents=True, exist_ok=True)
    wanted = {f"{src['members']}/{i}.stl": i for i in ids}
    archive = args.archive or f"{base}/{src['archive']}"
    found = 0
    stream = (
        open(archive, "rb")
        if pathlib.Path(str(archive)).exists()
        else _open(archive, timeout=3600)
    )
    with stream, tarfile.open(fileobj=stream, mode="r|gz") as tar:
        for member in tar:
            fid = wanted.get(member.name)
            if fid is None or not member.isfile():
                continue
            with tar.extractfile(member) as f, open(dest / f"{fid}.stl", "wb") as out:
                shutil.copyfileobj(f, out)
            found += 1
            if found == len(wanted):
                break
    (dest / "sample.json").write_text(canonical({"ids": ids}), encoding="utf-8")
    known, record, bad = manifest["models"], {}, []
    for fid in ids:
        path = dest / f"{fid}.stl"
        if not path.exists():
            print(f"tet_corpus: {fid}: not in the archive", file=sys.stderr)
            continue
        digest = sha256(path)
        if str(fid) in known and known[str(fid)] != digest:
            bad.append(fid)
            path.unlink()
        elif str(fid) not in known:
            record[str(fid)] = digest
    for fid in bad:
        print(f"tet_corpus: {fid}: hash mismatch, removed", file=sys.stderr)
    if record and args.record:
        pathlib.Path(args.record).write_text(canonical(record), encoding="utf-8")
        print(f"tet_corpus: {len(record)} hashes recorded", file=sys.stderr)
    return 1 if bad else 0


def measure_one(path, backend, max_cells):
    """Read and tetrahedralize one model; runs in its own process."""
    import numpy as np

    import meshioplusplus as mio

    out = {"status": "ok"}
    try:
        mesh = mio.read(path, file_format="stl")
        tri = [c.data for c in mesh.cells if c.type == "triangle"]
        surface = mio.Mesh(mesh.points, [("triangle", np.concatenate(tri))])
        surface = mio.clean(surface, weld=True)
    except Exception as exc:  # noqa: BLE001 - every failure is a result here
        return {"status": "read_failed", "error": f"{type(exc).__name__}: {exc}"[:500]}
    out["num_triangles"] = int(sum(len(c.data) for c in surface.cells))
    t0 = time.perf_counter()
    try:
        _, report = mio.tetrahedralize(
            surface, backend=backend, max_cells=max_cells, return_report=True
        )
    except ValueError as exc:  # the surface is refused by name
        out.update(status="refused", error=str(exc)[:500])
        return out
    except Exception as exc:  # noqa: BLE001
        out.update(status="failed", error=f"{type(exc).__name__}: {exc}"[:500])
        return out
    out["mesh_s"] = round(time.perf_counter() - t0, 3)
    for key in (
        "num_tetrahedra",
        "num_regions",
        "boundary_preserved",
        "hausdorff_to_surface",
        "min_dihedral_deg",
        "max_dihedral_deg",
        "backend",
    ):
        out[key] = report[key]
    return out


def run(args):
    import meshioplusplus as mio

    dest = pathlib.Path(args.dest)
    ids = json.loads((dest / "sample.json").read_text(encoding="utf-8"))["ids"]
    if args.models:
        only = {int(i) for i in args.models.split(",")}
        ids = [i for i in ids if i in only]
    results, hung = [], []
    for fid in ids:
        path = dest / f"{fid}.stl"
        row = {"id": fid}
        if not path.exists():
            row["status"] = "unavailable"
            results.append(row)
            continue
        cmd = [
            sys.executable,
            str(pathlib.Path(__file__).resolve()),
            "_one",
            str(path),
            args.backend,
            str(args.max_cells),
        ]
        t0 = time.perf_counter()
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=args.timeout
            )
            lines = proc.stdout.strip().splitlines()
            if proc.returncode == 0 and lines:
                row.update(json.loads(lines[-1]))
            else:
                tail = (proc.stderr or "").strip().splitlines()[-3:]
                row.update(
                    status="crashed",
                    error=f"exit {proc.returncode}: {' '.join(tail)}"[:500],
                )
        except subprocess.TimeoutExpired:
            row.update(status="timeout", error=f"over {args.timeout} s")
        row["wall_s"] = round(time.perf_counter() - t0, 3)
        if row["status"] in ("crashed", "timeout"):
            hung.append(fid)
        results.append(row)
        print(f"tet_corpus: {fid}: {row['status']}, {row['wall_s']} s", file=sys.stderr)
    doc = {
        "Version": 1,
        "meshioplusplus": mio.__version__,
        "settings": {
            "backend": args.backend,
            "max_cells": args.max_cells,
            "timeout_s": args.timeout,
        },
        "models": results,
    }
    pathlib.Path(args.out).write_text(canonical(finite(doc)), encoding="utf-8")
    for fid in hung:
        print(f"tet_corpus: model {fid} crashed or hung", file=sys.stderr)
    return 1 if hung else 0


def render(doc):
    """The Markdown summary of a run."""
    rows = [r for r in doc["models"] if r["status"] != "unavailable"]
    accepted = [r for r in rows if r["status"] not in ("read_failed", "refused")]
    ok = [r for r in accepted if r["status"] == "ok"]
    kept = [r for r in ok if r.get("boundary_preserved")]

    def share(part, whole):
        return f"{len(part)} of {len(whole)}" + (
            f" ({100.0 * len(part) / len(whole):.1f}%)" if whole else ""
        )

    lines = [
        f"### Thingi10K tetrahedralize ({doc['settings']['backend']}, "
        f"meshio++ {doc['meshioplusplus']})",
        "",
        "| measure | value |",
        "|---|---|",
        f"| models drawn | {len(doc['models'])} |",
        f"| accepted (closed, manifold, oriented as read) | {share(accepted, rows)} |",
        f"| meshed | {share(ok, accepted)} |",
        f"| boundary kept | {share(kept, accepted)} |",
    ]
    if ok:
        lo = [r["min_dihedral_deg"] for r in ok]
        hi = [r["max_dihedral_deg"] for r in ok]
        lines += [
            f"| least dihedral, worst / median model | {min(lo):.2f} / "
            f"{sorted(lo)[len(lo) // 2]:.2f} deg |",
            f"| greatest dihedral, worst / median model | {max(hi):.2f} / "
            f"{sorted(hi)[len(hi) // 2]:.2f} deg |",
            f"| tetrahedra | {sum(r['num_tetrahedra'] for r in ok):,} |",
            f"| time, meshed models | {sum(r.get('mesh_s', 0.0) for r in ok):.1f} s |",
        ]
    by = {}
    for r in rows:
        by[r["status"]] = by.get(r["status"], 0) + 1
    lines += ["", "Status: " + ", ".join(f"{k} {v}" for k, v in sorted(by.items()))]
    bad = [r for r in rows if r["status"] in ("crashed", "timeout", "failed")]
    if bad:
        lines += ["", "| model | status | error |", "|---|---|---|"]
        lines += [
            f"| {r['id']} | {r['status']} | {r.get('error', '')[:120]} |" for r in bad
        ]
    return "\n".join(lines) + "\n"


def summary(args):
    doc = json.loads(pathlib.Path(args.results).read_text(encoding="utf-8"))
    sys.stdout.write(render(doc))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("fetch", help="draw the sample and extract its models")
    p.add_argument("--manifest", default=str(MANIFEST))
    p.add_argument("--dest", required=True)
    p.add_argument("--base", help="the dataset's URL in place of the manifest's")
    p.add_argument("--archive", help="a local copy (or URL) of the dataset archive")
    p.add_argument("--record", help="write the hashes the manifest lacks here")
    p = sub.add_parser("run", help="tetrahedralize and measure every fetched model")
    p.add_argument("--dest", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--backend", default="auto")
    p.add_argument("--max-cells", type=int, default=5_000_000)
    p.add_argument("--timeout", type=float, default=300.0)
    p.add_argument("--models", help="a comma-separated subset of ids")
    p = sub.add_parser("summary", help="a Markdown summary of a run")
    p.add_argument("results")
    if argv is None:
        argv = sys.argv[1:]
    if argv and argv[0] == "_one":
        row = measure_one(argv[1], argv[2], int(argv[3]))
        print(json.dumps(finite(row), allow_nan=False, default=float))
        return 0
    args = parser.parse_args(argv)
    return {"fetch": fetch, "run": run, "summary": summary}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
