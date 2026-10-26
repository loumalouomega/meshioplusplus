"""Mesh a corpus of real CAD models and measure the result (roadmap 6.5).

The corpus is a manifest of public STEP files, never the files: each model
is named by its source (a git repository at a pinned commit, or an https
URL) and a SHA-256, and fetched where the run happens (the scheduled
``cad_corpus`` workflow). Subcommands:

``fetch``
    download every source into ``--dest`` and check the hashes; a model whose
    manifest entry has no hash yet gets one recorded in ``--record`` (commit
    it), and a source that cannot be reached is reported, not fatal.
``run``
    read each model (OpenCASCADE), mesh it with ``mesh_cad_surface`` at a
    size relative to its diagonal, and measure it: faces meshed and failed by
    status and by class (``classify``), open, non-manifold and misoriented
    edges, ``max_edge_gap``, the angle extremes, ``compute_quality``'s
    triangle metrics, the Hausdorff distance to the read's own ``BRepMesh``
    triangulation, points, time and peak memory. Each model runs in its own
    process under a time limit, so a crash or a hang is a result, not the end
    of the run.
``compare``
    exit non-zero when a model that meshed in the baseline no longer does,
    naming it (6.5.5's probe).
``baseline``
    write the stable part of a run (no times) as the next baseline.
``table``
    regenerate the corpus section of ``doc/benchmarks.md`` from the baseline;
    ``--check`` fails when it is stale.
``summary``
    a Markdown summary of a run, for the job summary.

Only the standard library, NumPy and meshioplusplus are used.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import pathlib
import random
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "tools" / "cad_corpus" / "manifest.json"
BASELINE = ROOT / "tools" / "cad_corpus" / "baseline.json"
DOC = ROOT / "doc" / "benchmarks.md"
BEGIN = "<!-- cad-corpus:begin (tools/cad_corpus.py table) -->"
END = "<!-- cad-corpus:end -->"

# The failure classes of roadmap 6.5.6, in report order.
CLASSES = (
    "gap_beyond_tolerance",
    "trim_outside_patch",
    "self_intersecting_loop",
    "open_loop",
    "degenerate_patch",
    "budget",
    "other",
)


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
    ids = set()
    for source in manifest["sources"]:
        if source.get("kind") not in ("git", "https"):
            raise SystemExit(
                f"{path}: source {source.get('id')!r}: kind is git or https"
            )
        for model in source["models"]:
            if model["id"] in ids:
                raise SystemExit(f"{path}: model id {model['id']!r} is listed twice")
            ids.add(model["id"])
    return manifest


def selected(source):
    """The source's models, or the seeded draw its ``sample`` asks for."""
    models = sorted(source["models"], key=lambda m: m["id"])
    sample = source.get("sample")
    if sample and sample["count"] < len(models):
        rng = random.Random(sample["seed"])
        models = sorted(rng.sample(models, sample["count"]), key=lambda m: m["id"])
    return models


def classify(status, reason):
    """The class (roadmap 6.5.6) of one face's failure, from its status and
    the reason the mesher gives."""
    if status == "projection_failed" and " apart" in reason:
        return "gap_beyond_tolerance"
    if status == "unsupported" and "parameter range" in reason:
        return "trim_outside_patch"
    if status == "self_intersecting":
        return "self_intersecting_loop"
    if status == "open_loop":
        return "open_loop"
    if status in ("degenerate", "empty"):
        return "degenerate_patch"
    if status == "budget_exceeded":
        return "budget"
    return "other"


# --- fetch -------------------------------------------------------------------------------


def git(*args, cwd=None):
    subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True, text=True)


def fetch_git(source, into):
    """A sparse checkout of the models' paths at the pinned commit."""
    if into.exists():
        shutil.rmtree(into)
    into.mkdir(parents=True)
    git("init", "-q", cwd=into)
    git("remote", "add", "origin", source["url"], cwd=into)
    git(
        "sparse-checkout",
        "set",
        "--no-cone",
        *[m["path"] for m in source["models"]],
        cwd=into,
    )
    git(
        "fetch",
        "-q",
        "--depth",
        "1",
        "--filter=blob:none",
        "origin",
        source["commit"],
        cwd=into,
    )
    git("checkout", "-q", "FETCH_HEAD", cwd=into)


def fetch_https(source, into):
    into.mkdir(parents=True, exist_ok=True)
    name = source.get("file", source["url"].rsplit("/", 1)[-1])
    target = into / name
    with urllib.request.urlopen(source["url"], timeout=300) as r, open(
        target, "wb"
    ) as f:
        shutil.copyfileobj(r, f)
    if source.get("archive") == "zip":
        with zipfile.ZipFile(target) as z:
            for model in source["models"]:
                z.extract(model["path"], into)


def fetch(args):
    manifest = load_manifest(args.manifest)
    dest = pathlib.Path(args.dest)
    recorded, problems = {}, []
    for source in manifest["sources"]:
        into = dest / source["id"]
        try:
            (fetch_git if source["kind"] == "git" else fetch_https)(source, into)
        except (OSError, subprocess.CalledProcessError) as exc:
            detail = getattr(exc, "stderr", "") or str(exc)
            problems.append(f"{source['id']}: not fetched ({detail.strip()[:300]})")
            continue
        for model in source["models"]:
            path = into / model["path"]
            if not path.exists():
                problems.append(f"{model['id']}: {model['path']} is not in the source")
                continue
            digest = sha256(path)
            if not model.get("sha256"):
                recorded[model["id"]] = digest
            elif digest != model["sha256"]:
                problems.append(
                    f"{model['id']}: SHA-256 {digest} is not the manifest's"
                )
                path.unlink()
    for line in problems:
        print("cad_corpus: " + line, file=sys.stderr)
    if recorded:
        print(
            f"cad_corpus: {len(recorded)} model(s) had no hash; recorded",
            file=sys.stderr,
        )
        if args.record:
            pathlib.Path(args.record).write_text(canonical(recorded), encoding="utf-8")
    # A source that cannot be reached is reported, not fatal (a mismatched
    # hash is: the file is not the one measured before).
    return 1 if any("SHA-256" in p for p in problems) else 0


# --- run ---------------------------------------------------------------------------------


def measure_one(
    path,
    divisor,
    min_angle,
    max_points,
    virtual_edges=False,
    metric="mean",
    reparam=0.0,
    model_tolerances=False,
    virtual_faces=False,
    lop_corners=False,
):
    """Read, mesh and measure one model; runs in its own process."""
    import resource

    import numpy as np

    import meshioplusplus as mio

    out = {"status": "ok"}
    t0 = time.perf_counter()
    try:
        model = mio.read(path, file_format="step")
    except Exception as exc:  # noqa: BLE001 - every failure is a result here
        return {"status": "read_failed", "error": f"{type(exc).__name__}: {exc}"[:500]}
    out["read_s"] = round(time.perf_counter() - t0, 3)
    brep = model.info["brep"]
    lo, hi = model.points.min(axis=0), model.points.max(axis=0)
    diag = float(np.linalg.norm(hi - lo))
    out["faces"] = len(brep["faces"])
    out["edges"] = len(brep.get("edges", []))
    out["diagonal"] = diag
    out["target_size"] = diag / divisor
    t0 = time.perf_counter()
    try:
        mesh, report = mio.mesh_cad_surface(
            model,
            out["target_size"],
            min_angle=min_angle,
            max_points=max_points,
            virtual_edges=virtual_edges,
            metric=metric,
            reparam_condition_number=reparam,
            model_tolerances=model_tolerances,
            virtual_faces=virtual_faces,
            lop_corners=lop_corners,
            return_report=True,
        )
    except Exception as exc:  # noqa: BLE001
        out.update(status="mesh_failed", error=f"{type(exc).__name__}: {exc}"[:500])
        return out
    out["mesh_s"] = round(time.perf_counter() - t0, 3)
    statuses = collections.Counter(f["status"] for f in report["faces"])
    failures = [
        {
            "face": f["face"],
            "status": f["status"],
            "class": classify(f["status"], f["reason"]),
            "reason": f["reason"][:200],
        }
        for f in report["faces"]
        if f["status"] != "ok"
    ]
    out["faces_ok"] = statuses.get("ok", 0)
    out["failures_by_status"] = {k: v for k, v in sorted(statuses.items()) if k != "ok"}
    out["failures_by_class"] = dict(
        sorted(collections.Counter(f["class"] for f in failures).items())
    )
    out["failures"] = failures
    for key in (
        "num_points",
        "num_triangles",
        "num_open_edges",
        "num_non_manifold_edges",
        "num_misoriented_edges",
        "num_faces_below_bound",
        "max_edge_gap",
        "min_angle_deg",
        "max_angle_deg",
    ):
        out[key] = report[key]
    out["num_warnings"] = len(report["warnings"])
    if virtual_edges:
        out["num_virtual_edges"] = len(report["virtual_edges"])
    if reparam > 0.0:
        out["num_reparametrised"] = sum(r["taken"] for r in report["reparametrised"])
        out["num_reparam_declined"] = sum(
            not r["taken"] for r in report["reparametrised"]
        )
    if virtual_faces:
        out["num_merged_faces"] = sum(m["merged"] for m in report["merged_faces"])
        out["num_slivers"] = len(report["merged_faces"])
    if model_tolerances:
        out["num_tolerance_uses"] = len(report["tolerance_uses"])
    if lop_corners:
        out["num_lopped_corners"] = len(report["lopped_corners"])
        out["num_lop_guaranteed"] = sum(
            c["guaranteed"] for c in report["lopped_corners"]
        )
    tri = [c.data for c in mesh.cells if c.type == "triangle"]
    if tri:
        surface = mio.Mesh(mesh.points, [("triangle", np.concatenate(tri))])
        quality = mio.compute_quality(surface)["metrics"]
        out["quality"] = {
            name: {k: quality["quality:" + name][k] for k in ("min", "mean", "max")}
            for name in ("aspect_ratio", "min_angle", "max_angle")
        }
        # Against BRepMesh's triangulation of the same model: both lie on
        # the surfaces, so the distance is their two chordal errors at most.
        h = mio.hausdorff_distance(surface, model, face_samples=2)
        out["hausdorff_to_brepmesh"] = h["distance"] / diag if diag > 0 else 0.0
        a = mesh.points[np.concatenate(tri)]
        area = (
            0.5
            * np.linalg.norm(
                np.cross(a[:, 1] - a[:, 0], a[:, 2] - a[:, 0]), axis=1
            ).sum()
        )
        b = model.points[model.cells[0].data]
        ref = (
            0.5
            * np.linalg.norm(
                np.cross(b[:, 1] - b[:, 0], b[:, 2] - b[:, 0]), axis=1
            ).sum()
        )
        out["area_ratio"] = float(area / ref) if ref > 0 else None
    out["peak_rss_mb"] = round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1
    )
    return out


def run(args):
    import meshioplusplus as mio
    from meshioplusplus import _core

    manifest = load_manifest(args.manifest)
    dest = pathlib.Path(args.dest)
    only = set(args.models.split(",")) if args.models else None
    results = []
    for source in manifest["sources"]:
        for model in selected(source):
            if only and model["id"] not in only:
                continue
            path = dest / source["id"] / model["path"]
            row = {"id": model["id"], "source": source["id"]}
            if not path.exists():
                row["status"] = "unavailable"
                results.append(row)
                print(f"cad_corpus: {model['id']}: unavailable", file=sys.stderr)
                continue
            cmd = [
                sys.executable,
                str(pathlib.Path(__file__).resolve()),
                "_one",
                str(path),
                str(args.size_divisor),
                str(args.min_angle),
                str(args.max_points),
                str(int(args.virtual_edges)),
                args.metric,
                str(args.reparam),
                str(int(args.model_tolerances)),
                str(int(args.virtual_faces)),
                str(int(args.lop_corners)),
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
            results.append(row)
            print(
                f"cad_corpus: {model['id']}: {row['status']}, "
                f"{row.get('faces_ok', 0)}/{row.get('faces', '?')} faces, {row['wall_s']} s",
                file=sys.stderr,
                flush=True,
            )
    doc = {
        "Version": 1,
        "meshioplusplus": mio.__version__,
        "occt": bool(getattr(_core, "__has_occt__", False)),
        "settings": {
            "size_divisor": args.size_divisor,
            "min_angle": args.min_angle,
            "max_points": args.max_points,
            "timeout_s": args.timeout,
            **({"virtual_edges": True} if args.virtual_edges else {}),
            **({"metric": args.metric} if args.metric != "mean" else {}),
            **({"reparam_condition_number": args.reparam} if args.reparam > 0 else {}),
            **({"model_tolerances": True} if args.model_tolerances else {}),
            **({"virtual_faces": True} if args.virtual_faces else {}),
            **({"lop_corners": True} if args.lop_corners else {}),
        },
        "models": results,
    }
    pathlib.Path(args.out).write_text(canonical(doc), encoding="utf-8")
    return 0


def meshed(row):
    """A model meshes when it reads and every one of its faces meshes."""
    return row.get("status") == "ok" and row.get("faces_ok") == row.get("faces")


# --- baseline, compare -------------------------------------------------------------------


def stable(results):
    """The part of a run that does not depend on the machine: no times, and
    the measured lengths rounded to three significant digits."""

    def sig(x):
        return float(f"{x:.3g}") if isinstance(x, float) else x

    models = {}
    for row in results["models"]:
        entry = {"status": row["status"], "meshed": meshed(row)}
        for key in (
            "faces",
            "faces_ok",
            "failures_by_class",
            "num_open_edges",
            "num_faces_below_bound",
            "num_triangles",
        ):
            if key in row:
                entry[key] = row[key]
        for key in ("min_angle_deg", "hausdorff_to_brepmesh"):
            if key in row:
                entry[key] = sig(row[key])
        models[row["id"]] = entry
    return {
        "Version": 1,
        "meshioplusplus": results["meshioplusplus"],
        "settings": results["settings"],
        "models": models,
    }


def baseline(args):
    pathlib.Path(args.out).write_text(
        canonical(
            stable(json.loads(pathlib.Path(args.results).read_text(encoding="utf-8")))
        ),
        encoding="utf-8",
    )
    return 0


def compare(args):
    now = json.loads(pathlib.Path(args.results).read_text(encoding="utf-8"))
    base = json.loads(pathlib.Path(args.baseline).read_text(encoding="utf-8"))
    rows = {r["id"]: r for r in now["models"]}
    lost, notes = [], []
    for model_id, entry in sorted(base["models"].items()):
        row = rows.get(model_id)
        if row is None or row.get("status") == "unavailable":
            notes.append(f"{model_id}: not measured this run (unavailable)")
            continue
        if entry["meshed"] and not meshed(row):
            lost.append(
                f"{model_id}: meshed in the baseline, now {row['status']}"
                f" ({row.get('faces_ok', 0)}/{row.get('faces', '?')} faces)"
                + (f": {row['error']}" if row.get("error") else "")
            )
        elif not entry["meshed"] and meshed(row):
            notes.append(f"{model_id}: meshes now (it did not in the baseline)")
        elif entry.get("faces_ok", 0) > row.get("faces_ok", 0):
            notes.append(
                f"{model_id}: {entry['faces_ok']} faces meshed in the baseline, "
                f"{row.get('faces_ok', 0)} now"
            )
    for line in notes:
        print("cad_corpus: " + line)
    for line in lost:
        print("cad_corpus: REGRESSION: " + line)
    return 1 if lost else 0


# --- table, summary ----------------------------------------------------------------------


def render(base):
    models = base["models"]
    measured = [m for m in models.values() if m["status"] != "unavailable"]
    faces = sum(m.get("faces", 0) for m in measured)
    faces_ok = sum(m.get("faces_ok", 0) for m in measured)
    classes = collections.Counter()
    for m in measured:
        classes.update(m.get("failures_by_class", {}))
    s = base["settings"]
    lines = [
        BEGIN,
        "",
        f"Measured with meshio++ {base['meshioplusplus']}: `mesh_cad_surface` at a target size of "
        f"the model's diagonal over {s['size_divisor']:g}, refined to {s['min_angle']:g} degrees, "
        f"each model in its own process under {s['timeout_s']:g} s. "
        f"**{sum(m['meshed'] for m in measured)} of {len(measured)} models** mesh every face; "
        f"**{faces_ok} of {faces} faces** ({100.0 * faces_ok / max(faces, 1):.1f}%) mesh.",
        "",
        "| model | faces | meshed | failures by class | open edges | below bound | min angle | Hausdorff / diagonal |",
        "|---|---:|---:|---|---:|---:|---:|---:|",
    ]
    for model_id, m in sorted(models.items()):
        if m["status"] not in ("ok", "mesh_failed"):
            lines.append(
                f"| `{model_id}` | {m.get('faces', '–')} | {m['status']} | | | | | |"
            )
            continue
        fails = (
            ", ".join(f"{k} {v}" for k, v in m.get("failures_by_class", {}).items())
            or "–"
        )
        hd = m.get("hausdorff_to_brepmesh")
        ang = m.get("min_angle_deg")
        lines.append(
            f"| `{model_id}` | {m.get('faces', 0)} | {m.get('faces_ok', 0)} | {fails} | "
            f"{m.get('num_open_edges', '–')} | {m.get('num_faces_below_bound', '–')} | "
            f"{'–' if ang is None else f'{ang:g}'} | {'–' if hd is None else f'{hd:.2g}'} |"
        )
    lines += [
        "",
        "Failures by class: "
        + (", ".join(f"{k} {classes[k]}" for k in CLASSES if classes[k]) or "none")
        + ".",
        "",
        END,
    ]
    return "\n".join(lines)


def table(args):
    base = json.loads(pathlib.Path(args.baseline).read_text(encoding="utf-8"))
    doc = pathlib.Path(args.doc)
    text = doc.read_text(encoding="utf-8")
    if BEGIN not in text or END not in text:
        raise SystemExit(f"{doc}: the corpus markers are missing")
    head, rest = text.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    new = head + render(base) + tail
    if args.check:
        if new != text:
            print(
                f"cad_corpus: {doc} is stale; run tools/cad_corpus.py table",
                file=sys.stderr,
            )
            return 1
        return 0
    doc.write_text(new, encoding="utf-8")
    return 0


def summary(args):
    results = json.loads(pathlib.Path(args.results).read_text(encoding="utf-8"))
    print("### CAD corpus\n")
    print(render(stable(results)).replace(BEGIN, "").replace(END, "").strip())
    print("\n| model | status | read s | mesh s | peak MB | error |")
    print("|---|---|---:|---:|---:|---|")
    for row in results["models"]:
        print(
            f"| `{row['id']}` | {row['status']} | {row.get('read_s', '')} | "
            f"{row.get('mesh_s', '')} | {row.get('peak_rss_mb', '')} | {row.get('error', '')} |"
        )
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("fetch", help="download the sources and check their hashes")
    p.add_argument("--manifest", default=str(MANIFEST))
    p.add_argument("--dest", required=True)
    p.add_argument(
        "--record", help="write the hashes of models that had none to this JSON"
    )
    p = sub.add_parser("run", help="mesh and measure every fetched model")
    p.add_argument("--manifest", default=str(MANIFEST))
    p.add_argument("--dest", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--size-divisor", type=float, default=40.0)
    p.add_argument("--min-angle", type=float, default=25.0)
    p.add_argument("--max-points", type=int, default=5_000_000)
    p.add_argument("--timeout", type=float, default=900.0)
    p.add_argument("--models", help="a comma-separated subset of model ids")
    p.add_argument(
        "--virtual-edges",
        action="store_true",
        help="mesh with virtual edges across slivers of edges (6.5.2)",
    )
    p.add_argument(
        "--metric",
        choices=("mean", "varying"),
        default="mean",
        help="the metric each face is meshed in (6.4.1)",
    )
    p.add_argument(
        "--reparam",
        type=float,
        default=0.0,
        metavar="K",
        help="flatten faces whose condition number passes K (6.4.4); 0 is off",
    )
    p.add_argument(
        "--model-tolerances",
        action="store_true",
        help="accept the gaps the models declare (6.5.1)",
    )
    p.add_argument(
        "--virtual-faces",
        action="store_true",
        help="mesh sliver faces with a neighbour (6.5.2)",
    )
    p.add_argument(
        "--lop-corners",
        action="store_true",
        help="lop acute corners in place of concentric shells (6.2.3)",
    )
    p = sub.add_parser("compare", help="fail when a model that meshed no longer does")
    p.add_argument("results")
    p.add_argument("--baseline", default=str(BASELINE))
    p = sub.add_parser("baseline", help="write a run's stable part as the baseline")
    p.add_argument("results")
    p.add_argument("--out", default=str(BASELINE))
    p = sub.add_parser(
        "table", help="regenerate the corpus section of doc/benchmarks.md"
    )
    p.add_argument("--baseline", default=str(BASELINE))
    p.add_argument("--doc", default=str(DOC))
    p.add_argument("--check", action="store_true")
    p = sub.add_parser("summary", help="a Markdown summary of a run")
    p.add_argument("results")
    if argv is None:
        argv = sys.argv[1:]
    if argv and argv[0] == "_one":
        path, divisor, min_angle, max_points = argv[1:5]
        virtual = len(argv) > 5 and argv[5] == "1"
        metric = argv[6] if len(argv) > 6 else "mean"
        reparam = float(argv[7]) if len(argv) > 7 else 0.0
        tolerances = len(argv) > 8 and argv[8] == "1"
        virtual_faces = len(argv) > 9 and argv[9] == "1"
        lop_corners = len(argv) > 10 and argv[10] == "1"
        row = measure_one(
            path,
            float(divisor),
            float(min_angle),
            int(max_points),
            virtual,
            metric,
            reparam,
            tolerances,
            virtual_faces,
            lop_corners,
        )
        print(json.dumps(finite(row), allow_nan=False, default=float))
        return 0
    args = parser.parse_args(argv)
    return {
        "fetch": fetch,
        "run": run,
        "compare": compare,
        "baseline": baseline,
        "table": table,
        "summary": summary,
    }[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
