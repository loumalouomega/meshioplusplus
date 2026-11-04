"""The surface-boolean spike on a corpus (roadmap 7.8, doc/spikes.md).

The spike's harness, ``meshioplusplus_boolean_spike`` (built with the C++
tests), runs union, intersection and difference on a pair of closed surfaces
and prints one JSON line per operation. This tool feeds it:

``synthetic``
    the generated corpus: the harness's own ``--synthetic SEED N`` (spheres,
    boxes and tori under seeded rigid motions, and the degenerate contacts),
    written as JSON lines.
``fetch``
    Thingi10K (Q. Zhou and A. Jacobson, 2016): draw a seeded sample of its
    closed, manifold, oriented models at the manifest's pinned revision and
    extract them, exactly as ``tools/tet_corpus.py fetch`` does (its code,
    this tool's manifest).
``run``
    pair each fetched model with a seeded rigid motion of itself (``self``)
    and with the next model in the sample, both scaled to the unit box
    (``cross``), and run the harness on each pair in its own process under a
    time limit, so a crash or a hang is a result naming the pair.
``summary``
    a Markdown table of a run's JSON lines, per class of pair: the share that
    succeeded, came out watertight and free of self-intersections, the
    failure reasons, the share with pieces that rounding turned over, the
    worst volume-identity error and the times.

Only the standard library, NumPy and meshioplusplus are used.
"""

from __future__ import annotations

import argparse
import collections
import importlib.util
import json
import math
import pathlib
import random
import re
import statistics
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
MANIFEST = HERE / "boolean_corpus" / "manifest.json"


def _tet_corpus():
    """tools/tet_corpus.py, whose fetch this tool shares."""
    spec = importlib.util.spec_from_file_location("tet_corpus", HERE / "tet_corpus.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def classify(case):
    """The class of a case name: ``seed78_12_sphere_box`` -> ``sphere_box``;
    ``thingi_self_123`` -> ``thingi_self``; a fixed case is its own class."""
    m = re.match(r"seed\d+_\d+_(.+)$", case)
    if m:
        return m.group(1)
    m = re.match(r"(thingi_(?:self|cross))_", case)
    if m:
        return m.group(1)
    return case


def rotation(rng):
    """A seeded rotation matrix (a uniform unit quaternion)."""
    while True:
        q = [rng.uniform(-1.0, 1.0) for _ in range(4)]
        n = math.sqrt(sum(c * c for c in q))
        if 1e-3 < n <= 1.0:
            break
    w, x, y, z = (c / n for c in q)
    return [
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ]


def pairs(ids, seed):
    """The pairs of a sample: for each model a ``self`` pair (the model and
    itself turned about its centre and moved by a fifth to a half of its
    size) and a ``cross`` pair (it and the next model, both in the unit box,
    the second turned and moved by up to a third); each with its motion."""
    rng = random.Random(seed)
    out = []
    for k, fid in enumerate(ids):
        d = rotation(rng)[0]
        out.append(
            {
                "case": f"thingi_self_{fid}",
                "a": fid,
                "b": fid,
                "unit": False,
                "rotation": rotation(rng),
                "shift": [c * rng.uniform(0.2, 0.5) for c in d],
            }
        )
        if len(ids) > 1:
            other = ids[(k + 1) % len(ids)]
            d = rotation(rng)[0]
            out.append(
                {
                    "case": f"thingi_cross_{fid}_{other}",
                    "a": fid,
                    "b": other,
                    "unit": True,
                    "rotation": rotation(rng),
                    "shift": [c * rng.uniform(0.0, 0.33) for c in d],
                }
            )
    return out


def place(points, rot, shift, unit):
    """The points scaled to the unit box (``unit``) or about their centre,
    turned by ``rot`` about the centre and moved by ``shift`` times the size."""
    import numpy as np

    p = np.asarray(points, dtype=float)
    lo, hi = p.min(axis=0), p.max(axis=0)
    centre, size = (lo + hi) / 2.0, float(np.max(hi - lo)) or 1.0
    q = (p - centre) / size if unit else p - centre
    q = q @ np.asarray(rot).T
    scale = 1.0 if unit else size
    return q + np.asarray(shift) * scale + (0.0 if unit else centre)


def synthetic(args):
    with open(args.out, "w", encoding="utf-8") as out:
        subprocess.run(
            [args.harness, "--synthetic", str(args.seed), str(args.count)],
            stdout=out,
            check=True,
        )
    return 0


def fetch(args):
    return _tet_corpus().fetch(args)


def run(args):
    import numpy as np

    import meshioplusplus as mio

    manifest = json.loads(pathlib.Path(args.manifest).read_text(encoding="utf-8"))
    dest = pathlib.Path(args.dest)
    ids = json.loads((dest / "sample.json").read_text(encoding="utf-8"))["ids"]
    ids = [i for i in ids if (dest / f"{i}.stl").exists()]
    if args.models:
        keep = {int(m) for m in args.models.split(",") if m}
        ids = [i for i in ids if i in keep]
    work = dest / "pairs"
    work.mkdir(exist_ok=True)
    crashed = 0
    with open(args.out, "w", encoding="utf-8") as out:
        for pair in pairs(ids, manifest["pairs"]["seed"]):
            a = mio.read(dest / f"{pair['a']}.stl")
            b = mio.read(dest / f"{pair['b']}.stl")
            tri_a = np.vstack([c.data for c in a.cells if c.type == "triangle"])
            tri_b = np.vstack([c.data for c in b.cells if c.type == "triangle"])
            ident = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            pa = place(a.points, ident, [0.0, 0.0, 0.0], pair["unit"])
            pb = place(b.points, pair["rotation"], pair["shift"], pair["unit"])
            fa, fb = work / f"{pair['case']}_a.vtu", work / f"{pair['case']}_b.vtu"
            mio.write(fa, mio.Mesh(pa, [("triangle", tri_a)]))
            mio.write(fb, mio.Mesh(pb, [("triangle", tri_b)]))
            try:
                res = subprocess.run(
                    [args.harness, str(fa), str(fb), "--case", pair["case"]],
                    capture_output=True,
                    text=True,
                    timeout=args.timeout,
                )
                lines = [ln for ln in res.stdout.splitlines() if ln.startswith("{")]
                if res.returncode != 0 or len(lines) != 3:
                    crashed += 1
                    reason = (
                        f"harness exit {res.returncode}: {res.stderr.strip()[-200:]}"
                    )
                    lines = [_failed(pair["case"], op, reason) for op in _OPS]
            except subprocess.TimeoutExpired:
                crashed += 1
                lines = [_failed(pair["case"], op, "timeout") for op in _OPS]
            for ln in lines:
                out.write(ln + "\n")
            out.flush()
            fa.unlink()
            fb.unlink()
    return 1 if crashed and args.strict else 0


_OPS = ("union", "intersection", "difference")


def _failed(case, op, why):
    return json.dumps({"case": case, "op": op, "ok": False, "failure": why})


def render(rows):
    """The Markdown summary of the JSON-line rows."""
    groups = collections.OrderedDict()
    for r in rows:
        groups.setdefault(classify(r["case"]), []).append(r)
    lines = [
        "| class | runs | ok | watertight | self-intersection free | "
        "rounding turned pieces | worst identity error | median / max s | failures |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]

    def share(rs, key):
        return f"{sum(1 for r in rs if r.get(key)) / len(rs):.0%}"

    for name, rs in groups.items():
        ok = [r for r in rs if r.get("ok")]
        errs = [
            max(
                r.get("identity_error_union", -1),
                r.get("identity_error_difference", -1),
            )
            for r in ok
        ]
        errs = [e for e in errs if e >= 0]
        times = [r["seconds"] for r in rs if "seconds" in r]
        reasons = collections.Counter(r["failure"] for r in rs if not r.get("ok"))
        turned = sum(1 for r in ok if r.get("rounding_inverted", 0) > 0)
        lines.append(
            f"| {name} | {len(rs)} | {share(rs, 'ok')} | {share(rs, 'watertight')} | "
            f"{share(rs, 'self_intersection_free')} | "
            f"{(turned / len(ok)) if ok else 0:.0%} | "
            f"{(max(errs) if errs else float('nan')):.1e} | "
            + (f"{statistics.median(times):.3f} / {max(times):.3f}" if times else "-")
            + " | "
            + (", ".join(f"{k} ({v})" for k, v in reasons.most_common()) or "-")
            + " |"
        )
    total = len(rows)
    lines.append("")
    lines.append(
        f"{total} runs: {sum(1 for r in rows if r.get('ok'))} succeeded, "
        f"{sum(1 for r in rows if r.get('ok') and r.get('watertight') and r.get('self_intersection_free'))}"
        " watertight and free of self-intersections."
    )
    return "\n".join(lines)


def summary(args):
    rows = [
        json.loads(ln)
        for ln in pathlib.Path(args.runs).read_text(encoding="utf-8").splitlines()
        if ln.strip()
    ]
    print(render(rows))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("synthetic", help="run the generated corpus")
    p.add_argument("--harness", required=True)
    p.add_argument("--seed", type=int, default=78)
    p.add_argument("--count", type=int, default=200)
    p.add_argument("--out", required=True)
    p.set_defaults(func=synthetic)
    p = sub.add_parser("fetch", help="draw the Thingi10K sample and extract it")
    p.add_argument("--manifest", default=str(MANIFEST))
    p.add_argument("--dest", required=True)
    p.add_argument("--base", default=None)
    p.add_argument("--archive", default=None)
    p.add_argument("--record", default=None)
    p.set_defaults(func=fetch)
    p = sub.add_parser("run", help="run the boolean on every pair of the sample")
    p.add_argument("--manifest", default=str(MANIFEST))
    p.add_argument("--dest", required=True)
    p.add_argument("--harness", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--timeout", type=float, default=300.0)
    p.add_argument("--models", default="")
    p.add_argument("--strict", action="store_true", help="exit 1 when a pair crashed")
    p.set_defaults(func=run)
    p = sub.add_parser("summary", help="a Markdown summary of a run")
    p.add_argument("runs")
    p.set_defaults(func=summary)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
