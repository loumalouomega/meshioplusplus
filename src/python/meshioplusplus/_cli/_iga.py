from .._helpers import _writer_map, read, reader_map, write
from .._iga import iga_evaluate, iga_refine, iga_sample
from ._json import emit_json


def _add_io(parser, what):
    parser.add_argument(
        "infile",
        type=str,
        help="NURBS model to be read from (MFEM NURBS, CAD JSON, STEP, IGES)",
    )
    if what:
        parser.add_argument("outfile", type=str, help=what)
        parser.add_argument(
            "--output-format",
            "-o",
            type=str,
            choices=sorted(list(_writer_map.keys())),
            help="output file format",
            default=None,
        )
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--model",
        choices=["auto", "nurbs", "brep"],
        default="auto",
        help="which model of the mesh to use when it carries both (default: auto)",
    )


def _ints(text, name):
    try:
        values = [int(v) for v in str(text).split(",") if v.strip()]
    except ValueError:
        raise SystemExit(f"{name}: expected integers, got {text!r}")
    return values[0] if len(values) == 1 else values


def _add_sampling(parser):
    parser.add_argument(
        "--trim-samples",
        type=int,
        default=16,
        help="B-rep: points per knot span a trimming curve is sampled at (default: 16)",
    )
    parser.add_argument(
        "--no-edges", action="store_true", help="B-rep: do not sample the edge curves"
    )
    parser.add_argument(
        "--no-trim",
        action="store_true",
        help="B-rep: keep every cell, ignoring the trimming loops",
    )


def add_refine_args(parser):
    _add_io(parser, "refined model, sampled, to be written to")
    parser.add_argument(
        "--subdivisions",
        "-s",
        type=str,
        default="1",
        help="parts every knot span is split into: N, or one per direction N,M,L (default: 1)",
    )
    parser.add_argument(
        "--elevate",
        "-e",
        type=str,
        default="0",
        help="degrees to add: N, or one per direction N,M,L (default: 0)",
    )
    parser.add_argument(
        "--sample-subdivisions",
        type=int,
        default=0,
        help="cells per knot span of the written sampling (default: 0, automatic)",
    )
    _add_sampling(parser)
    parser.add_argument("--json", action="store_true", help="emit a summary as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def refine_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out = iga_refine(
        mesh,
        subdivisions=_ints(args.subdivisions, "--subdivisions"),
        elevate=_ints(args.elevate, "--elevate"),
        sample_subdivisions=args.sample_subdivisions,
        model=args.model,
        trim_samples=args.trim_samples,
        edges=not args.no_edges,
        trim=not args.no_trim,
    )
    write(args.outfile, out, file_format=args.output_format)
    _summary(out, args, "iga_refine")
    return 0


def add_sample_args(parser):
    _add_io(parser, "sampled model to be written to")
    parser.add_argument(
        "--subdivisions",
        "-s",
        type=int,
        default=0,
        help="cells per knot span and direction (default: 0, automatic)",
    )
    _add_sampling(parser)
    parser.add_argument("--json", action="store_true", help="emit a summary as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def sample_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out = iga_sample(
        mesh,
        subdivisions=args.subdivisions,
        model=args.model,
        trim_samples=args.trim_samples,
        edges=not args.no_edges,
        trim=not args.no_trim,
    )
    write(args.outfile, out, file_format=args.output_format)
    _summary(out, args, "iga_sample")
    return 0


def _summary(mesh, args, name):
    cells = sum(len(c.data) for c in mesh.cells)
    kind = "nurbs" if "nurbs" in (mesh.info or {}) else "brep"
    patches = len((mesh.info or {}).get(kind, {}).get("patches", []))
    if args.json:
        emit_json(
            {
                "operation": name,
                "model": kind,
                "num_patches": patches,
                "num_points": len(mesh.points),
                "num_cells": cells,
                "cell_blocks": [
                    {"type": c.type, "num_cells": len(c.data)} for c in mesh.cells
                ],
            }
        )
    elif not args.quiet:
        print(f"{name}: {patches} patch(es), {len(mesh.points)} points, {cells} cells")


def add_evaluate_args(parser):
    _add_io(parser, None)
    parser.add_argument(
        "--param",
        "-p",
        action="append",
        required=True,
        help="parametric coordinates u[,v[,w]] to evaluate at; repeat for more points",
    )
    parser.add_argument(
        "--patch",
        type=int,
        default=0,
        help="patch of a NURBS model, face of a B-rep (default: 0)",
    )
    parser.add_argument(
        "--derivatives",
        type=int,
        choices=[0, 1],
        default=0,
        help="1: also the first derivatives and unit normals (default: 0)",
    )
    parser.add_argument(
        "--trim-samples",
        type=int,
        default=16,
        help="B-rep: points per knot span a trimming curve is sampled at (default: 16)",
    )
    # The verb prints JSON regardless -- emit_json in evaluate_cmd -- so the
    # flag changes nothing; the other verbs take it, and test_cli_json.py runs
    # every verb both CLIs share with --json and holds the exit codes equal.
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit a summary as JSON (always; accepted for uniformity)",
    )


def evaluate_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    params = [[float(v) for v in text.split(",")] for text in args.param]
    res = iga_evaluate(
        mesh,
        params,
        patch=args.patch,
        derivatives=args.derivatives,
        model=args.model,
        trim_samples=args.trim_samples,
    )
    # The native CLI's shape: arrays flattened row-major, widths beside them.
    out = {
        "num_points": len(params),
        "space_dim": res["points"].shape[1],
        "points": res["points"].ravel(),
        "fields": [
            {"components": f.shape[1], "values": f.ravel()} for f in res["fields"]
        ],
    }
    for key in ("derivatives", "normals", "inside"):
        if key in res:
            out[key] = res[key].ravel()
    emit_json(out)
    return 0
