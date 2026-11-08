"""The ``metric`` verb: intersect, normalise, clamp and grade metric fields
(``doc/metric.md``)."""

from .._helpers import read, write
from .._metric import prepare_metric
from ._hull import _io_args
from ._json import emit_json


def add_args(parser):
    _io_args(parser)
    parser.add_argument(
        "--arrays",
        type=str,
        required=True,
        metavar="A,B",
        help="the metric fields (point data) to intersect",
    )
    parser.add_argument("--output", type=str, default="metric")
    parser.add_argument("--target-vertices", type=float, default=0.0)
    parser.add_argument("--norm", type=float, default=2.0)
    parser.add_argument("--hmin", type=float, default=0.0)
    parser.add_argument("--hmax", type=float, default=0.0)
    parser.add_argument("--gradation", type=float, default=0.0)
    parser.add_argument("--max-sweeps", type=int, default=100)


def metric_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = prepare_metric(
        mesh,
        [a for a in args.arrays.split(",") if a],
        output=args.output,
        target_vertices=args.target_vertices,
        norm=args.norm,
        hmin=args.hmin,
        hmax=args.hmax,
        gradation=args.gradation,
        max_sweeps=args.max_sweeps,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return
    print(
        f"metric complexity {report['complexity_before']:.6g} -> "
        f"{report['complexity_after']:.6g}; {report['num_graded']} point(s) graded"
    )
