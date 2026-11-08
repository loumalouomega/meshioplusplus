"""The ``untangle`` verb: move nodes until no tetrahedron or hexahedron
corner is inverted (``doc/untangle.md``)."""

from .._helpers import read, write
from .._untangle import untangle
from ._hull import _io_args
from ._json import emit_json


def add_args(parser):
    _io_args(parser)
    parser.add_argument(
        "--no-fix-boundary", action="store_true", help="let the boundary nodes move"
    )
    parser.add_argument(
        "--rings",
        type=int,
        default=2,
        help="edge rings around the inverted cells whose nodes may move",
    )
    parser.add_argument(
        "--smooth", action="store_true", help="untangle and smooth the whole mesh"
    )
    parser.add_argument("--max-sweeps", type=int, default=100)


def untangle_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = untangle(
        mesh,
        fix_boundary=not args.no_fix_boundary,
        rings=args.rings,
        smooth=args.smooth,
        max_sweeps=args.max_sweeps,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return
    print(
        f"inverted cells {report['num_inverted_before']} -> "
        f"{report['num_inverted_after']}; least corner scaled Jacobian "
        f"{report['min_quality_before']:.4g} -> {report['min_quality_after']:.4g} "
        f"({report['sweeps']} sweep(s))"
    )
