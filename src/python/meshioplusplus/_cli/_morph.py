"""The ``morph`` verb: move the boundary (or chosen points) by a
displacement, harmonically extended into the mesh (``doc/morph.md``)."""

from .._helpers import read, write
from .._morph import morph
from ._hull import _io_args
from ._json import emit_json


def add_args(parser):
    _io_args(parser)
    parser.add_argument(
        "--array", type=str, required=True, help="point data: the displacement"
    )
    parser.add_argument(
        "--no-fix-boundary",
        action="store_true",
        help="do not constrain the boundary points",
    )
    parser.add_argument(
        "--constrained-array",
        type=str,
        default="",
        help="point data whose nonzero entries constrain their points too",
    )
    parser.add_argument("--stiffening", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument(
        "--allow-inverted",
        action="store_true",
        help="write a result with inverted cells instead of failing",
    )
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--max-iterations", type=int, default=0)


def morph_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = morph(
        mesh,
        args.array,
        fix_boundary=not args.no_fix_boundary,
        constrained=args.constrained_array or None,
        stiffening=args.stiffening,
        steps=args.steps,
        allow_inverted=args.allow_inverted,
        tolerance=args.tolerance,
        max_iterations=args.max_iterations,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return
    text = (
        f"morphed: {report['num_constrained']} constrained, {report['num_free']} free "
        f"point(s); largest displacement {report['max_displacement']:.9g}"
    )
    if report["num_inverted"]:
        text += f"; {report['num_inverted']} inverted cell(s)"
    print(text)
