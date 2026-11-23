"""The ``smooth-field`` verb: screened-Poisson smoothing of a point array on a
triangle surface (``doc/mesh_operator.md``)."""

from .._helpers import read, write
from .._smooth_field import smooth_field
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(
        parser,
        "output mesh with the smoothed array",
        in_help="triangle surface mesh",
        output_format=False,
    )
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument("--array", type=str, required=True, help="the point array")
    parser.add_argument(
        "--length",
        type=float,
        required=True,
        help="the smoothing length in the mesh's units (0 returns the field)",
    )
    parser.add_argument(
        "--output-name", type=str, default=None, help="default: <array>:smooth"
    )
    parser.add_argument("--tolerance", type=float, default=1e-12)
    parser.add_argument("--max-iterations", type=int, default=0)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def smooth_field_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = smooth_field(
        mesh,
        args.array,
        args.length,
        output=args.output_name,
        tolerance=args.tolerance,
        max_iterations=args.max_iterations,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return 0
    print(
        f"smoothed '{args.array}' on the scale {args.length} "
        f"({report['iterations']} iteration(s))"
    )
    return 0
