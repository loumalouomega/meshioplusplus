"""The ``parameterize`` verb: a UV atlas of a triangle surface
(``doc/parameterize.md``)."""

from .._helpers import read, write
from .._parameterize import parameterize
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(
        parser,
        "output mesh with uv (OBJ, glTF and USD write texture coordinates)",
        in_help="triangle surface mesh",
        output_format=False,
    )
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument(
        "--padding",
        type=float,
        default=0.02,
        help="gap between charts, a fraction of sqrt(surface area) (default 0.02)",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=8,
        help="how many times a piece that does not flatten may be halved (default 8)",
    )
    parser.add_argument(
        "--max-stretch",
        type=float,
        default=2.0,
        help="a chart stretched more is halved again (default 2; inf keeps it)",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def parameterize_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = parameterize(
        mesh,
        padding=args.padding,
        max_depth=args.max_depth,
        max_stretch=args.max_stretch,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return
    print(
        f"{report['num_charts']} chart(s), {report['num_flipped']} flipped "
        f"triangle(s), stretch {report['stretch']:.6g}, packing "
        f"{report['packing_efficiency']:.3g}, {report['num_seam_points']} seam point(s)"
    )
