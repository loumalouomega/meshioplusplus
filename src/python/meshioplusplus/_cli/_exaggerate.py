from .._deviation import exaggerate
from .._helpers import read, write
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(parser, out_help="write the exaggerated mesh here")
    parser.add_argument(
        "--field",
        type=str,
        default="deviation:distance",
        help="the signed deviation point data (default: deviation:distance)",
    )
    parser.add_argument(
        "--scale", type=float, default=10.0, help="displacement per unit deviation (10)"
    )
    parser.add_argument(
        "--clamp",
        type=float,
        default=0.0,
        help="clamp the deviation to +-C first (0: the suggested colour range)",
    )
    parser.add_argument(
        "--smooth-fraction",
        type=float,
        default=0.01,
        help="smooth over this fraction of the diagonal (0.01; 0 off)",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def exaggerate_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, r = exaggerate(
        mesh,
        field=args.field,
        scale=args.scale,
        clamp=args.clamp,
        smooth_fraction=args.smooth_fraction,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(r)
    else:
        print(
            f"{r['num_moved']} points moved; clamp +-{r['clamp']:.6g}, "
            f"smoothing radius {r['radius']:.6g}"
        )
        for w in r["warnings"]:
            print(f"  warning: {w}")
    return 0
