from .._deviation import deviation_regions
from .._helpers import read, write
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(
        parser,
        out_help="write the mesh with one cell region per deviation region here",
        in_help="measured surface carrying deviation:distance (from deviation)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        required=True,
        metavar="T",
        help="cells beyond +-T are excess or missing material",
    )
    parser.add_argument(
        "--field",
        type=str,
        default="deviation:distance",
        help="the signed deviation point data (default: deviation:distance)",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="deviation:",
        help="region names are <prefix>excess_<k> and <prefix>missing_<k>",
    )
    parser.add_argument(
        "--min-area",
        type=float,
        default=0.0,
        metavar="A",
        help="drop regions smaller than this area",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def deviation_regions_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, r = deviation_regions(
        mesh,
        args.tolerance,
        field=args.field,
        prefix=args.prefix,
        min_area=args.min_area,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(r)
    else:
        print(
            f"{r['num_excess']} excess and {r['num_missing']} missing region(s) beyond "
            f"+-{args.tolerance:g}"
        )
        for g in r["regions"]:
            c = g["centroid"]
            print(
                f"  {g['name']}: {g['num_cells']} cells, area {g['area']:.6g}, "
                f"extreme {g['extreme']:.6g}, volume {g['volume']:.6g}, "
                f"centroid ({c[0]:.6g}, {c[1]:.6g}, {c[2]:.6g})"
            )
        for w in r["warnings"]:
            print(f"  warning: {w}")
    return 0
