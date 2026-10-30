from .._helpers import _writer_map, read, reader_map, write
from .._thickness import thickness
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="surface mesh file to be measured")
    parser.add_argument(
        "outfile",
        type=str,
        nargs="?",
        default=None,
        help="optional: write the mesh with thickness:distance and thickness:hit_cell here",
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
        "--output-format",
        "-o",
        type=str,
        choices=sorted(list(_writer_map.keys())),
        help="output file format",
        default=None,
    )
    parser.add_argument(
        "--other",
        type=str,
        default=None,
        metavar="FILE",
        help="a second surface: measure the clearance to it",
    )
    parser.add_argument(
        "--other-format",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input format of --other",
        default=None,
    )
    parser.add_argument(
        "--direction",
        type=str,
        choices=("auto", "inward", "outward"),
        default="auto",
        help="inward (wall thickness), outward (clearance) or auto: inward for one "
        "surface, outward for two (default)",
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=0.0,
        metavar="D",
        help="longest ray (default: 0, the bounding-box diagonal)",
    )
    parser.add_argument(
        "--min",
        type=float,
        default=None,
        metavar="T",
        help="exit with status 1 when any distance is below T",
    )
    parser.add_argument(
        "--region", type=str, default="", help="named cell region to measure"
    )
    parser.add_argument(
        "--region-other",
        type=str,
        default="",
        help="named cell region of the target surface (--other, else the input)",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def thickness_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    other = None
    if args.other is not None:
        other = read(args.other, file_format=args.other_format)
    out, r = thickness(
        mesh,
        other,
        direction=args.direction,
        max_distance=args.max_distance,
        threshold=args.min,
        region=args.region,
        region_other=args.region_other,
        attach=args.outfile is not None,
        return_report=True,
    )
    if args.outfile is not None:
        write(args.outfile, out, file_format=args.output_format)
    failed = args.min is not None and r["num_below"] > 0
    if args.json:
        obj = dict(r)
        if args.min is not None:
            obj["threshold"] = args.min
            obj["passed"] = not failed
        emit_json(obj)
    else:
        print(
            f"{r['direction']} rays from {r['num_cells']} cells: {r['num_hit']} hit, "
            f"{r['num_missed']} missed"
        )
        if r["num_hit"] > 0:
            print(
                f"  min {r['min']:.9g} (cell {r['min_cell']}), max {r['max']:.9g}, "
                f"mean {r['mean']:.9g}"
            )
        if args.min is not None:
            print(
                f"  {'FAIL' if failed else 'pass'}: {r['num_below']} cells below {args.min:g}"
            )
        for w in r["warnings"]:
            print(f"  warning: {w}")
    return 1 if failed else 0
