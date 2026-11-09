from .._helpers import _reader_names, _writer_names, read, write
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
        choices=_reader_names(),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--output-format",
        "-o",
        type=str,
        choices=_writer_names(),
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
        choices=_reader_names(),
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
    parser.add_argument(
        "--method",
        choices=("ray", "sphere", "cone"),
        default="ray",
        help="one ray, the sphere at its midpoint, or a cone of rays (default: ray)",
    )
    parser.add_argument(
        "--cone-angle",
        type=float,
        default=15.0,
        metavar="DEG",
        help="the cone's half-angle (default: 15)",
    )
    parser.add_argument(
        "--cone-rays", type=int, default=16, help="the cone's rays (default: 16)"
    )
    parser.add_argument(
        "--two-sided",
        action="store_true",
        help="rays both ways, either face of a triangle counting (scan winding)",
    )
    parser.add_argument(
        "--max-normal-deviation",
        type=float,
        default=90.0,
        metavar="DEG",
        help="step over faces more edge-on to the ray than this (default: 90)",
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
        method=args.method,
        cone_angle=args.cone_angle,
        cone_rays=args.cone_rays,
        two_sided=args.two_sided,
        max_normal_deviation=args.max_normal_deviation,
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
            f"{r['direction']} {r['method']} from {r['num_cells']} cells: {r['num_hit']} hit, "
            f"{r['num_missed']} missed"
        )
        if r["num_hit"] > 0:
            print(
                f"  min {r['min']:.9g} (cell {r['min_cell']}), max {r['max']:.9g}, "
                f"mean {r['mean']:.9g}, p2 {r['p2']:.6g}, p95 {r['p95']:.6g}"
            )
        if args.min is not None:
            print(
                f"  {'FAIL' if failed else 'pass'}: {r['num_below']} cells below {args.min:g}"
            )
        for w in r["warnings"]:
            print(f"  warning: {w}")
    return 1 if failed else 0
