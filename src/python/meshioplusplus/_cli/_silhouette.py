"""The ``silhouette`` verb: the outline of a part along a direction
(``doc/silhouette.md``)."""

from .._helpers import read, write
from .._silhouette import silhouette
from ._json import emit_json
from ._sweep_common import vec3


def add_args(parser):
    parser.add_argument("infile", type=str, help="a surface or volume mesh")
    parser.add_argument(
        "outfile",
        type=str,
        nargs="?",
        default=None,
        help="write the outline (or, with --extrude, the prism) here",
    )
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument(
        "--normal",
        type=str,
        default="0,0,1",
        help="the viewing direction 'x,y,z' (default: 0,0,1)",
    )
    parser.add_argument(
        "--clearance",
        type=float,
        default=0.0,
        help="grow the outline by this distance (negative shrinks it)",
    )
    parser.add_argument(
        "--round-angle",
        type=float,
        default=180.0,
        metavar="DEG",
        help="round the corners turning by less than DEG (default: 180, all)",
    )
    parser.add_argument("--arc-segments", type=int, default=64)
    parser.add_argument(
        "--extrude",
        action="store_true",
        help="the closed surface of the prism over the outline instead",
    )
    parser.add_argument("--region", type=str, default=None)
    parser.add_argument("--merge-tolerance", type=float, default=1e-9)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def silhouette_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, r = silhouette(
        mesh,
        vec3(args.normal, "silhouette", "normal"),
        clearance=args.clearance,
        round_angle=args.round_angle,
        arc_segments=args.arc_segments,
        extrude=args.extrude,
        region=args.region,
        merge_tolerance=args.merge_tolerance,
        return_report=True,
    )
    if args.outfile is not None:
        write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(r)
        return 0
    print(
        f"{r['num_loops']} loop(s): {r['num_outer']} outer, {r['num_holes']} hole(s); "
        f"area {r['area']:.9g}, perimeter {r['perimeter']:.9g}; heights "
        f"{r['min_height']:.9g} to {r['max_height']:.9g}"
    )
    for w in r["warnings"]:
        print(f"warning: {w}")
    return 0
