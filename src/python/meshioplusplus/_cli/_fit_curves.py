"""The ``fit-curves`` verb: lines, arcs and circles from chained sections
(``doc/fit_curves.md``)."""

from .._fit_curves import fit_curves
from .._helpers import read, write
from ._json import emit_json
from ._sweep_common import vec3


def add_args(parser):
    parser.add_argument(
        "infile", type=str, help="a chained line mesh (slice --chain, sections)"
    )
    parser.add_argument(
        "outfile",
        type=str,
        nargs="?",
        default=None,
        help="write the annotated mesh here",
    )
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.0,
        help="the largest residual of a piece (0: 1e-3 of the diagonal)",
    )
    parser.add_argument(
        "--snap-angle",
        type=float,
        default=0.0,
        metavar="DEG",
        help="snap lines within DEG of the axes or diagonals (default: 0)",
    )
    parser.add_argument(
        "--normal",
        type=str,
        default=None,
        help="the sketch normal 'x,y,z' (default: from the points)",
    )
    parser.add_argument("--min-arc-points", type=int, default=4)
    parser.add_argument("--min-arc-sweep", type=float, default=5.0, metavar="DEG")
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def fit_curves_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, r = fit_curves(
        mesh,
        tolerance=args.tolerance,
        snap_angle=args.snap_angle,
        normal=(
            None if args.normal is None else vec3(args.normal, "fit-curves", "normal")
        ),
        min_arc_points=args.min_arc_points,
        min_arc_sweep=args.min_arc_sweep,
        return_report=True,
    )
    if args.outfile is not None:
        write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(r)
        return 0
    print(
        f"{r['num_lines']} line(s), {r['num_arcs']} arc(s), {r['num_circles']} circle(s); "
        f"max residual {r['max_residual']:.6g} (tolerance {r['tolerance']:.6g})"
    )
    for p in r["pieces"]:
        extra = (
            ""
            if p["type"] == "line"
            else f" centre ({p['centre'][0]:.9g}, {p['centre'][1]:.9g}) radius "
            f"{p['radius']:.9g} sweep {p['sweep']:.6g}"
        )
        print(
            f"  chain {p['chain']} {p['type']}: ({p['start'][0]:.9g}, {p['start'][1]:.9g}) -> "
            f"({p['end'][0]:.9g}, {p['end'][1]:.9g}){extra}"
        )
    for w in r["warnings"]:
        print(f"warning: {w}")
    return 0
