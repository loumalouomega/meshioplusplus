"""The ``subdivide-surface`` verb: Loop or Catmull-Clark subdivision with
Hoppe crease and boundary rules (``doc/subdivide_surface.md``)."""

from .._helpers import read, write
from .._subdivide_surface import subdivide_surface
from ._hull import _io_args
from ._json import emit_json


def add_args(parser):
    _io_args(parser)
    parser.add_argument(
        "--scheme",
        choices=("loop", "catmull_clark"),
        default="loop",
        help="'loop' (triangles, default) or 'catmull_clark' (polygons to quads)",
    )
    parser.add_argument("--levels", type=int, default=1, help="subdivision levels")
    parser.add_argument(
        "--crease-angle",
        type=float,
        default=60.0,
        help="dihedral angle in degrees above which an edge is sharp (0 disables)",
    )
    parser.add_argument(
        "--no-preserve-boundary",
        action="store_true",
        help="smooth open boundaries like interior edges",
    )
    parser.add_argument(
        "--no-interpolate",
        action="store_true",
        help="do not interpolate floating-point point data",
    )
    parser.add_argument(
        "--record-parent-ids",
        action="store_true",
        help="attach subdivide:parent_cell",
    )


def subdivide_surface_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = subdivide_surface(
        mesh,
        args.scheme,
        levels=args.levels,
        crease_angle=args.crease_angle,
        preserve_boundary=not args.no_preserve_boundary,
        interpolate_point_data=not args.no_interpolate,
        record_parent_ids=args.record_parent_ids,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(
            {
                "output": args.outfile,
                "num_cells": int(sum(len(b.data) for b in out.cells)),
                "num_crease_edges": int(report["num_crease_edges"]),
                "num_boundary_edges": int(report["num_boundary_edges"]),
            }
        )
        return 0
    for w in report.get("warnings", []):
        print(f"warning: {w}")
    print(
        f"subdivided: {report['num_crease_edges']} crease edge(s), "
        f"{report['num_boundary_edges']} boundary edge(s)"
    )
    return 0
