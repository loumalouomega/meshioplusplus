"""The ``locate`` verb: the cell containing each point of a file
(``doc/ray_queries.md``)."""

from .._cast_rays import _located_mesh, _xyz, locate_points
from .._helpers import write
from ._common import add_transfer_io_args, read_source_target
from ._json import emit_json


def add_args(parser):
    add_transfer_io_args(
        parser,
        "vertex mesh of the located points (locate:cell, locate:nodes, locate:barycentric)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-12,
        help="a barycentric weight at least -tolerance counts as inside (default 1e-12)",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def locate_cmd(args):
    mesh, points = read_source_target(args)
    xyz = _xyz(points.points)
    res = locate_points(mesh, xyz, tolerance=args.tolerance)
    write(args.outfile, _located_mesh(xyz, res), file_format=args.output_format)
    if args.json:
        emit_json(
            {
                "output": args.outfile,
                "num_points": len(xyz),
                "num_located": res["num_located"],
            }
        )
        return
    print(f"{res['num_located']} of {len(xyz)} points lie inside the mesh")
