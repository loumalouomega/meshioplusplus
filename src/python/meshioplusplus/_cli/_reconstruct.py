"""The ``reconstruct`` verb: a triangle surface through an oriented point cloud
(``doc/reconstruct_surface.md``)."""

from .._helpers import read, write
from .._reconstruct_surface import reconstruct_surface
from ._hull import _io_args
from ._json import emit_json


def add_args(parser):
    _io_args(parser)
    parser.add_argument(
        "--normals", default="normals", help="the normals' point data name"
    )
    parser.add_argument(
        "--k", type=int, default=8, help="tangent planes blended per lattice node"
    )
    parser.add_argument(
        "--cell-size", type=float, default=0.0, help="lattice cell size (0: automatic)"
    )
    parser.add_argument(
        "--resolution",
        type=int,
        default=0,
        help="cells along the longest axis, with --cell-size 0",
    )
    parser.add_argument(
        "--max-distance", type=float, default=0.0, help="trim distance (0: automatic)"
    )
    parser.add_argument(
        "--no-estimate",
        action="store_true",
        help="fail instead of estimating missing normals",
    )
    parser.add_argument(
        "--normals-k", type=int, default=16, help="neighbours per estimated normal"
    )
    parser.add_argument(
        "--max-cells", type=int, default=20_000_000, help="the band's cell budget"
    )


def reconstruct_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = reconstruct_surface(
        mesh,
        normals=args.normals,
        k=args.k,
        cell_size=args.cell_size,
        resolution=args.resolution,
        max_distance=args.max_distance,
        estimate_normals=not args.no_estimate,
        normals_k=args.normals_k,
        max_cells=args.max_cells,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return
    text = (
        f"{report['num_triangles']} triangles in {report['num_components']} piece(s) "
        f"through {report['num_samples']} samples"
    )
    if report["estimated_normals"]:
        text += " (normals estimated)"
    print(text)
