"""The ``point-normals`` verb: normals of a point cloud by local plane fits
(``doc/point_normals.md``)."""

from .._helpers import read, write
from .._point_normals import estimate_normals
from ._hull import _io_args
from ._json import emit_json


def add_args(parser):
    _io_args(parser)
    parser.add_argument("--k", type=int, default=16, help="neighbours per plane fit")
    parser.add_argument(
        "--orientation",
        choices=("consistent", "viewpoint", "none"),
        default="consistent",
    )
    parser.add_argument(
        "--viewpoint",
        type=str,
        default="0,0,0",
        metavar="X,Y,Z",
        help="the point normals face, with --orientation viewpoint",
    )
    parser.add_argument(
        "--variation", action="store_true", help="also attach normals:variation"
    )


def point_normals_cmd(args):
    try:
        viewpoint = [float(v) for v in args.viewpoint.split(",")]
    except ValueError:
        raise SystemExit("point-normals: --viewpoint expects three numbers") from None
    if len(viewpoint) != 3:
        raise SystemExit("point-normals: --viewpoint expects three numbers")
    mesh = read(args.infile, file_format=args.input_format)
    out, report = estimate_normals(
        mesh,
        k=args.k,
        orientation=args.orientation,
        viewpoint=viewpoint,
        record_variation=args.variation,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return 0
    text = f"normals of {len(mesh.points)} points"
    if report["num_components"]:
        text += f" in {report['num_components']} piece(s)"
    if report["num_degenerate"]:
        text += f", {report['num_degenerate']} degenerate (NaN)"
    print(text)
    return 0
