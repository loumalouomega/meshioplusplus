import json

from .._helpers import read, write
from .._mirror import mirror
from ._sweep_common import add_io_args, add_report_args, vec3


def add_args(parser):
    add_io_args(parser, "mirrored mesh to be written to")
    parser.add_argument(
        "--origin",
        type=str,
        default="0,0,0",
        help="a point of the plane 'x,y,z' (negatives need --origin=-1,0,0)",
    )
    parser.add_argument(
        "--normal",
        type=str,
        default="1,0,0",
        help="the plane's normal 'x,y,z' (negatives need --normal=0,0,-1)",
    )
    parser.add_argument(
        "--keep-original",
        action="store_true",
        help="write the input followed by its image",
    )
    parser.add_argument(
        "--no-weld",
        action="store_true",
        help="with --keep-original, do not share the points on the plane",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-9,
        help="on the plane: fraction of the bounding-box diagonal (default 1e-9)",
    )
    parser.add_argument(
        "--reflect-data",
        action="store_true",
        help="reflect float vector and tensor data on the reflected copies",
    )
    add_report_args(parser)


def mirror_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = mirror(
        mesh,
        vec3(args.origin, "mirror", "origin"),
        vec3(args.normal, "mirror", "normal"),
        keep_original=args.keep_original,
        weld=not args.no_weld,
        tolerance=args.tolerance,
        reflect_vector_data=args.reflect_data,
        return_report=True,
    )
    if args.json:
        print(json.dumps(report))
    elif not args.quiet:
        print(
            f"{report['num_flipped_cells']} cells reflected, "
            f"{report['num_points_on_plane']} points on the plane, "
            f"{report['num_dropped_cells']} cells and {report['num_dropped_sides']} "
            "side entries dropped"
        )
        for w in report["warnings"]:
            print(f"warning: {w}")
    write(args.outfile, out, file_format=args.output_format)
