"""The ``filter-points`` verb: outlier removal and thinning of a point cloud
(``doc/filter_points.md``)."""

from .._filter_points import filter_points
from .._helpers import read, write
from ._hull import _io_args
from ._json import emit_json


def add_args(parser):
    _io_args(parser)
    parser.add_argument(
        "--outliers",
        choices=("statistical",),
        default=None,
        help="statistical outlier removal (uses --k and --std-ratio)",
    )
    parser.add_argument("--k", type=int, default=16, help="neighbours of --outliers")
    parser.add_argument(
        "--std-ratio",
        type=float,
        default=1.0,
        help="standard deviations above the mean for --outliers (default 1)",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=0.0,
        help="radius outlier removal: points need --min-neighbors within this radius",
    )
    parser.add_argument(
        "--min-neighbors", type=int, default=4, help="neighbours needed within --radius"
    )
    parser.add_argument(
        "--voxel",
        type=float,
        default=0.0,
        help="voxel thinning: one point per cube of this side",
    )
    parser.add_argument(
        "--voxel-mode",
        choices=("nearest", "average"),
        default="nearest",
        help="voxel representative: 'nearest' real point (default) or cube 'average'",
    )
    parser.add_argument(
        "--poisson",
        type=float,
        default=0.0,
        help="Poisson-disk thinning: kept points are at least this far apart",
    )


def filter_points_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    try:
        out, report = filter_points(
            mesh,
            statistical_k=args.k if args.outliers == "statistical" else 0,
            std_ratio=args.std_ratio,
            radius=args.radius,
            min_neighbors=args.min_neighbors,
            voxel_size=args.voxel,
            min_distance=args.poisson,
            voxel_mode=args.voxel_mode,
            return_report=True,
        )
    except ValueError as exc:
        raise SystemExit(f"filter-points: {exc}") from None
    write(args.outfile, out, file_format=args.output_format)
    kept = len(out.points)
    if args.json:
        report.pop("point_map")
        emit_json({"output": args.outfile, "num_kept": kept, **report})
        return 0
    print(
        f"kept {kept} of {report['num_input']} points "
        f"(statistical {report['num_removed_statistical']}, radius "
        f"{report['num_removed_radius']}, voxel {report['num_removed_voxel']}, "
        f"poisson {report['num_removed_poisson']} removed)"
    )
    for w in report["warnings"]:
        print(f"  warning: {w}")
    return 0
