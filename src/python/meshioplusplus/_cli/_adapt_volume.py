"""The ``adapt-volume`` verb: size-field tetrahedral remeshing
(``doc/adapt_volume.md``)."""

from .._adapt_volume import adapt_volume
from .._helpers import read, write
from ._hull import _io_args
from ._json import emit_json


def add_args(parser):
    _io_args(parser)
    parser.add_argument(
        "--hmin", type=float, default=0.0, help="smallest edge (<= 0: automatic)"
    )
    parser.add_argument(
        "--hmax", type=float, default=0.0, help="largest edge (<= 0: automatic)"
    )
    parser.add_argument(
        "--hgrad", type=float, default=1.3, help="gradation (default: 1.3)"
    )
    parser.add_argument(
        "--hausd",
        type=float,
        default=0.0,
        help="how far the boundary may move (<= 0: 1e-3 of the diagonal)",
    )
    parser.add_argument(
        "--hausd-relative",
        action="store_true",
        help="read --hausd as a fraction of the bounding-box diagonal",
    )
    parser.add_argument("--size-array", type=str, default=None)
    parser.add_argument("--curvature-sizing", action="store_true")
    parser.add_argument("--feature-angle", type=float, default=30.0)
    parser.add_argument("--corner-angle", type=float, default=30.0)
    parser.add_argument("--max-normal-deviation", type=float, default=30.0)
    parser.add_argument(
        "--preserve-boundary",
        action="store_true",
        help="leave the boundary vertices and triangles bit for bit",
    )
    parser.add_argument("--no-preserve-features", action="store_true")
    parser.add_argument(
        "--quality-metric",
        choices=("min_dihedral", "mean_ratio", "amips"),
        default="min_dihedral",
    )
    parser.add_argument("--quality-min", type=float, default=-1.0)
    parser.add_argument("--max-iterations", type=int, default=10)
    parser.add_argument("--target-min-dihedral", type=float, default=0.0)
    parser.add_argument(
        "--boundary-sizing",
        choices=("input", "geometry"),
        default="input",
        help="input: keep the boundary's resolution; geometry: size it by curvature",
    )
    for name in ("collapse", "split", "swap", "smooth"):
        parser.add_argument(f"--no-{name}", action="store_true")


def adapt_volume_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = adapt_volume(
        mesh,
        hmin=args.hmin,
        hmax=args.hmax,
        hgrad=args.hgrad,
        hausd=args.hausd,
        hausd_relative=args.hausd_relative,
        size_array=args.size_array,
        curvature_sizing=args.curvature_sizing,
        feature_angle=args.feature_angle,
        corner_angle=args.corner_angle,
        max_normal_deviation=args.max_normal_deviation,
        preserve_boundary=args.preserve_boundary,
        preserve_features=not args.no_preserve_features,
        quality_metric=args.quality_metric,
        quality_min=args.quality_min,
        max_iterations=args.max_iterations,
        target_min_dihedral=args.target_min_dihedral,
        collapse=not args.no_collapse,
        split=not args.no_split,
        swap=not args.no_swap,
        smooth=not args.no_smooth,
        boundary_sizing=args.boundary_sizing,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        # Timings differ from run to run; the JSON holds what is reproducible.
        emit_json(
            {
                "output": args.outfile,
                **{k: v for k, v in report.items() if k != "timings"},
            }
        )
        return 0
    print(
        f"tets {report['num_tets_before']} -> {report['num_tets_after']}, "
        f"points {report['num_points_before']} -> {report['num_points_after']}"
    )
    print(
        f"collapsed {report['num_collapsed']}, split {report['num_split']}, "
        f"swapped {report['num_swapped']}, smoothed {report['num_smoothed']} "
        f"in {report['iterations']} iteration(s)"
    )
    print(
        f"dihedral [{report['min_dihedral_before']:.2f}, "
        f"{report['max_dihedral_before']:.2f}] -> "
        f"[{report['min_dihedral_after']:.2f}, {report['max_dihedral_after']:.2f}] "
        f"degrees; boundary within {report['max_boundary_deviation']:.4g} "
        f"(hausd {report['hausd']:.4g})"
    )
    return 0
