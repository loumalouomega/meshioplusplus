from .._fit_primitives import TYPES, fit_primitives
from .._helpers import read, write
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(
        parser,
        "mesh annotated with the fits to be written to",
        in_help="segmented surface mesh to be read from",
    )
    parser.add_argument(
        "--patch-data",
        type=str,
        default="segment:patch",
        help="integer cell-data patch label; '' fits the whole surface "
        "(default: segment:patch)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-3,
        help="acceptance tolerance on the orthogonal distance (default: 1e-3)",
    )
    parser.add_argument(
        "--absolute",
        action="store_true",
        help="the tolerance is a length, not a fraction of each patch's diagonal",
    )
    parser.add_argument(
        "--min-inlier-fraction",
        type=float,
        default=0.9,
        help="share of points that must lie within tolerance (default: 0.9)",
    )
    parser.add_argument(
        "--types",
        type=str,
        default=",".join(TYPES[1:]),
        help="comma-separated primitives to try (default: all)",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=50,
        help="Levenberg-Marquardt iteration cap per fit (default: 50)",
    )
    parser.add_argument(
        "--region", type=str, default="", help="fit this named cell region as one patch"
    )
    parser.add_argument(
        "--seed",
        type=str,
        action="append",
        default=None,
        metavar="X,Y,Z",
        help="grow and fit a feature from the facet nearest this point (repeatable)",
    )
    parser.add_argument(
        "--feature-angle",
        type=float,
        default=30.0,
        help="seed growth: edges sharper than this are not crossed (default: 30)",
    )
    parser.add_argument(
        "--rim-angle",
        type=float,
        default=15.0,
        help="seed growth: peel rim faces turning more from the fit (default: 15)",
    )
    parser.add_argument(
        "--output-region",
        type=str,
        default="fit:feature",
        help="seed growth: the cell region written (default: fit:feature)",
    )
    parser.add_argument(
        "--used-points",
        type=str,
        choices=("auto", "all", "3sigma", "2sigma", "1sigma"),
        default="auto",
        help="sigma trimming (default: auto, 3sigma with --region or --seed)",
    )
    parser.add_argument("--json", action="store_true", help="emit the counts as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def fit_primitives_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = fit_primitives(
        mesh,
        patch_data=args.patch_data,
        tolerance=args.tolerance,
        relative=not args.absolute,
        min_inlier_fraction=args.min_inlier_fraction,
        types=[t.strip() for t in args.types.split(",") if t.strip()],
        max_iterations=args.max_iterations,
        return_report=True,
        region=args.region or None,
        seeds=(
            [[float(v) for v in s.split(",")] for s in args.seed] if args.seed else None
        ),
        feature_angle=args.feature_angle,
        rim_angle=args.rim_angle,
        output_region=args.output_region,
        used_points=args.used_points,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(report)
    elif not args.quiet:
        print(f"fit_primitives: {report['num_patches']} patch(es)")
        for t in TYPES[1:]:
            print(f"  {t + ':':10s} {report['num_' + t]}")
        print(f"  {'unfitted:':10s} {report['num_unfitted']}")
        if report["num_feature_cells"]:
            print(f"  feature:   {report['num_feature_cells']} cells")
        for w in report["warnings"]:
            print(f"  warning: {w}")
    return 0
