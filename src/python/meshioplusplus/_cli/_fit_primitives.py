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
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(report)
    elif not args.quiet:
        print(f"fit_primitives: {report['num_patches']} patch(es)")
        for t in TYPES[1:]:
            print(f"  {t + ':':10s} {report['num_' + t]}")
        print(f"  {'unfitted:':10s} {report['num_unfitted']}")
    return 0
