from .._helpers import _writer_map, read, reader_map, write
from .._smooth import smooth


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="smoothed mesh to be written to")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--output-format",
        "-o",
        type=str,
        choices=sorted(list(_writer_map.keys())),
        help="output file format",
        default=None,
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=["taubin", "laplacian", "odt"],
        default="taubin",
        help=(
            "smoothing operator: taubin (shrink-free, default), "
            "laplacian (stronger per pass, shrinks volume), or "
            "odt (optimal-Delaunay-triangulation smoothing -- tet-only, "
            "C++-core only with no pure-Python fallback)"
        ),
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=10,
        help=(
            "number of iterations; for taubin one iteration is two passes "
            "(+lambda then mu) (default: 10)"
        ),
    )
    parser.add_argument(
        "--lambda",
        dest="lambda_",
        type=float,
        default=-1.0,
        help=(
            "relaxation factor of the smoothing pass, in (0, 1); "
            "the default depends on --method (0.5 for laplacian, "
            "0.33 for taubin, 0.9 for odt)"
        ),
    )
    parser.add_argument(
        "--mu",
        type=float,
        default=-0.34,
        help=(
            "taubin only: un-shrinking factor, must satisfy mu < -lambda < 0; "
            "negative values need the --mu=-0.34 form (default: -0.34)"
        ),
    )
    parser.add_argument(
        "--no-fix-boundary",
        dest="fix_boundary",
        action="store_false",
        help="let boundary nodes move too (they are pinned by default)",
    )
    parser.add_argument(
        "--no-preserve-features",
        dest="preserve_features",
        action="store_false",
        help="do not pin feature nodes (sharp corners/creases are kept by default)",
    )
    parser.add_argument(
        "--feature-angle",
        type=float,
        default=30.0,
        help=(
            "angle in degrees between boundary facet normals above which their "
            "shared nodes are pinned as features (default: 30)"
        ),
    )
    parser.add_argument(
        "--no-guard-inversion",
        dest="guard_inversion",
        action="store_false",
        help="do not reject moves that would invert an incident cell",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary output"
    )


def smooth_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    out, report = smooth(
        mesh,
        method=args.method,
        iterations=args.iterations,
        lambda_=args.lambda_,
        mu=args.mu,
        fix_boundary=args.fix_boundary,
        preserve_features=args.preserve_features,
        feature_angle=args.feature_angle,
        guard_inversion=args.guard_inversion,
        return_report=True,
    )

    if not args.quiet:
        print("smoothed mesh")
        print(f"  nodes moved:              {report['num_nodes_moved']}")
        print(f"  max displacement:         {report['max_displacement']:.6g}")
        print(f"  skipped (inversion):      {report['num_skipped_inversion']}")

    write(args.outfile, out, file_format=args.output_format)
    return 0
