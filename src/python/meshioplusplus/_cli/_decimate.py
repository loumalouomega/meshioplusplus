from .._decimate import decimate
from .._helpers import _writer_map, read, reader_map, write


def add_args(parser):
    parser.add_argument("infile", type=str, help="surface mesh file to be read from")
    parser.add_argument("outfile", type=str, help="decimated mesh to be written to")
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
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--ratio",
        type=float,
        default=None,
        help="fraction of the (triangulated) faces to KEEP, in (0, 1]",
    )
    group.add_argument(
        "--target-faces",
        type=int,
        default=None,
        help="absolute face count to stop at (within one collapse)",
    )
    group.add_argument(
        "--max-error",
        type=float,
        default=None,
        help="collapse only while the cheapest quadric error is at most this",
    )
    parser.add_argument(
        "--placement",
        choices=["optimal", "midpoint", "endpoint"],
        default="optimal",
        help=(
            "where the surviving vertex goes: the quadric minimizer (midpoint "
            "when ill-conditioned), the edge midpoint, or the cheaper endpoint "
            "(default: optimal)"
        ),
    )
    parser.add_argument(
        "--no-preserve-boundary",
        dest="preserve_boundary",
        action="store_false",
        help="allow boundary vertices to collapse (the outline may change)",
    )
    parser.add_argument(
        "--no-preserve-features",
        dest="preserve_features",
        action="store_false",
        help="allow feature vertices (sharp corners/creases) to collapse",
    )
    parser.add_argument(
        "--feature-angle",
        type=float,
        default=30.0,
        help="degrees between face normals above which a vertex is a feature "
        "(default: 30)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="do not print the collapse summary",
    )


def decimate_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = decimate(
        mesh,
        ratio=args.ratio,
        target_faces=args.target_faces,
        max_error=args.max_error,
        placement=args.placement,
        preserve_boundary=args.preserve_boundary,
        preserve_features=args.preserve_features,
        feature_angle=args.feature_angle,
        return_report=True,
    )
    if not args.quiet:
        faces_out = sum(len(cb.data) for cb in out.cells)
        print(
            f"decimated to {faces_out} faces "
            f"({report['faces_removed']} removed, "
            f"{report['points_removed']} points removed, "
            f"{report['collapses_rejected']} collapses rejected, "
            f"max error {report['max_error_applied']:.6g})"
        )
    write(args.outfile, out, file_format=args.output_format)
    return 0
