from .._decimate_volume import decimate_volume
from .._helpers import _writer_map, read, reader_map, write


def add_args(parser):
    parser.add_argument("infile", type=str, help="tet mesh file to be read from")
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
        help="fraction of the tets to KEEP, in (0, 1]",
    )
    group.add_argument(
        "--target-cells",
        type=int,
        default=None,
        help="absolute tet count to stop at (within one collapse)",
    )
    group.add_argument(
        "--max-error",
        type=float,
        default=None,
        help="collapse only while the cheapest boundary-touching quadric "
        "error is at most this",
    )
    parser.add_argument(
        "--placement",
        choices=["optimal", "midpoint", "endpoint"],
        default="optimal",
        help=(
            "where the surviving vertex goes: the quadric minimizer (midpoint "
            "when ill-conditioned or purely interior), the edge midpoint, or "
            "the cheaper endpoint (default: optimal)"
        ),
    )
    parser.add_argument(
        "--preserve-boundary",
        action="store_true",
        help="pin every boundary vertex outright, reproducing `decimate`'s own "
        "default instead of letting boundary vertices participate",
    )
    parser.add_argument(
        "--no-preserve-features",
        dest="preserve_features",
        action="store_false",
        help="allow boundary feature vertices (sharp corners/creases) to collapse",
    )
    parser.add_argument(
        "--feature-angle",
        type=float,
        default=30.0,
        help="degrees between boundary-triangle normals above which a vertex "
        "is a feature (default: 30)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="do not print the collapse summary",
    )


def decimate_volume_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = decimate_volume(
        mesh,
        ratio=args.ratio,
        target_cells=args.target_cells,
        max_error=args.max_error,
        placement=args.placement,
        preserve_boundary=args.preserve_boundary,
        preserve_features=args.preserve_features,
        feature_angle=args.feature_angle,
        return_report=True,
    )
    if not args.quiet:
        tets_out = sum(len(cb.data) for cb in out.cells)
        print(
            f"decimated to {tets_out} tets "
            f"({report['tets_removed']} removed, "
            f"{report['points_removed']} points removed, "
            f"{report['collapses_rejected']} collapses rejected, "
            f"max error {report['max_error_applied']:.6g})"
        )
    write(args.outfile, out, file_format=args.output_format)
    return 0
