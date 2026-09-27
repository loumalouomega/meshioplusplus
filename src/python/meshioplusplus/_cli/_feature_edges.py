from .._feature_edges import feature_edges
from .._helpers import _writer_map, read, reader_map, write
from ._json import emit_json


def add_args(parser):
    parser.add_argument(
        "infile", type=str, help="surface or volume mesh to be read from"
    )
    parser.add_argument(
        "outfile", type=str, help="line mesh of the reported edges to be written to"
    )
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
        "--angle",
        type=float,
        default=30.0,
        metavar="DEG",
        help="largest dihedral angle still treated as smooth, 0-180 (default: 30)",
    )
    for flag, what in (
        ("feature", "sharp edges"),
        ("boundary", "open edges"),
        ("non-manifold", "edges used by three or more faces"),
        ("inconsistent", "face pairs wound the same way"),
    ):
        parser.add_argument(
            f"--no-{flag}", action="store_true", help=f"do not report {what}"
        )
    parser.add_argument(
        "--region",
        type=str,
        default="",
        help="restrict to this named cell region (default: every cell)",
    )
    parser.add_argument("--json", action="store_true", help="emit the counts as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def feature_edges_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = feature_edges(
        mesh,
        feature_angle=args.angle,
        feature=not args.no_feature,
        boundary=not args.no_boundary,
        non_manifold=not args.no_non_manifold,
        inconsistent=not args.no_inconsistent,
        region=args.region,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"edges": len(out.cells[0].data), **report})
    elif not args.quiet:
        print(
            f"feature edges ({args.angle:g} degrees): {len(out.cells[0].data)} written"
        )
        print(f"  sharp:        {report['num_feature']}")
        print(f"  boundary:     {report['num_boundary']}")
        print(f"  non-manifold: {report['num_non_manifold']}")
        print(f"  inconsistent: {report['num_inconsistent']}")
    return 0
