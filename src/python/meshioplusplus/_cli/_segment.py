from .._helpers import _writer_map, read, reader_map, write
from .._segment import segment
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="surface mesh to be read from")
    parser.add_argument(
        "outfile", type=str, help="segmented mesh (with patch regions) to be written to"
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
        help="largest dihedral angle still crossed, 0-180 (default: 30)",
    )
    parser.add_argument(
        "--normal-tolerance",
        type=float,
        default=0.0,
        metavar="DEG",
        help="largest angle between a cell's normal and its seed's (default: off)",
    )
    parser.add_argument(
        "--curvature-tolerance",
        type=float,
        default=0.0,
        metavar="K",
        help="largest curvature-estimate difference from the patch mean, 1/length "
        "(default: off)",
    )
    parser.add_argument(
        "--min-patch-cells",
        type=int,
        default=1,
        metavar="N",
        help="merge patches with fewer cells into a neighbour (default: 1)",
    )
    parser.add_argument(
        "--region",
        type=str,
        default="",
        help="restrict to this named cell region (default: every cell)",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="patch",
        help="patch regions are named PREFIX:k (default: patch)",
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=["angle", "primitive"],
        default="angle",
        help="end patches by angles alone, or where the plane, sphere, cylinder, "
        "cone or torus the facets lie on changes (default: angle)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-3,
        metavar="T",
        help="primitive method: how far a point may lie from its patch's "
        "primitive, relative to the diagonal (default: 1e-3)",
    )
    parser.add_argument("--json", action="store_true", help="emit the counts as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def segment_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = segment(
        mesh,
        feature_angle=args.angle,
        normal_tolerance=args.normal_tolerance,
        curvature_tolerance=args.curvature_tolerance,
        min_patch_cells=args.min_patch_cells,
        region=args.region,
        prefix=args.prefix,
        return_report=True,
        method=args.method,
        tolerance=args.tolerance,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(report)
    elif not args.quiet:
        how = "primitive" if args.method == "primitive" else f"{args.angle:g} degrees"
        print(f"segment ({how}): {report['num_patches']} patch(es)")
        print(f"  merged:        {report['num_merged']}")
        print(f"  feature edges: {report['num_feature_edges']}")
        print(f"  open edges:    {report['num_boundary_edges']}")
    return 0
