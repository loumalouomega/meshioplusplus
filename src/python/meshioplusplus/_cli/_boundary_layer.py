import json

from .._boundary_layer import boundary_layer
from .._helpers import _writer_map, read, reader_map, write


def add_args(parser):
    parser.add_argument("infile", type=str, help="triangle and quad surface")
    parser.add_argument("outfile", type=str, help="layers to be written to")
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
        "--first-height", type=float, required=True, help="the first layer's thickness"
    )
    parser.add_argument(
        "--num-layers", type=int, default=3, help="how many layers (default 3)"
    )
    parser.add_argument(
        "--growth",
        type=float,
        default=1.2,
        help="each layer over the one before (default 1.2)",
    )
    parser.add_argument(
        "--reverse",
        action="store_true",
        help="grow against the faces' normals",
    )
    parser.add_argument(
        "--collision-fraction",
        type=float,
        default=0.45,
        help="of the distance a point's ray hits (default 0.45)",
    )
    parser.add_argument(
        "--smooth",
        type=int,
        default=5,
        help="passes smoothing each point's share over its neighbours (default 5)",
    )
    parser.add_argument(
        "--max-shrink",
        type=int,
        default=6,
        help="halvings of a point's layers before they collapse (default 6)",
    )
    parser.add_argument("--json", action="store_true", help="print the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="do not print the summary"
    )


def boundary_layer_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = boundary_layer(
        mesh,
        args.first_height,
        num_layers=args.num_layers,
        growth=args.growth,
        reverse=args.reverse,
        collision_fraction=args.collision_fraction,
        smooth_iterations=args.smooth,
        max_shrink=args.max_shrink,
        return_report=True,
    )
    if args.json:
        print(json.dumps(report))
    elif not args.quiet:
        print(
            f"{report['num_cells']} cells ({report['num_wedges']} wedges, "
            f"{report['num_hexahedra']} hexahedra, {report['num_pyramids']} pyramids, "
            f"{report['num_tetrahedra']} tetrahedra), least scaled Jacobian "
            f"{report['min_scaled_jacobian']:.3f}; thickness "
            f"{report['min_thickness']:.6g} to {report['max_thickness']:.6g}, "
            f"{report['num_shrunk']} points thinner, {report['num_collapsed']} "
            f"collapsed; cap {'closed' if report['cap_closed'] else 'open'}"
        )
    write(args.outfile, out, file_format=args.output_format)
