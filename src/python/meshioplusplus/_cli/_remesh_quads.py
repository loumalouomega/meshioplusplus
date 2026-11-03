import json

from .._helpers import _writer_map, read, reader_map, write
from .._remesh_quads import remesh_quads


def add_args(parser):
    parser.add_argument("infile", type=str, help="triangle surface")
    parser.add_argument("outfile", type=str, help="quads to be written to")
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
        "--faces",
        "-f",
        type=int,
        default=0,
        help="the number of faces asked for (default 0: half the triangles)",
    )
    parser.add_argument(
        "--backend",
        type=str,
        default="auto",
        metavar="NAME",
        help="the quad remesher: auto (QuadriFlow, else Instant Meshes), quadriflow "
        "or instant_meshes, the installed programs (default auto); "
        "remesh_quads_backends() lists the names this install takes",
    )
    parser.add_argument(
        "--sharp", action="store_true", help="align with and keep the sharp edges"
    )
    parser.add_argument(
        "--no-boundary",
        action="store_true",
        help="do not keep the boundary of an open surface",
    )
    parser.add_argument(
        "--adaptive",
        action="store_true",
        help="let the face size follow the curvature (QuadriFlow)",
    )
    parser.add_argument(
        "--dominant",
        action="store_true",
        help="allow triangles where quads fit badly (Instant Meshes)",
    )
    parser.add_argument(
        "--seed", type=int, default=0, help="the program's random seed (default 0)"
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=50_000_000,
        help="refuse a target or result of more faces (default 50000000)",
    )
    parser.add_argument("--json", action="store_true", help="print the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="do not print the summary"
    )


def remesh_quads_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = remesh_quads(
        mesh,
        target_faces=args.faces,
        backend=args.backend,
        sharp=args.sharp,
        boundary=not args.no_boundary,
        adaptive=args.adaptive,
        dominant=args.dominant,
        seed=args.seed,
        max_cells=args.max_cells,
        return_report=True,
    )
    if args.json:
        print(json.dumps(report))
    elif not args.quiet:
        print(
            f"{report['num_quads']} quads, {report['num_triangles']} triangles, "
            f"{report['num_points']} points; scaled Jacobian "
            f"{report['min_scaled_jacobian']:.3f} (mean "
            f"{report['mean_scaled_jacobian']:.3f}), "
            f"{report['hausdorff_to_surface']:.6g} from the surface; {report['backend']}"
        )
    write(args.outfile, out, file_format=args.output_format)
