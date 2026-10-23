from .._helpers import _writer_map, read, reader_map, write
from .._mesh_cad import mesh_cad_surface
from ._json import emit_json


def add_args(parser):
    parser.add_argument(
        "infile",
        type=str,
        help="CAD model to be read from (Kratos CAD JSON, STEP, IGES)",
    )
    parser.add_argument("outfile", type=str, help="surface mesh to be written to")
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
        "--size",
        "-s",
        type=float,
        required=True,
        help="target edge length h (required, > 0)",
    )
    parser.add_argument(
        "--deflection",
        type=float,
        default=0.0,
        help="largest distance of a boundary segment from its edge (0: off)",
    )
    parser.add_argument(
        "--min-points-per-edge",
        type=int,
        default=0,
        help="at least this many points inside every edge (default: 0)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.0,
        help="geometric tolerance (default: 0, 1e-7 of the control nets' diagonal)",
    )
    parser.add_argument(
        "--no-face-ids",
        action="store_true",
        help="no cad:face cell data and no face:<k> regions",
    )
    parser.add_argument(
        "--no-edges", action="store_true", help="no line block of the edges"
    )
    parser.add_argument(
        "--no-match",
        action="store_true",
        help="do not pair trims that no edge references",
    )
    parser.add_argument(
        "--strict", action="store_true", help="fail when any face cannot be meshed"
    )
    parser.add_argument(
        "--max-points",
        type=int,
        default=50_000_000,
        help="the most points in the result (default: 50000000)",
    )
    parser.add_argument(
        "--min-angle",
        type=float,
        default=0.0,
        help="refine until every angle on the surface is at least this many "
        "degrees, under 60 (default: 0, no refinement)",
    )
    parser.add_argument(
        "--smooth",
        type=int,
        default=0,
        help="passes of guarded smoothing (default: 0)",
    )
    parser.add_argument(
        "--smooth-method",
        choices=("laplacian", "odt"),
        default="laplacian",
        help="how smoothing moves a point (default: laplacian)",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def mesh_cad_surface_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = mesh_cad_surface(
        mesh,
        args.size,
        max_chordal_deviation=args.deflection,
        min_points_per_edge=args.min_points_per_edge,
        tolerance=args.tolerance,
        record_face_ids=not args.no_face_ids,
        edges=not args.no_edges,
        match_unreferenced_trims=not args.no_match,
        strict=args.strict,
        max_points=args.max_points,
        min_angle=args.min_angle,
        smooth_iterations=args.smooth,
        smooth_method=args.smooth_method,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(report)
    elif not args.quiet:
        print("mesh_cad_surface")
        print(
            f"  faces:          {report['num_faces']} ({report['num_failed']} failed)"
        )
        print(f"  points:         {report['num_points']}")
        print(f"  triangles:      {report['num_triangles']}")
        print(f"  edge segments:  {report['num_edge_segments']}")
        print(f"  open edges:     {report['num_open_edges']}")
        print(
            f"  angles:         {report['min_angle_deg']:.2f} to "
            f"{report['max_angle_deg']:.2f} degrees"
        )
        print(f"  edge lengths:   {report['min_edge']:g} to {report['max_edge']:g}")
        if args.min_angle > 0.0 or args.smooth > 0:
            print(
                f"  refinement:     {report['num_refinement_points']} points, "
                f"{report['num_smoothed']} smoothing moves, "
                f"{report['num_faces_below_bound']} faces below {args.min_angle:g} degrees"
            )
        for f in report["faces"]:
            if f["status"] != "ok":
                print(f"  face {f['face']}: {f['status']}: {f['reason']}")
            for c in f.get("corners", []):
                what = (
                    "pole" if c["pole"] else f"corner of {c['angle_deg']:.2f} degrees"
                )
                print(
                    f"  face {f['face']}: {what} at point {c['point']} keeps "
                    f"triangles down to {c['min_forced_deg']:.2f} degrees"
                )
    return 0
