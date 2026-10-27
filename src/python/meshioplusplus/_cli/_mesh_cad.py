import argparse

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
    parser.add_argument(
        "--order",
        type=int,
        choices=(1, 2, 3),
        default=1,
        help="1: linear triangles; 2 or 3: curved triangle6 or triangle10 whose "
        "nodes lie on the surface and edges (default: 1)",
    )
    parser.add_argument(
        "--curvature-deviation",
        type=float,
        default=0.0,
        help="size curved regions so a triangle sits at most this far from the "
        "surface (default: 0, off)",
    )
    parser.add_argument(
        "--cells-across-gaps",
        type=float,
        default=0.0,
        help="at least this many elements across a thin wall or a gap "
        "(default: 0, off)",
    )
    parser.add_argument(
        "--grading",
        type=float,
        default=0.0,
        help="how fast the size may grow, a ratio above 1 (default: 0, 1.5 when "
        "a size source is on)",
    )
    parser.add_argument(
        "--pole-rings",
        type=int,
        default=0,
        help="rings of points around each pole, 0 to 8 (default: 0, the fan)",
    )
    parser.add_argument(
        "--face-sizes",
        type=_size_pairs,
        default=None,
        help="sizes given per face, K=H pairs comma-separated (0-based faces)",
    )
    parser.add_argument(
        "--edge-sizes",
        type=_size_pairs,
        default=None,
        help="sizes given per edge, K=H pairs comma-separated (0-based edges)",
    )
    parser.add_argument(
        "--background",
        default=None,
        help="a mesh whose point_data --background-array is the size at its points",
    )
    parser.add_argument(
        "--background-array",
        default="size",
        help="the background's size array (default: size)",
    )
    parser.add_argument(
        "--metric-background",
        default=None,
        help="a mesh whose 9-component point_data --metric-array is a field's "
        "Hessian (or whose 1-component one is the field itself): an anisotropic "
        "metric from it sizes the triangles, at any order (roadmap 6.10.2)",
    )
    parser.add_argument(
        "--metric-array",
        default="",
        help="the metric background's Hessian array (e.g. f:hessian) or field",
    )
    parser.add_argument(
        "--metric-error",
        type=float,
        default=0.0,
        help="the target interpolation error of --metric-background (positive)",
    )
    parser.add_argument(
        "--metric",
        choices=("mean", "varying"),
        default="mean",
        help="each face's mean metric, or its metric as it varies (default: mean)",
    )
    parser.add_argument(
        "--virtual-edges",
        action="store_true",
        help="sample chains of short edges meeting smoothly as one",
    )
    parser.add_argument(
        "--reparam-condition-number",
        type=float,
        default=0.0,
        metavar="K",
        help="mesh a face whose metric's condition number passes K in a "
        "flattening of its coarse mesh (default: 0, off)",
    )
    parser.add_argument(
        "--model-tolerances",
        action="store_true",
        help="accept the gaps the model itself declares",
    )
    parser.add_argument(
        "--virtual-faces",
        action="store_true",
        help="mesh a sliver face with the neighbour it meets tangentially",
    )
    parser.add_argument(
        "--lop-corners",
        action="store_true",
        help="lop acute corners (Pav and Walkington 2005) in place of "
        "concentric shells, and state the angle bound they prove",
    )
    parser.add_argument(
        "--untangle",
        action="store_true",
        help="with --order 2 or 3, move a face's own nodes (quads' too) to "
        "untangle a folded curved element before it is straightened (roadmap 6.6.3)",
    )
    parser.add_argument(
        "--quads",
        choices=("none", "dominant", "all"),
        default="none",
        help="recombine the triangles into a quad-dominant or all-quad mesh "
        "(roadmap 6.8; with --order 2 or 3 writes quad9 or quad16; all takes no "
        "--virtual-edges)",
    )
    parser.add_argument(
        "--quad-pairs",
        choices=("matching", "greedy"),
        default="matching",
        help="pair the triangles by a maximum matching (the default) or "
        "greedily (roadmap 6.8.1)",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def _size_pairs(text):
    """``"0=0.1,3=0.2"`` as ``[(0, 0.1), (3, 0.2)]``."""
    pairs = []
    for item in text.split(","):
        index, sep, size = item.partition("=")
        if not sep:
            raise argparse.ArgumentTypeError("K=H pairs, comma-separated")
        pairs.append((int(index), float(size)))
    return pairs


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
        order=args.order,
        curvature_deviation=args.curvature_deviation,
        cells_across_gaps=args.cells_across_gaps,
        grading=args.grading,
        pole_rings=args.pole_rings,
        face_sizes=args.face_sizes,
        edge_sizes=args.edge_sizes,
        background=None if args.background is None else read(args.background),
        background_array=args.background_array,
        metric=args.metric,
        virtual_edges=args.virtual_edges,
        reparam_condition_number=args.reparam_condition_number,
        model_tolerances=args.model_tolerances,
        virtual_faces=args.virtual_faces,
        lop_corners=args.lop_corners,
        metric_background=(
            None if args.metric_background is None else read(args.metric_background)
        ),
        metric_array=args.metric_array,
        metric_error=args.metric_error,
        untangle=args.untangle,
        quads=args.quads,
        quad_pairs=args.quad_pairs,
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
        if args.order > 1:
            untangled = (
                f", {report['num_untangled']} untangled" if args.untangle else ""
            )
            print(
                f"  curved:         order {args.order}, "
                f"{report['num_curved_straightened']} kept straight-sided{untangled}"
            )
        if args.quads != "none":
            print(
                f"  quads:          {report['num_quads']} "
                f"({report['num_leftover_triangles']} triangles left), least scaled "
                f"Jacobian {report['min_quad_scaled_jacobian']:.3f}"
            )
        if args.metric_background is not None:
            print(f"  metric angles:  from {report['min_metric_angle']:.2f} degrees")
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
        for c in report["lopped_corners"]:
            stated = f"bound {c['stated_bound_deg']:.2f} degrees " + (
                "guaranteed" if c["guaranteed"] else f"not stated: {c['reason']}"
            )
            print(
                f"  face {c['face']}: corner of {c['angle_deg']:.2f} degrees at point "
                f"{c['point']} lopped at radius {c['lop_radius']:g}; {stated}"
            )
    return 0
