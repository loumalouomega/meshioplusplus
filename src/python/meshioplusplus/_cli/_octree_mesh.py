from .._helpers import _writer_map, read, reader_map, write
from .._octree_mesh import CONFORMITIES, FILLS, OUTPUTS, SKIN_GRADINGS, octree_mesh
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="skin mesh file to be read from")
    parser.add_argument("outfile", type=str, help="volume mesh to be written to")
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
        "--output",
        type=str,
        default="tet",
        choices=list(OUTPUTS),
        help="castellated hexahedra, conforming tetrahedra (default) or all-hexahedra",
    )
    parser.add_argument(
        "--min-size",
        type=float,
        default=0.0,
        help="smallest leaf edge (0: max-size / 32)",
    )
    parser.add_argument(
        "--max-size",
        type=float,
        default=0.0,
        help="largest leaf edge (0: an eighth of the skin's diagonal)",
    )
    parser.add_argument(
        "--max-level",
        type=int,
        default=-1,
        help="cap on the finest level (-1: from the sizes)",
    )
    parser.add_argument(
        "--bounds",
        type=str,
        default=None,
        help="'xlo,ylo,zlo,xhi,yhi,zhi'; the padded bounding box by default "
        "(negatives need --bounds=)",
    )
    parser.add_argument(
        "--padding", type=float, default=0.0, help="grow the box by this on every side"
    )
    parser.add_argument(
        "--padding-relative",
        type=float,
        default=0.1,
        help="grow the box by this fraction of its diagonal on every side",
    )
    parser.add_argument(
        "--curvature-angle",
        type=float,
        default=15.0,
        help="degrees a leaf edge may subtend on the skin (0: off)",
    )
    parser.add_argument(
        "--feature-angle",
        type=float,
        default=30.0,
        help="edges sharper than this many degrees are features (0: off)",
    )
    parser.add_argument(
        "--corner-angle",
        type=float,
        default=30.0,
        help="a feature turning by more than this many degrees has a corner",
    )
    parser.add_argument(
        "--feature-size",
        type=float,
        default=0.0,
        help="size along features (0: min-size)",
    )
    parser.add_argument(
        "--gap-cells",
        type=float,
        default=2.0,
        help="leaves across the local thickness (0: off)",
    )
    parser.add_argument(
        "--gap-outside", action="store_true", help="also resolve gaps outside the skin"
    )
    parser.add_argument(
        "--size-array", type=str, default=None, help="point_data array of target sizes"
    )
    parser.add_argument(
        "--surface-region",
        type=str,
        default=None,
        help="mesh only this cell region's skin",
    )
    parser.add_argument(
        "--grading-layers",
        type=int,
        default=1,
        help="rings of leaves around the skin at its level",
    )
    parser.add_argument(
        "--skin-grading",
        type=str,
        default="graded",
        choices=list(SKIN_GRADINGS),
        help="let the size vary along the skin, or keep it uniform (angle guarantee)",
    )
    parser.add_argument(
        "--surface-size",
        type=float,
        default=0.0,
        help="skin size under --skin-grading uniform (0: smallest asked for)",
    )
    parser.add_argument(
        "--winding-threshold",
        type=float,
        default=0.5,
        help="inside where the winding number exceeds this",
    )
    parser.add_argument(
        "--winding-beta",
        type=float,
        default=2.0,
        help="fast winding number accuracy (0: exact)",
    )
    parser.add_argument(
        "--conformity",
        type=str,
        default="stuffing",
        choices=list(CONFORMITIES),
        help="how tetrahedra conform to the skin",
    )
    parser.add_argument(
        "--alpha-long",
        type=float,
        default=0.24799,
        help="stuffing threshold, long edges",
    )
    parser.add_argument(
        "--alpha-short",
        type=float,
        default=0.40173,
        help="stuffing threshold, short edges",
    )
    parser.add_argument(
        "--no-preserve-features",
        action="store_true",
        help="do not snap the boundary to feature corners and edges",
    )
    parser.add_argument(
        "--fit-iterations",
        type=int,
        default=5,
        help="rounds of fitting and optimisation (0: off)",
    )
    parser.add_argument(
        "--no-dihedral-floor-guard",
        action="store_true",
        help="let optimisation lower the minimum dihedral",
    )
    parser.add_argument(
        "--fill",
        type=str,
        default="inside",
        choices=list(FILLS),
        help="castellated output: all leaves, the inside ones, or inside and cut",
    )
    parser.add_argument(
        "--keep-outside",
        action="store_true",
        help="tetrahedral output: also mesh the outside as the 'outside' region",
    )
    parser.add_argument(
        "--no-distance",
        action="store_true",
        help="do not attach sdf:distance / sdf:closest_point",
    )
    parser.add_argument(
        "--attach-quality", action="store_true", help="attach the quality:* cell data"
    )
    parser.add_argument(
        "--no-pillow", action="store_true", help="hexahedral output: no boundary layer"
    )
    parser.add_argument(
        "--untangle-iterations",
        type=int,
        default=50,
        help="hexahedral output: untangling iterations",
    )
    parser.add_argument(
        "--max-cells", type=int, default=20000000, help="refuse above this many leaves"
    )
    parser.add_argument(
        "--max-tets",
        type=int,
        default=100000000,
        help="refuse above this many output cells",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def octree_mesh_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    bounds = None
    if args.bounds is not None:
        bounds = [float(x) for x in args.bounds.split(",")]
        if len(bounds) != 6:
            raise ValueError("octree-mesh: --bounds expects 'xlo,ylo,zlo,xhi,yhi,zhi'")

    out, report = octree_mesh(
        mesh,
        output=args.output,
        min_size=args.min_size,
        max_size=args.max_size,
        max_level=args.max_level,
        bounds=bounds,
        padding=args.padding,
        padding_relative=args.padding_relative,
        curvature_angle=args.curvature_angle,
        feature_angle=args.feature_angle,
        corner_angle=args.corner_angle,
        feature_size=args.feature_size,
        gap_cells=args.gap_cells,
        gap_outside=args.gap_outside,
        size_array=args.size_array,
        surface_region=args.surface_region,
        grading_layers=args.grading_layers,
        skin_grading=args.skin_grading,
        surface_size=args.surface_size,
        winding_threshold=args.winding_threshold,
        winding_beta=args.winding_beta,
        conformity=args.conformity,
        alpha_long=args.alpha_long,
        alpha_short=args.alpha_short,
        preserve_features=not args.no_preserve_features,
        fit_iterations=args.fit_iterations,
        dihedral_floor_guard=not args.no_dihedral_floor_guard,
        fill=args.fill,
        keep_outside=args.keep_outside,
        attach_distance=not args.no_distance,
        attach_quality=args.attach_quality,
        pillow=not args.no_pillow,
        untangle_iterations=args.untangle_iterations,
        max_cells=args.max_cells,
        max_tets=args.max_tets,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)

    if args.json:
        # Timings differ from run to run; the JSON holds what is reproducible.
        emit_json({k: v for k, v in report.items() if k != "timings"})
    elif not args.quiet:
        print(f"octree mesh ({report['output']})")
        print(f"  root size:      {report['root_size']:g}")
        print(f"  levels:         {report['finest_level']} of {report['max_level']}")
        print(f"  leaves:         {report['num_leaves']}")
        print(f"  cells:          {report['num_cells']}")
        if report["output"] == "castellated":
            print(f"  inside:         {report['num_inside']}")
            print(f"  surrogate faces:{report['num_surrogate_faces']:>8}")
            print(f"  hanging nodes:  {report['num_hanging']}")
        else:
            print(
                f"  dihedral:       {report['min_dihedral']:.2f} .. "
                f"{report['max_dihedral']:.2f} deg"
            )
        print(f"  min scaled Jacobian: {report['min_scaled_jacobian']:.4f}")
        for w in report["warnings"]:
            print(f"  warning: {w}")
    return 0
