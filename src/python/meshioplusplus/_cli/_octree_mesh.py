from .._helpers import read, write
from .._octree_mesh import (
    CONFORMITIES,
    FILLS,
    HEX_GRADINGS,
    OUTPUTS,
    SKIN_GRADINGS,
    octree_mesh,
)
from ._common import add_box_args, add_io_args, parse_bounds
from ._json import emit_json


def add_args(parser):
    add_io_args(
        parser, "volume mesh to be written to", in_help="skin mesh file to be read from"
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
    add_box_args(parser, bounds_of="padded", padding_relative=0.1)
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
        help="how tetrahedra conform to the skin; cleaving: one material per nested "
        "shell (tetrahedral and hexahedral outputs)",
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
        "--hex-grading",
        type=str,
        default="uniform",
        choices=list(HEX_GRADINGS),
        help="hexahedral output: one grid at the skin's size, or the octree's grading",
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
    parser.add_argument(
        "--no-volume-fix",
        action="store_true",
        help="skip the volume fix of the tetrahedral output",
    )
    parser.add_argument(
        "--zero-volume-factor",
        type=float,
        default=1.0e-6,
        help="degenerate below this fraction of the longest edge cubed",
    )
    parser.add_argument(
        "--cap-flatness",
        type=float,
        default=0.01,
        help="a flat tetrahedron is a cap under this height over its face",
    )
    parser.add_argument(
        "--crossing-margin",
        type=float,
        default=0.01,
        help="a sliver's edges cross this far inside both, as a fraction",
    )
    parser.add_argument(
        "--crossing-gap",
        type=float,
        default=0.05,
        help="... and come this close, relative to the shorter edge",
    )
    parser.add_argument(
        "--post-adapt",
        action="store_true",
        help="remesh the tetrahedra with adapt_volume (doc/adapt_volume.md)",
    )
    parser.add_argument("--adapt-hmin", type=float, default=0.0)
    parser.add_argument("--adapt-hmax", type=float, default=0.0)
    parser.add_argument("--adapt-hgrad", type=float, default=1.3)
    parser.add_argument("--adapt-hausd", type=float, default=0.0)
    parser.add_argument("--adapt-max-iterations", type=int, default=10)
    parser.add_argument("--adapt-preserve-boundary", action="store_true")
    parser.add_argument(
        "--adapt-boundary-sizing", choices=("input", "geometry"), default="input"
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def octree_mesh_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    bounds = parse_bounds(args.bounds, "octree-mesh")

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
        hex_grading=args.hex_grading,
        untangle_iterations=args.untangle_iterations,
        max_cells=args.max_cells,
        max_tets=args.max_tets,
        volume_fix=not args.no_volume_fix,
        zero_volume_factor=args.zero_volume_factor,
        cap_flatness=args.cap_flatness,
        crossing_margin=args.crossing_margin,
        crossing_gap=args.crossing_gap,
        post_adapt=args.post_adapt,
        adapt_hmin=args.adapt_hmin,
        adapt_hmax=args.adapt_hmax,
        adapt_hgrad=args.adapt_hgrad,
        adapt_hausd=args.adapt_hausd,
        adapt_max_iterations=args.adapt_max_iterations,
        adapt_preserve_boundary=args.adapt_preserve_boundary,
        adapt_boundary_sizing=args.adapt_boundary_sizing,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)

    if args.json:
        # Timings differ from run to run; the JSON holds what is reproducible.
        doc = {k: v for k, v in report.items() if k != "timings"}
        if "adapt" in doc:
            doc["adapt"] = {k: v for k, v in doc["adapt"].items() if k != "timings"}
        emit_json(doc)
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
