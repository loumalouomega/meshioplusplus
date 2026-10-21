from .._fit_cad import DEFAULT_RANSAC_SEED
from .._fit_primitives import TYPES
from .._helpers import read, reader_map, write
from .._refine_patches import refine_patches
from .._segment import segment
from ._json import emit_json


def parse_seed(text):
    """A RANSAC seed: decimal, or ``0x`` hexadecimal."""
    return int(text, 0)


def add_args(parser):
    parser.add_argument("infile", type=str, help="surface mesh to be read from")
    parser.add_argument(
        "outfile", type=str, help="mesh with the refined patches to be written to"
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
        "--tolerance",
        type=float,
        default=1e-3,
        metavar="T",
        help="how far a point may lie from its patch's primitive, relative to the "
        "diagonal (default: 1e-3)",
    )
    parser.add_argument(
        "--absolute",
        action="store_true",
        help="the tolerance is a length, not a fraction of the diagonal",
    )
    parser.add_argument(
        "--normal-weight",
        type=float,
        default=1.0,
        metavar="W",
        help="the weight of the normal term in a cell's cost (default: 1)",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=1e-3,
        metavar="E",
        help="stop when an iteration adds less than this to the area share "
        "primitives explain (default: 1e-3)",
    )
    parser.add_argument(
        "--feature-angle",
        type=float,
        default=30.0,
        metavar="DEG",
        help="dihedral angle above which an edge is a feature no patch crosses, "
        "also segment's when the input has no segment:patch (default: 30)",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=10,
        metavar="N",
        help="the most iterations (default: 10)",
    )
    parser.add_argument(
        "--min-patch-cells",
        type=int,
        default=3,
        metavar="N",
        help="a freeform patch with fewer cells joins a neighbour (default: 3)",
    )
    parser.add_argument(
        "--ransac-seed",
        type=parse_seed,
        default=DEFAULT_RANSAC_SEED,
        metavar="S",
        help="the seed of the RANSAC split (default: 0x6d657368696f2b2b)",
    )
    parser.add_argument(
        "--types",
        type=str,
        default=",".join(TYPES[1:]),
        help="comma-separated primitives a patch may be (default: all)",
    )
    parser.add_argument(
        "--no-split",
        action="store_true",
        help="do not split a primitive out of a freeform patch",
    )
    parser.add_argument(
        "--no-merge",
        action="store_true",
        help="do not merge neighbouring patches one primitive explains",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="patch",
        help="patch regions are named PREFIX:k (default: patch)",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def refine_patches_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    if "segment:patch" not in mesh.cell_data:
        mesh = segment(mesh, feature_angle=args.feature_angle)
    out, report = refine_patches(
        mesh,
        tolerance=args.tolerance,
        relative=not args.absolute,
        normal_weight=args.normal_weight,
        epsilon=args.epsilon,
        feature_angle=args.feature_angle,
        max_iterations=args.max_iterations,
        min_patch_cells=args.min_patch_cells,
        ransac_seed=args.ransac_seed,
        types=[t.strip() for t in args.types.split(",") if t.strip()],
        allow_split=not args.no_split,
        allow_merge=not args.no_merge,
        prefix=args.prefix,
        return_report=True,
    )
    write(args.outfile, out)
    if args.json:
        emit_json(report)
    elif not args.quiet:
        print(
            f"refine_patches: {report['num_patches']} patch(es) after "
            f"{report['num_iterations']} iteration(s) ({report['stop']})"
        )
        print(f"  primitive: {report['num_primitive']}")
        print(f"  freeform:  {report['num_freeform']}")
        print(f"  facets:    {report['num_facets']}")
        print(
            f"  splits {report['num_splits']}, merges {report['num_merges']}, "
            f"absorbed {report['num_absorbed']}"
        )
        print(
            f"  explained area {report['initial_explained_area_fraction']:.6g} -> "
            f"{report['explained_area_fraction']:.6g}, Hausdorff "
            f"{report['hausdorff_max']:.6g} (tolerance {report['tolerance']:.6g})"
        )
    return 0
