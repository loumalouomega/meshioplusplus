from .._helpers import read, write
from .._optimize_volume import optimize_volume
from ._common import add_io_args


def add_args(parser):
    add_io_args(
        parser,
        "ODT-remeshed mesh to be written to",
        in_help="tetrahedral mesh to be read from",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=10,
        help="optimisation sweeps (relocation + flips); stops early at a fixed point (default: 10)",
    )
    parser.add_argument(
        "--no-relocate",
        action="store_true",
        help="skip the ODT vertex-relocation half (flips only)",
    )
    parser.add_argument(
        "--no-flip",
        action="store_true",
        help="skip the topological-flip half (reduces to ODT smoothing)",
    )
    parser.add_argument(
        "--no-preserve-boundary",
        action="store_true",
        help="allow boundary vertices to move during relocation (may drift off the surface)",
    )
    parser.add_argument(
        "--min-improvement",
        type=float,
        default=1e-6,
        help="strict gain a move must deliver: degrees of dihedral angle under "
        "--objective worst, scaled Jacobian under energy (default: 1e-6)",
    )
    parser.add_argument(
        "--objective",
        type=str,
        choices=["worst", "energy"],
        default="worst",
        help="worst (default): raise the least dihedral angle, never lowering it; "
        "energy: the ODT energy with scaled-Jacobian flips (v18.0.0 and before)",
    )
    parser.add_argument(
        "--max-ring-size",
        type=int,
        default=7,
        help="worst only: the largest ring an edge or multi-face removal takes (default: 7)",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def optimize_volume_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    out, report = optimize_volume(
        mesh,
        max_iterations=args.max_iterations,
        relocate=not args.no_relocate,
        flip=not args.no_flip,
        preserve_boundary=not args.no_preserve_boundary,
        min_improvement=args.min_improvement,
        return_report=True,
        objective=args.objective,
        max_ring_size=args.max_ring_size,
    )

    if not args.quiet:
        print("ODT-remeshed")
        print(
            f"  flips (2-3 / 3-2):        {report['num_flips']} "
            f"({report['num_23_flips']} / {report['num_32_flips']})"
        )
        if args.objective == "worst":
            print(
                f"  edge / multi-face removals: {report['num_edge_removals']} / "
                f"{report['num_multi_face_removals']}"
            )
            print(f"  worst-cell smoothing moves: {report['num_worst_smoothed']}")
        print(f"  vertices moved:           {report['num_vertices_moved']}")
        print(f"  tets:                     {report['num_tets']}")
        print(
            f"  min quality before/after: {report['min_quality_before']:.4f} / "
            f"{report['min_quality_after']:.4f}"
        )
        print(
            f"  min dihedral before/after: {report['min_dihedral_before']:.2f} / "
            f"{report['min_dihedral_after']:.2f} degrees"
        )

    write(args.outfile, out, file_format=args.output_format)
    return 0
