from .._helpers import _writer_map, read, reader_map, write
from .._optimize_volume import optimize_volume


def add_args(parser):
    parser.add_argument("infile", type=str, help="tetrahedral mesh to be read from")
    parser.add_argument("outfile", type=str, help="ODT-remeshed mesh to be written to")
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
        help="strict scaled-Jacobian gain a flip must deliver to be accepted (default: 1e-6)",
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
    )

    if not args.quiet:
        print("ODT-remeshed")
        print(
            f"  flips (2-3 / 3-2):        {report['num_flips']} "
            f"({report['num_23_flips']} / {report['num_32_flips']})"
        )
        print(f"  vertices moved:           {report['num_vertices_moved']}")
        print(f"  tets:                     {report['num_tets']}")
        print(
            f"  min quality before/after: {report['min_quality_before']:.4f} / "
            f"{report['min_quality_after']:.4f}"
        )

    write(args.outfile, out, file_format=args.output_format)
    return 0
