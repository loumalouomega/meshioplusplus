from .._helpers import _writer_map, read, reader_map, write
from .._interpolate import interpolate


def add_args(parser):
    parser.add_argument("source", type=str, help="mesh file whose data is sampled")
    parser.add_argument("target", type=str, help="mesh file receiving the samples")
    parser.add_argument("outfile", type=str, help="interpolated mesh to be written to")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format (both source and target)",
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
        "--method",
        type=str,
        choices=["nearest", "barycentric"],
        default="nearest",
        help=(
            "how target samples draw their values: nearest source point "
            "(default; dtype-preserving) or barycentric in a simplexified "
            "source (linear-exact; Float64). cell_data always transfers by "
            "nearest source-cell centroid"
        ),
    )
    parser.add_argument(
        "--arrays",
        type=str,
        default=None,
        help=(
            "comma-separated source array names to transfer "
            "(default: every source point_data array; cell_data transfers "
            "only when named explicitly)"
        ),
    )
    parser.add_argument(
        "--extrapolate",
        action="store_true",
        help=(
            "barycentric only: give a target point outside the source domain "
            "the nearest source point's value instead of --default-value"
        ),
    )
    parser.add_argument(
        "--default-value",
        type=float,
        default=0.0,
        help=(
            "barycentric only: fill value for target points outside the "
            "source domain; negative values need the --default-value=-1 form "
            "(default: 0)"
        ),
    )
    parser.add_argument(
        "--on-conflict",
        type=str,
        choices=["error", "overwrite", "suffix"],
        default="error",
        help=(
            "what to do when a transferred name already exists on the target: "
            "error (default), overwrite, or suffix (writes to NAME_interp)"
        ),
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="do not print the transfer summary",
    )


def interpolate_cmd(args):
    source = read(args.source, file_format=args.input_format)
    target = read(args.target, file_format=args.input_format)
    arrays = None
    if args.arrays is not None:
        arrays = [a.strip() for a in args.arrays.split(",") if a.strip()]
    out = interpolate(
        source,
        target,
        method=args.method,
        arrays=arrays,
        extrapolate=args.extrapolate,
        default_value=args.default_value,
        on_conflict=args.on_conflict,
    )
    if not args.quiet:
        n = len(arrays) if arrays is not None else len(source.point_data)
        print(f"interpolated {n} array(s) onto {len(out.points)} target points")
    write(args.outfile, out, file_format=args.output_format)
