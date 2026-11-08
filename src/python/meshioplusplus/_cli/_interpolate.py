from .._helpers import write
from .._interpolate import interpolate
from ._common import (
    add_on_conflict_arg,
    add_transfer_io_args,
    add_transfer_quiet_arg,
    read_source_target,
    split_csv,
)


def add_args(parser):
    add_transfer_io_args(parser, "interpolated mesh to be written to")
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
    add_on_conflict_arg(parser)
    add_transfer_quiet_arg(parser)


def interpolate_cmd(args):
    source, target = read_source_target(args)
    arrays = None if args.arrays is None else split_csv(args.arrays)
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
