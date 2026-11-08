from .._conservative_interpolate import conservative_interpolate
from .._helpers import write
from ._common import (
    add_on_conflict_arg,
    add_transfer_io_args,
    add_transfer_quiet_arg,
    read_source_target,
    split_csv,
)


def add_args(parser):
    add_transfer_io_args(parser, "conservatively remapped mesh to be written to")
    parser.add_argument(
        "--arrays",
        type=str,
        default=None,
        help=(
            "comma-separated source array names to transfer "
            "(default: every source point_data AND cell_data array)"
        ),
    )
    parser.add_argument(
        "--default-value",
        type=float,
        default=0.0,
        help=(
            "fill value for a target cell with no source overlap; negative "
            "values need the --default-value=-1 form (default: 0)"
        ),
    )
    add_on_conflict_arg(parser)
    add_transfer_quiet_arg(parser)


def conservative_interpolate_cmd(args):
    source, target = read_source_target(args)
    arrays = None if args.arrays is None else split_csv(args.arrays)
    out = conservative_interpolate(
        source,
        target,
        arrays=arrays,
        default_value=args.default_value,
        on_conflict=args.on_conflict,
    )
    if not args.quiet:
        n = (
            len(arrays)
            if arrays is not None
            else len(source.point_data) + len(source.cell_data)
        )
        print(
            f"conservatively interpolated {n} array(s) onto "
            f"{len(out.points)} target points"
        )
    write(args.outfile, out, file_format=args.output_format)
