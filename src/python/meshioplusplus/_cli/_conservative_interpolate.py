from .._conservative_interpolate import conservative_interpolate
from .._helpers import _writer_map, read, reader_map, write


def add_args(parser):
    parser.add_argument("source", type=str, help="mesh file whose data is sampled")
    parser.add_argument("target", type=str, help="mesh file receiving the samples")
    parser.add_argument(
        "outfile", type=str, help="conservatively remapped mesh to be written to"
    )
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


def conservative_interpolate_cmd(args):
    source = read(args.source, file_format=args.input_format)
    target = read(args.target, file_format=args.input_format)
    arrays = None
    if args.arrays is not None:
        arrays = [a.strip() for a in args.arrays.split(",") if a.strip()]
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
