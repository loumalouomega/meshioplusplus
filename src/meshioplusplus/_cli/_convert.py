import numpy as np

from .._helpers import _writer_map, read, reader_map, write


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
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
        "--ascii",
        "-a",
        action="store_true",
        help="write in ASCII format variant (where applicable, default: binary)",
    )
    parser.add_argument("outfile", type=str, help="mesh file to be written to")
    parser.add_argument(
        "--float-format",
        "-f",
        type=str,
        help="float format used in output ASCII files (default: .16e)",
    )
    parser.add_argument(
        "--sets-to-int-data",
        "-s",
        action="store_true",
        help="if possible, convert sets to integer data (useful if the output type does not support sets)",
    )
    parser.add_argument(
        "--int-data-to-sets",
        "-d",
        action="store_true",
        help="if possible, convert integer data to sets (useful if the output type does not support integer data)",
    )
    parser.add_argument(
        "--points-only",
        action="store_true",
        help=(
            "read geometry only, dropping every data array. Readers with a "
            "native selective path skip the arrays outright rather than "
            "reading and discarding them."
        ),
    )
    parser.add_argument(
        "--arrays",
        type=str,
        default=None,
        help=(
            "comma-separated data-array names to keep; everything else is "
            "skipped. Pass an empty string to keep none. Names absent from the "
            "file are ignored."
        ),
    )


def convert(args):
    if args.points_only and args.arrays is not None:
        raise ValueError("--points-only and --arrays are mutually exclusive")

    arrays = None
    if args.arrays is not None:
        arrays = [name for name in args.arrays.split(",") if name]

    # read mesh data
    mesh = read(
        args.infile,
        file_format=args.input_format,
        points_only=args.points_only,
        arrays=arrays,
    )

    # Some converters (like VTK) require `points` to be contiguous.
    mesh.points = np.ascontiguousarray(mesh.points)

    if (args.points_only or arrays is not None) and (
        args.sets_to_int_data or args.int_data_to_sets
    ):
        raise ValueError(
            "--points-only/--arrays cannot be combined with -s/-d: the set "
            "conversions operate on the very data arrays that were skipped"
        )

    if args.sets_to_int_data:
        mesh.point_sets_to_data()
        mesh.cell_sets_to_data()

    if args.int_data_to_sets:
        # Snapshot the keys: *_data_to_sets(key) deletes that key from the dict,
        # so iterating the live view raises "dictionary changed size".
        for key in list(mesh.point_data):
            mesh.point_data_to_sets(key)
        for key in list(mesh.cell_data):
            mesh.cell_data_to_sets(key)

    # write it out
    kwargs = {"file_format": args.output_format}
    if args.float_format is not None:
        kwargs["float_fmt"] = args.float_format
    if args.ascii:
        kwargs["binary"] = False

    write(args.outfile, mesh, **kwargs)
