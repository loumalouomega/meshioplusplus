from .._helpers import _writer_map, read, reader_map, write
from .._undo_green import undo_green


def add_args(parser):
    parser.add_argument(
        "coarse", type=str, help="mesh a prior refine() call was run on"
    )
    parser.add_argument("fine", type=str, help="that call's output")
    parser.add_argument("outfile", type=str, help="undone mesh to be written to")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format (both coarse and fine)",
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
        "--quiet",
        "-q",
        action="store_true",
        help="do not print the undo summary",
    )


def undo_green_cmd(args):
    coarse = read(args.coarse, file_format=args.input_format)
    fine = read(args.fine, file_format=args.input_format)
    out, report = undo_green(coarse, fine, return_report=True)
    if not args.quiet:
        print(
            f"undid {report['num_groups_undone']} green group(s), removing "
            f"{report['num_cells_removed']} cell(s)"
        )
    write(args.outfile, out, file_format=args.output_format)
