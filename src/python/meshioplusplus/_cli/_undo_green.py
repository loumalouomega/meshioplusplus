from .._helpers import _reader_names, _writer_names, read, write
from .._undo_green import undo_green


def add_args(parser):
    parser.add_argument(
        "meshes",
        nargs="+",
        metavar="MESH",
        help=(
            "COARSE FINE OUTFILE (parents read from COARSE), or FINE OUTFILE "
            "(parents from FINE's refine:cell_id/refine:parent_id hierarchy)"
        ),
    )
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=_reader_names(),
        help="input file format (both coarse and fine)",
        default=None,
    )
    parser.add_argument(
        "--output-format",
        "-o",
        type=str,
        choices=_writer_names(),
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
    if len(args.meshes) not in (2, 3):
        raise SystemExit("undo-green: expects COARSE FINE OUTFILE or FINE OUTFILE")
    *inputs, outfile = args.meshes
    meshes = [read(p, file_format=args.input_format) for p in inputs]
    if len(meshes) == 2:
        out, report = undo_green(meshes[0], meshes[1], return_report=True)
    else:
        out, report = undo_green(meshes[0], return_report=True)
    if not args.quiet:
        print(
            f"undid {report['num_groups_undone']} green group(s), removing "
            f"{report['num_cells_removed']} cell(s)"
        )
    write(outfile, out, file_format=args.output_format)
