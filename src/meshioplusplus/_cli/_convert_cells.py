from .._convert_cells import MODES, convert_cells
from .._helpers import _writer_map, read, reader_map, write


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="converted mesh to be written to")
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
        "--mode",
        type=str,
        choices=list(MODES),
        default="linearize",
        help=(
            "linearize: higher-order cells -> their linear base; "
            "simplexify: decompose into same-dimension simplices; "
            "elevate: linear -> serendipity quadratic (default: linearize)"
        ),
    )
    parser.add_argument(
        "--record-parent-ids",
        action="store_true",
        help="attach convert:parent_cell cell_data of the source cell indices",
    )


def convert_cells_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out = convert_cells(mesh, mode=args.mode, record_parent_ids=args.record_parent_ids)
    write(args.outfile, out, file_format=args.output_format)
    return 0
