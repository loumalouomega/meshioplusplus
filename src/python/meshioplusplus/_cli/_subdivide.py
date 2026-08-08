from .._helpers import _writer_map, read, reader_map, write
from .._subdivide import subdivide


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="subdivided mesh to be written to")
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
        "--record-parent-ids",
        action="store_true",
        help="attach subdivide:parent_cell cell_data of the source cell indices",
    )


def subdivide_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out = subdivide(mesh, record_parent_ids=args.record_parent_ids)
    write(args.outfile, out, file_format=args.output_format)
    return 0
