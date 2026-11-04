from .._helpers import _reader_names, _writer_names, read, write
from .._surface import extract_surface


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="surface mesh file to be written to")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=_reader_names(),
        help="input file format",
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
        "--parent-ids",
        "-p",
        action="store_true",
        help="record each boundary facet's parent cell id as cell_data",
    )


def extract_surface_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    surface = extract_surface(mesh, record_parent_ids=args.parent_ids)
    write(args.outfile, surface, file_format=args.output_format)
    return 0
