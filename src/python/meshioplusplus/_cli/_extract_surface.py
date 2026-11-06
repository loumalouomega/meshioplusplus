from .._helpers import read, write
from .._surface import extract_surface
from ._common import add_io_args


def add_args(parser):
    add_io_args(parser, "surface mesh file to be written to")
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
