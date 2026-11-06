from .._helpers import read, write
from .._subdivide import subdivide
from ._common import add_io_args


def add_args(parser):
    add_io_args(parser, "subdivided mesh to be written to")
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
