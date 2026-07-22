from .._helpers import _writer_map, read, reader_map, write
from .._refine import refine


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="refined mesh to be written to")
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
        "--levels",
        type=int,
        default=1,
        help=(
            "how many times to subdivide; each level splits a triangle/quad "
            "into 4 and a tetra/wedge/hexahedron into 8 (default: 1)"
        ),
    )
    parser.add_argument(
        "--record-parent-ids",
        action="store_true",
        help="attach refine:parent_cell cell_data of the original cell indices",
    )


def refine_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out = refine(mesh, levels=args.levels, record_parent_ids=args.record_parent_ids)
    write(args.outfile, out, file_format=args.output_format)
    return 0
