from .._agglomerate import agglomerate
from .._helpers import _writer_map, read, reader_map, write


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="coarsened mesh to be written to")
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
        "--target-group-size",
        type=int,
        default=8,
        help="approximate member cells per output group (default: 8)",
    )


def agglomerate_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out = agglomerate(mesh, target_group_size=args.target_group_size)
    write(args.outfile, out, file_format=args.output_format)
    return 0
