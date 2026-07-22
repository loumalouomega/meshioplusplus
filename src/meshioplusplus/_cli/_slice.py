from .._helpers import _writer_map, read, reader_map, write
from .._slice import slice as slice_mesh


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="cross-section mesh to be written to")
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
        "--origin",
        type=str,
        default="0,0,0",
        help="a point on the plane 'x,y,z' (negatives need --origin=-1,0,0)",
    )
    parser.add_argument(
        "--normal",
        type=str,
        default="0,0,1",
        help="the plane normal 'x,y,z' (negatives need --normal=0,0,-1)",
    )
    parser.add_argument(
        "--record-parent-ids",
        action="store_true",
        help="attach the originating input cell id as slice:parent_cell cell_data",
    )


def _floats(text):
    return [float(x) for x in text.split(",")]


def slice_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    origin = _floats(args.origin)
    normal = _floats(args.normal)
    if len(origin) != 3:
        raise ValueError("slice: --origin expects 'x,y,z'")
    if len(normal) != 3:
        raise ValueError("slice: --normal expects 'x,y,z'")

    out = slice_mesh(
        mesh, origin=origin, normal=normal, record_parent_ids=args.record_parent_ids
    )

    write(args.outfile, out, file_format=args.output_format)
    return 0
