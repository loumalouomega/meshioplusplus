from .._helpers import read, write
from .._slice import sections
from .._slice import slice as slice_mesh
from ._common import add_io_args


def add_args(parser):
    add_io_args(parser, "cross-section mesh to be written to")
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
    parser.add_argument(
        "--chain",
        action="store_true",
        help="walk a surface's line section into polylines (section:chain, section:closed)",
    )
    parser.add_argument(
        "--min-length",
        type=float,
        default=0.0,
        help="with --chain, drop chains shorter than this",
    )


def add_sections_args(parser):
    add_io_args(parser, "merged line sections to be written to")
    parser.add_argument(
        "--direction",
        type=str,
        default="0,0,1",
        help="the planes' normal 'x,y,z' (negatives need --direction=0,0,-1)",
    )
    parser.add_argument(
        "--offsets",
        type=str,
        required=True,
        help="the planes' distances along the unit direction 'a,b,...' "
        "(negatives need --offsets=-1,0,1)",
    )
    parser.add_argument(
        "--min-length", type=float, default=0.0, help="drop chains shorter than this"
    )
    parser.add_argument(
        "--record-parent-ids",
        action="store_true",
        help="attach the originating input cell id as slice:parent_cell cell_data",
    )


def sections_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    direction = _floats(args.direction)
    if len(direction) != 3:
        raise ValueError("sections: --direction expects 'x,y,z'")
    out = sections(
        mesh,
        direction,
        _floats(args.offsets),
        min_length=args.min_length,
        record_parent_ids=args.record_parent_ids,
    )
    write(args.outfile, out, file_format=args.output_format)
    return 0


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
        mesh,
        origin=origin,
        normal=normal,
        record_parent_ids=args.record_parent_ids,
        chain=args.chain,
        min_length=args.min_length,
    )

    write(args.outfile, out, file_format=args.output_format)
    return 0
