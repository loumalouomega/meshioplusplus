from .._helpers import _writer_map, read, reader_map, write
from .._split import split


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument(
        "outpattern",
        type=str,
        help="output filename pattern containing '{key}' (e.g. out_{key}.vtu)",
    )
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
        "--by",
        type=str,
        choices=["type", "region", "regions", "component"],
        default="type",
        help=(
            "split criterion (default: type). 'regions' (plural) is one piece "
            "per named Cell region and is NOT a partition -- a cell in "
            "several regions lands in several pieces. See 'meshioplusplus "
            "regions' to list a mesh's regions first."
        ),
    )
    parser.add_argument(
        "--tag",
        type=str,
        default=None,
        help="for --by region: the integer cell_data name to split on",
    )


def split_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    pieces = split(mesh, by=args.by, tag=args.tag)

    if "{key}" not in args.outpattern:
        raise ValueError(
            "split: output pattern must contain '{key}' (e.g. out_{key}.vtu)"
        )

    print(f"split into {len(pieces)} piece(s) by {args.by}")
    for key, piece in pieces.items():
        path = args.outpattern.format(key=key)
        ncells = sum(len(cb.data) for cb in piece.cells)
        print(f"  {key}: {len(piece.points)} points, {ncells} cells -> {path}")
        write(path, piece, file_format=args.output_format)
    return 0
