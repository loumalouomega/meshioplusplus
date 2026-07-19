from .._helpers import _writer_map, read, reader_map, write
from .._reorder import compute_bandwidth, reorder


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument(
        "outfile", type=str, help="renumbered mesh file to be written to"
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
        "--method",
        "-m",
        type=str,
        choices=["rcm", "morton", "hilbert"],
        default="rcm",
        help="renumbering strategy (default: rcm)",
    )
    parser.add_argument(
        "--report",
        "-r",
        action="store_true",
        help="print the connectivity bandwidth before and after",
    )


def reorder_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    before = compute_bandwidth(mesh) if args.report else None
    out = reorder(mesh, method=args.method)
    if args.report:
        after = compute_bandwidth(out)
        print(f"bandwidth ({args.method}): {before} -> {after}")

    write(args.outfile, out, file_format=args.output_format)
    return 0
