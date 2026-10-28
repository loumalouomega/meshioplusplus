from .._extrude import extrude
from .._helpers import read, write
from ._sweep_common import add_io_args, add_report_args, print_sweep, vec3


def add_args(parser):
    add_io_args(parser, "swept mesh to be written to")
    parser.add_argument(
        "--direction",
        type=str,
        default="0,0,1",
        help="the whole displacement 'x,y,z' (negatives need --direction=0,0,-1)",
    )
    parser.add_argument(
        "--num-layers", type=int, default=1, help="how many layers (default 1)"
    )
    parser.add_argument(
        "--ratio",
        type=float,
        default=1.0,
        help="each layer's thickness over the one below's (default 1, uniform)",
    )
    parser.add_argument(
        "--offsets",
        type=str,
        default=None,
        help="explicit level fractions 'a,b,...', increasing and ending at 1",
    )
    parser.add_argument(
        "--record-layers",
        action="store_true",
        help="attach the layer of each cell as extrude:layer cell data",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=20_000_000,
        help="refuse more profile cells times layers (default 20000000)",
    )
    add_report_args(parser)


def extrude_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    offsets = None
    if args.offsets is not None:
        offsets = [float(x) for x in args.offsets.split(",")]
    out, report = extrude(
        mesh,
        vec3(args.direction, "extrude", "direction"),
        num_layers=args.num_layers,
        ratio=args.ratio,
        offsets=offsets,
        record_layers=args.record_layers,
        max_cells=args.max_cells,
        return_report=True,
    )
    print_sweep(args, report, revolve=False)
    write(args.outfile, out, file_format=args.output_format)
