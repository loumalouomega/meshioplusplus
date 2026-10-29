from .._extrude import revolve
from .._helpers import read, write
from ._sweep_common import add_io_args, add_report_args, print_sweep, vec3


def add_args(parser):
    add_io_args(parser, "swept mesh to be written to")
    parser.add_argument(
        "--axis",
        type=str,
        default="0,1,0",
        help="the axis direction 'x,y,z' (negatives need --axis=0,-1,0)",
    )
    parser.add_argument(
        "--axis-point",
        type=str,
        default="0,0,0",
        help="a point of the axis 'x,y,z'",
    )
    parser.add_argument(
        "--angle",
        type=float,
        default=360.0,
        help="degrees, right-hand rule about the axis (default 360, a full turn)",
    )
    parser.add_argument(
        "--num-layers", type=int, default=16, help="steps of the turn (default 16)"
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-9,
        help="on the axis: fraction of the bounding-box diagonal (default 1e-9)",
    )
    parser.add_argument(
        "--record-layers",
        action="store_true",
        help="attach the layer of each cell as revolve:layer cell data",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=20_000_000,
        help="refuse more profile cells times layers (default 20000000)",
    )
    add_report_args(parser)


def revolve_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = revolve(
        mesh,
        vec3(args.axis, "revolve", "axis"),
        axis_point=vec3(args.axis_point, "revolve", "axis-point"),
        angle=args.angle,
        num_layers=args.num_layers,
        tolerance=args.tolerance,
        record_layers=args.record_layers,
        max_cells=args.max_cells,
        return_report=True,
    )
    print_sweep(args, report, revolve=True)
    write(args.outfile, out, file_format=args.output_format)
