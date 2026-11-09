from .._colorize import colorize
from .._helpers import read, write
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(parser, out_help="write the mesh with red/green/blue UInt8 data here")
    parser.add_argument(
        "--field", type=str, required=True, help="point or cell data to colour by"
    )
    parser.add_argument("--cmap", type=str, default="viridis", help="colormap name")
    parser.add_argument("--vmin", type=float, default=None, help="range minimum")
    parser.add_argument("--vmax", type=float, default=None, help="range maximum")
    parser.add_argument(
        "--symmetric",
        action="store_true",
        help="with the automatic range, +-max|v| so zero sits mid-map",
    )
    parser.add_argument(
        "--steps", type=int, default=0, help="bands at fixed steps (0: continuous)"
    )
    parser.add_argument(
        "--nan-color",
        type=str,
        default="#808080",
        help="#rrggbb for a non-finite value",
    )
    parser.add_argument(
        "--component", type=int, default=None, help="colour by this vector component"
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def colorize_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, r = colorize(
        mesh,
        args.field,
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
        symmetric=args.symmetric,
        steps=args.steps,
        nan_color=args.nan_color,
        component=args.component,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(r)
    else:
        print(
            f"{args.field} coloured on the {r['location']}s over "
            f"[{r['vmin']:.6g}, {r['vmax']:.6g}]"
        )
    return 0
