from .._helpers import read, write
from .._voxelize import voxelize
from ._common import add_box_args, add_io_args, parse_bounds, parse_ints


def add_args(parser):
    add_io_args(parser, "grid to be written to")
    parser.add_argument(
        "--resolution",
        type=str,
        default=None,
        help="cell counts 'nx,ny,nz', or 'nx,ny' for 2-D pixels around a planar "
        "polyline (exactly one of --resolution/--cell-size)",
    )
    parser.add_argument(
        "--cell-size",
        type=float,
        default=None,
        help="cubic cell size (exactly one of --resolution/--cell-size)",
    )
    add_box_args(parser)
    parser.add_argument(
        "--fill",
        type=str,
        default="all",
        choices=["all", "surface", "inside"],
        help="which cells to keep: the whole box, only those a triangle passes "
        "through, or only those inside the surface",
    )
    parser.add_argument(
        "--sign",
        type=str,
        default="pseudonormal",
        choices=["pseudonormal", "winding-number"],
        help="how --fill=inside decides what is inside",
    )
    parser.add_argument(
        "--plane",
        type=str,
        default="auto",
        choices=["auto", "xy", "xz", "yz"],
        help="plane of a 2-D skin ('auto': the axis of zero extent)",
    )
    parser.add_argument(
        "--attach-occupancy",
        action="store_true",
        help="attach the Int64 voxel:occupancy cell_data array",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=20000000,
        help="refuse above this many cells (default ~256^3, about 1.5 GB)",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def voxelize_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    resolution = parse_ints(
        args.resolution,
        "voxelize",
        "--resolution",
        counts=(2, 3),
        form="'nx,ny,nz' (or 'nx,ny')",
    )
    bounds = parse_bounds(args.bounds, "voxelize")

    out, report = voxelize(
        mesh,
        resolution=resolution,
        cell_size=args.cell_size,
        bounds=bounds,
        padding=args.padding,
        padding_relative=args.padding_relative,
        fill=args.fill,
        attach_occupancy=args.attach_occupancy,
        max_cells=args.max_cells,
        sign=args.sign,
        watertight_check="warn" if args.fill == "inside" else "off",
        plane=args.plane,
        return_report=True,
    )

    if not args.quiet:
        d = report["dims"]
        s = report["spacing"]
        print(f"voxelized ({args.fill})")
        print(f"  grid:           {d[0]} x {d[1]} x {d[2]}")
        print(f"  cell size:      {s[0]:g}, {s[1]:g}, {s[2]:g}")
        print(f"  cells kept:     {report['num_occupied']}")

    write(args.outfile, out, file_format=args.output_format)
    return 0
