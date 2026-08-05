from .._helpers import _writer_map, read, reader_map, write
from .._voxelize import voxelize


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="grid to be written to")
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
        "--resolution",
        type=str,
        default=None,
        help="cell counts 'nx,ny,nz' (exactly one of --resolution/--cell-size)",
    )
    parser.add_argument(
        "--cell-size",
        type=float,
        default=None,
        help="cubic cell size (exactly one of --resolution/--cell-size)",
    )
    parser.add_argument(
        "--bounds",
        type=str,
        default=None,
        help="'xlo,ylo,zlo,xhi,yhi,zhi'; the mesh's bounding box by default "
        "(negatives need --bounds=)",
    )
    parser.add_argument(
        "--padding", type=float, default=0.0, help="grow the box by this on every side"
    )
    parser.add_argument(
        "--padding-relative",
        type=float,
        default=0.0,
        help="grow the box by this fraction of its diagonal on every side",
    )
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

    resolution = None
    if args.resolution is not None:
        resolution = [int(x) for x in args.resolution.split(",")]
        if len(resolution) != 3:
            raise ValueError("voxelize: --resolution expects 'nx,ny,nz'")
    bounds = None
    if args.bounds is not None:
        bounds = [float(x) for x in args.bounds.split(",")]
        if len(bounds) != 6:
            raise ValueError("voxelize: --bounds expects 'xlo,ylo,zlo,xhi,yhi,zhi'")

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
