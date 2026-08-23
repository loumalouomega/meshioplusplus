from .._helpers import _writer_map, read, reader_map, write
from .._remesh_volume import remesh_volume


def add_args(parser):
    parser.add_argument("infile", type=str, help="volume mesh or closed surface file to be read from")
    parser.add_argument("outfile", type=str, help="retetrahedralized mesh to be written to")
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
        help="cell counts 'nx,ny,nz' of the root lattice (exactly one of "
        "--resolution/--cell-size)",
    )
    parser.add_argument(
        "--cell-size",
        type=float,
        default=None,
        help="cubic cell size of the root lattice (exactly one of "
        "--resolution/--cell-size)",
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
        default=0.1,
        help="grow the box by this fraction of its diagonal on every side (default: 0.1)",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=20000000,
        help="refuse to generate a root lattice above this many cells",
    )
    parser.add_argument(
        "--max-tets",
        type=int,
        default=20000000,
        help="refuse an output with more tets than this (checked after cutting)",
    )
    parser.add_argument(
        "--warp-fraction",
        type=float,
        default=0.35,
        help="fraction of a lattice vertex's shortest incident edge within "
        "which it may be warped onto the surface; 0 disables warping "
        "(default: 0.35, negatives need --warp-fraction=)",
    )
    parser.add_argument(
        "--sign",
        type=str,
        default="pseudonormal",
        choices=["pseudonormal", "winding-number"],
        help="how lattice vertices are classified inside/outside",
    )
    parser.add_argument(
        "--watertight-check",
        type=str,
        default="warn",
        choices=["off", "warn", "error"],
        help="what to do about an input surface that is not watertight",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def remesh_volume_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    resolution = None
    if args.resolution is not None:
        resolution = [int(x) for x in args.resolution.split(",")]
        if len(resolution) != 3:
            raise ValueError("remesh-volume: --resolution expects 'nx,ny,nz'")
    bounds = None
    if args.bounds is not None:
        bounds = [float(x) for x in args.bounds.split(",")]
        if len(bounds) != 6:
            raise ValueError("remesh-volume: --bounds expects 'xlo,ylo,zlo,xhi,yhi,zhi'")

    out, report = remesh_volume(
        mesh,
        resolution=resolution,
        cell_size=args.cell_size,
        bounds=bounds,
        padding=args.padding,
        padding_relative=args.padding_relative,
        max_cells=args.max_cells,
        max_tets=args.max_tets,
        warp_fraction=args.warp_fraction,
        sign=args.sign,
        watertight_check=args.watertight_check,
        return_report=True,
    )

    if not args.quiet:
        print("retetrahedralized")
        print(f"  tets emitted:             {report['num_tets']}")
        print(f"  vertices warped:          {report['num_vertices_warped']}")
        print(f"  tets rejected (degen.):   {report['num_tets_rejected']}")
        print(f"  non-manifold edges:       {report['num_non_manifold_edges']}")

    write(args.outfile, out, file_format=args.output_format)
    return 0
