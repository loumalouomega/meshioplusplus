from .._helpers import read, write
from .._remesh_volume import remesh_volume
from ._common import add_box_args, add_io_args, parse_bounds, parse_ints


def add_args(parser):
    add_io_args(
        parser,
        "retetrahedralized mesh to be written to",
        in_help="volume mesh or closed surface file to be read from",
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
    add_box_args(
        parser,
        padding_relative=0.1,
        padding_relative_help=(
            "grow the box by this fraction of its diagonal on every side (default: 0.1)"
        ),
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

    resolution = parse_ints(args.resolution, "remesh-volume", "--resolution")
    bounds = parse_bounds(args.bounds, "remesh-volume")

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
