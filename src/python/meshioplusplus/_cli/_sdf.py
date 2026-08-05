from .._helpers import _writer_map, read, reader_map, write
from .._sdf import compute_sdf


def add_args(parser):
    parser.add_argument("infile", type=str, help="surface mesh to be read from")
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
        "--structure",
        type=str,
        default="voxel",
        choices=["voxel", "octree"],
        help="a dense lattice, or one refined near the surface. An octree's "
        "output is 1-irregular (it has hanging nodes)",
    )
    parser.add_argument(
        "--resolution",
        type=str,
        default=None,
        help="cell counts 'nx,ny,nz' (voxel only; exactly one of "
        "--resolution/--cell-size)",
    )
    parser.add_argument(
        "--cell-size",
        type=float,
        default=None,
        help="cubic cell size (voxel only)",
    )
    parser.add_argument(
        "--bounds",
        type=str,
        default=None,
        help="'xlo,ylo,zlo,xhi,yhi,zhi'; the surface's bounding box by default "
        "(negatives need --bounds=)",
    )
    parser.add_argument(
        "--padding", type=float, default=0.0, help="grow the box by this on every side"
    )
    parser.add_argument(
        "--padding-relative",
        type=float,
        default=0.1,
        help="grow the box by this fraction of its diagonal on every side; the "
        "default is non-zero because a field that stops at the surface is not "
        "much use",
    )
    parser.add_argument(
        "--root-resolution",
        type=int,
        default=8,
        help="octree: cell count per axis of the root lattice",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=4,
        help="octree: how many refinement passes",
    )
    parser.add_argument(
        "--band-cells",
        type=float,
        default=1.0,
        help="octree: refine while |distance| <= this * the cell's own diagonal",
    )
    parser.add_argument(
        "--sign",
        type=str,
        default="pseudonormal",
        choices=["unsigned", "pseudonormal", "winding-number"],
        help="how to decide which side of the surface a point is on",
    )
    parser.add_argument(
        "--location",
        type=str,
        default="corner",
        choices=["corner", "center"],
        help="evaluate at the grid's points (point_data) or its cell centres "
        "(cell_data). isosurface needs point_data",
    )
    parser.add_argument(
        "--band",
        type=float,
        default=0.0,
        help="clamp distances beyond this and mark them in sdf:band",
    )
    parser.add_argument(
        "--watertight-check",
        type=str,
        default="warn",
        choices=["off", "warn", "error"],
        help="what to do about a surface that is not closed",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=20000000,
        help="refuse above this many cells, re-checked after every octree pass",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def sdf_cmd(args):
    surface = read(args.infile, file_format=args.input_format)

    resolution = None
    if args.resolution is not None:
        resolution = [int(x) for x in args.resolution.split(",")]
        if len(resolution) != 3:
            raise ValueError("sdf: --resolution expects 'nx,ny,nz'")
    bounds = None
    if args.bounds is not None:
        bounds = [float(x) for x in args.bounds.split(",")]
        if len(bounds) != 6:
            raise ValueError("sdf: --bounds expects 'xlo,ylo,zlo,xhi,yhi,zhi'")

    out, report = compute_sdf(
        surface,
        structure=args.structure,
        resolution=resolution,
        cell_size=args.cell_size,
        bounds=bounds,
        padding=args.padding,
        padding_relative=args.padding_relative,
        root_resolution=args.root_resolution,
        max_depth=args.max_depth,
        band_cells=args.band_cells,
        max_cells=args.max_cells,
        sign=args.sign,
        location=args.location,
        band=args.band,
        watertight_check=args.watertight_check,
        return_report=True,
    )

    if not args.quiet:
        d = report["dims"]
        s = report["spacing"]
        q = report["quality"]
        cells = sum(len(cb.data) for cb in out.cells)
        print(f"signed distance field ({args.structure})")
        print(f"  root grid:      {d[0]} x {d[1]} x {d[2]}")
        print(f"  finest cell:    {s[0]:g}, {s[1]:g}, {s[2]:g}")
        if report["max_depth"]:
            print(f"  octree depth:   {report['max_depth']}")
        print(f"  cells:          {cells}")
        if args.band > 0.0:
            print(f"  banded:         {report['num_banded']}")
        if not q["watertight"]:
            print(
                f"  surface:        NOT watertight "
                f"({q['boundary_edges']} boundary, "
                f"{q['non_manifold_edges']} non-manifold, "
                f"{q['inconsistent_pairs']} inconsistent, "
                f"{q['degenerate_triangles']} degenerate)"
            )

    write(args.outfile, out, file_format=args.output_format)
    return 0
