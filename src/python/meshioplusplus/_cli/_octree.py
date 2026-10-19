from .._helpers import _writer_map, read, reader_map, write
from .._octree import octree


def add_args(parser):
    parser.add_argument("infile", type=str, help="skin mesh file to be read from")
    parser.add_argument("outfile", type=str, help="tree to be written to")
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
        "--root-cells",
        type=str,
        default="1,1,1",
        help="root cells per axis 'nx,ny,nz' (a forest of cubic roots)",
    )
    parser.add_argument(
        "--min-level", type=int, default=0, help="split every leaf down to this level"
    )
    parser.add_argument(
        "--max-level",
        type=int,
        default=4,
        help="never split beyond this level (at most 20)",
    )
    parser.add_argument(
        "--bounds",
        type=str,
        default=None,
        help="'xlo,ylo,zlo,xhi,yhi,zhi'; the skin's bounding box by default "
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
        "--dimension",
        type=int,
        default=0,
        choices=[0, 2, 3],
        help="2 for a planar skin (quadtree), 3 for a surface (octree), 0 to decide "
        "from the cell types",
    )
    parser.add_argument(
        "--sample-distance",
        type=float,
        default=0.0,
        help="refine further within this distance of the skin",
    )
    parser.add_argument(
        "--balance",
        type=str,
        default="face",
        choices=["none", "face", "full"],
        help="2:1 neighbour balance: across faces, or faces, edges and corners",
    )
    parser.add_argument(
        "--fill",
        type=str,
        default="all",
        choices=["all", "surface", "inside"],
        help="which leaves to keep: all, those the skin passes through, or inside",
    )
    parser.add_argument(
        "--sign",
        type=str,
        default="pseudonormal",
        choices=["pseudonormal", "winding-number"],
        help="how --fill=inside decides what is inside (3-D)",
    )
    parser.add_argument(
        "--attach-distance",
        action="store_true",
        help="attach the signed distance at leaf centres as sdf:distance",
    )
    parser.add_argument(
        "--even-odd",
        action="store_true",
        help="2-D: use the even-odd rule instead of the winding number",
    )
    parser.add_argument(
        "--transition",
        type=str,
        default="hanging",
        choices=["hanging", "polyhedra", "templates"],
        help="output where a fine leaf meets a coarse one: hanging nodes, conforming "
        "polyhedra, or conforming hexahedra/pyramids/tetrahedra",
    )
    parser.add_argument(
        "--snap",
        action="store_true",
        help="move the outer nodes onto the skin without inverting a leaf "
        "(needs --fill surface|inside)",
    )
    parser.add_argument(
        "--surface-level",
        type=int,
        default=None,
        help="split cut leaves down to this level; the criteria below refine further",
    )
    parser.add_argument(
        "--curvature-threshold",
        type=float,
        default=0.0,
        help="refine to --max-level where |mean curvature| exceeds this",
    )
    parser.add_argument(
        "--feature-angle",
        type=float,
        default=0.0,
        help="refine to --max-level along edges sharper than this many degrees",
    )
    parser.add_argument(
        "--size-array",
        type=str,
        default=None,
        help="point_data array of target edge lengths on the skin",
    )
    parser.add_argument(
        "--plane",
        type=str,
        default="auto",
        choices=["auto", "xy", "xz", "yz"],
        help="plane of a 2-D skin",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=20000000,
        help="refuse above this many leaves",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def octree_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    roots = [int(x) for x in args.root_cells.split(",")]
    if len(roots) not in (1, 2, 3):
        raise ValueError("octree: --root-cells expects 'nx,ny,nz'")
    roots = (roots * 3)[:3] if len(roots) == 1 else roots + [1] * (3 - len(roots))
    bounds = None
    if args.bounds is not None:
        bounds = [float(x) for x in args.bounds.split(",")]
        if len(bounds) != 6:
            raise ValueError("octree: --bounds expects 'xlo,ylo,zlo,xhi,yhi,zhi'")

    out, report = octree(
        mesh,
        root_cells=roots,
        min_level=args.min_level,
        max_level=args.max_level,
        bounds=bounds,
        padding=args.padding,
        padding_relative=args.padding_relative,
        dimension=args.dimension,
        sample_distance=args.sample_distance,
        balance=args.balance,
        fill=args.fill,
        attach_distance=args.attach_distance,
        even_odd=args.even_odd,
        max_cells=args.max_cells,
        sign=args.sign,
        transition=args.transition,
        snap=args.snap,
        surface_level=args.surface_level,
        curvature_threshold=args.curvature_threshold,
        feature_angle=args.feature_angle,
        size_array=args.size_array,
        plane=args.plane,
        watertight_check="warn" if args.fill == "inside" else "off",
        return_report=True,
    )

    if not args.quiet:
        kind = "quadtree" if report["dimension"] == 2 else "octree"
        print(f"{kind} ({args.fill})")
        print(f"  root size:      {report['root_size']:g}")
        print(f"  finest level:   {report['finest_level']}")
        print(f"  leaves kept:    {report['num_leaves']}")
        print(f"  hanging nodes:  {report['num_hanging']}")
        if args.snap:
            print(
                f"  snapped nodes:  {report['num_snapped']} "
                f"({report['num_refused']} refused)"
            )

    write(args.outfile, out, file_format=args.output_format)
    return 0
