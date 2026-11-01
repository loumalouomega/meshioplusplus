from .._helpers import _writer_map, read, reader_map, write
from .._select import select
from ._crop import _parse_where
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="mesh with the region, written to")
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
        "--name", "-n", type=str, required=True, help="the region's name"
    )
    parser.add_argument(
        "--kind",
        type=str,
        choices=["cell", "point", "side"],
        default="cell",
        help="what the region groups (default: cell)",
    )
    parser.add_argument("--tag", type=int, default=-1, help="the region's integer tag")
    parser.add_argument("--box", type=str, help="'xmin,ymin,zmin,xmax,ymax,zmax'")
    parser.add_argument(
        "--plane", type=str, help="'px,py,pz,nx,ny,nz' (point + normal)"
    )
    parser.add_argument("--sphere", type=str, help="'cx,cy,cz,r'")
    parser.add_argument(
        "--where",
        type=str,
        help="a data predicate 'NAME OP VALUE' on point_data (points) or "
        "cell_data (cells; a side takes its cell's value)",
    )
    parser.add_argument(
        "--quality",
        type=str,
        help="a quality predicate 'METRIC OP VALUE', e.g. 'scaled_jacobian < 0.3'",
    )
    parser.add_argument(
        "--component", type=int, default=None, help="the k-th node-connected component"
    )
    parser.add_argument(
        "--seed",
        type=str,
        help="'x,y,z': grow a smooth patch from the facet nearest it",
    )
    parser.add_argument(
        "--feature-angle",
        type=float,
        default=30.0,
        help="the largest dihedral angle, in degrees, the growth crosses (default: 30)",
    )
    parser.add_argument(
        "--normal-tolerance",
        type=float,
        default=0.0,
        help="the largest angle to the seed's normal, in degrees; 0 disables (default)",
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["all", "any", "centroid"],
        default="all",
        help="how a cell or side reduces its nodes for --box/--plane/--sphere",
    )
    parser.add_argument(
        "--all-facets",
        action="store_true",
        help="side regions: consider interior facets too, not only boundary ones",
    )
    parser.add_argument("--json", action="store_true", help="emit the counts as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def _floats(text, n, what):
    vals = [float(x) for x in text.split(",")]
    if len(vals) != n:
        raise ValueError(f"select: --{what} expects {n} comma-separated numbers")
    return vals


def select_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    kwargs = {}
    if args.box is not None:
        kwargs["box"] = _floats(args.box, 6, "box")
    if args.plane is not None:
        vals = _floats(args.plane, 6, "plane")
        kwargs["plane"] = (vals[:3], vals[3:])
    if args.sphere is not None:
        vals = _floats(args.sphere, 4, "sphere")
        kwargs["sphere"] = (vals[:3], vals[3])
    if args.where is not None:
        kwargs["where"] = _parse_where(args.where, "select")
    if args.quality is not None:
        kwargs["quality"] = _parse_where(args.quality, "select")
    if args.seed is not None:
        kwargs["seed"] = _floats(args.seed, 3, "seed")
    out, report = select(
        mesh,
        args.name,
        args.kind,
        component=args.component,
        feature_angle=args.feature_angle,
        normal_tolerance=args.normal_tolerance,
        mode=args.mode,
        boundary_only=not args.all_facets,
        tag=args.tag,
        return_report=True,
        **kwargs,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(report)
    elif not args.quiet:
        print(
            f"select: {report['num_selected']} of {report['num_candidates']} "
            f"{args.kind}(s) in region '{args.name}'"
        )
        if report["num_grown"] >= 0:
            print(f"  grown patch:  {report['num_grown']} facet(s)")
        if report["num_components"] >= 0:
            print(f"  components:   {report['num_components']}")
    return 0
