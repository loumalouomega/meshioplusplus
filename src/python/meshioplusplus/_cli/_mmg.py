"""The MMG backend's subcommands (roadmap 3.2): ``adapt``,
``discretize-levelset``, ``mesh-2d`` and ``move``."""

from .._helpers import _writer_map, read, reader_map, write
from .._mmg import adapt, discretize_levelset, mesh_domain, move


def _add_io_args(parser, what):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help=f"{what} mesh to be written to")
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


def add_mmg_args(parser):
    """The parameters every MMG operation shares."""
    for name, text in (
        ("hmin", "minimal edge size"),
        ("hmax", "maximal edge size"),
        ("hausd", "how far the boundary may move (Hausdorff distance)"),
        ("hgrad", "gradation"),
        ("hgradreq", "gradation from required entities"),
    ):
        parser.add_argument(f"--{name}", type=float, default=None, help=text)
    parser.add_argument(
        "--angle-detection",
        type=float,
        default=45.0,
        help="ridge angle in degrees; 0 turns the detection off (default: 45)",
    )
    for name, text in (
        ("optim", "keep the size, only improve the quality"),
        ("noinsert", "no point insertion or removal"),
        ("noswap", "no edge or face flip"),
        ("nomove", "no point relocation"),
        ("nosurf", "leave the boundary untouched (not mmgs)"),
        ("nreg", "regularise the normals"),
        ("opnbdy", "keep internal open boundaries (mmg3d, mmg2d)"),
    ):
        parser.add_argument(f"--{name}", action="store_true", help=text)
    parser.add_argument(
        "--frozen", type=str, default=None, help="a point region MMG must keep"
    )
    parser.add_argument(
        "--required",
        type=str,
        action="append",
        default=None,
        help="a region whose entities MMG must keep (repeatable)",
    )
    parser.add_argument(
        "--verbose-mmg", type=int, default=-1, help="MMG's verbosity (default: -1)"
    )
    parser.add_argument("--quiet", action="store_true", help="no run summary")


def mmg_kwargs(args):
    return dict(
        hmin=args.hmin,
        hmax=args.hmax,
        hausd=args.hausd,
        hgrad=args.hgrad,
        hgradreq=args.hgradreq,
        angle_detection=args.angle_detection if args.angle_detection > 0 else None,
        optim=args.optim,
        noinsert=args.noinsert,
        noswap=args.noswap,
        nomove=args.nomove,
        nosurf=args.nosurf,
        nreg=args.nreg,
        opnbdy=args.opnbdy,
        frozen=args.frozen,
        required=args.required,
        verbose=args.verbose_mmg,
    )


def summary(what, report):
    q = report["quality"]
    text = (
        f"{what} with {report['library']} ({report['status']}): "
        f"{report['num_points_in']} -> {report['num_points_out']} points, "
        f"{report['num_elements_in']} -> {report['num_elements_out']} elements, "
        f"aspect ratio min {q['min']:.3g} mean {q['mean']:.3g}"
    )
    for warning in report["warnings"]:
        text += f"\n  warning: {warning}"
    return text


def add_adapt_args(parser):
    _add_io_args(parser, "adapted")
    parser.add_argument("--hsiz", type=float, default=None, help="constant edge size")
    parser.add_argument(
        "--size", type=str, default=None, help="point_data name of a size field"
    )
    parser.add_argument(
        "--metric",
        type=str,
        default=None,
        help="point_data name of a metric (lower triangle, Medit .sol order)",
    )
    add_mmg_args(parser)


def adapt_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = adapt(
        mesh,
        hsiz=args.hsiz,
        size=args.size,
        metric=args.metric,
        return_report=True,
        **mmg_kwargs(args),
    )
    if not args.quiet:
        print(summary("adapted", report))
    write(args.outfile, out, file_format=args.output_format)
    return 0


def add_discretize_levelset_args(parser):
    _add_io_args(parser, "cut")
    parser.add_argument(
        "--level-set",
        type=str,
        required=True,
        help="point_data name of the level set (a signed distance, say)",
    )
    parser.add_argument(
        "--isovalue", type=float, default=0.0, help="the value to cut at (default: 0)"
    )
    parser.add_argument(
        "--rmc",
        type=float,
        default=None,
        help="remove cut components smaller than this volume fraction",
    )
    parser.add_argument(
        "--material",
        type=str,
        action="append",
        default=None,
        metavar="REF,SPLIT,IN,OUT",
        help="a material: the cells of medit:ref REF cut (SPLIT 1) into IN and "
        "OUT, or kept whole (SPLIT 0); repeatable",
    )
    add_mmg_args(parser)


def _material(text):
    parts = text.split(",")
    if len(parts) != 4:
        raise SystemExit(f"--material expects REF,SPLIT,IN,OUT, got {text!r}")
    ref, split, rin, rout = (int(p) for p in parts)
    return ref, bool(split), rin, rout


def discretize_levelset_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = discretize_levelset(
        mesh,
        args.level_set,
        args.isovalue,
        materials=[_material(m) for m in args.material or []],
        remove_small_components=args.rmc,
        return_report=True,
        **mmg_kwargs(args),
    )
    if not args.quiet:
        print(summary("cut", report))
    write(args.outfile, out, file_format=args.output_format)
    return 0


def add_mesh_2d_args(parser):
    _add_io_args(parser, "triangulated")
    parser.add_argument("--hsiz", type=float, default=None, help="constant edge size")
    parser.add_argument(
        "--size", type=str, default=None, help="point_data name of a size field"
    )
    parser.add_argument(
        "--hole",
        type=str,
        action="append",
        default=None,
        metavar="X,Y",
        help="a point inside a region to leave out (repeatable)",
    )
    add_mmg_args(parser)


def _hole(text):
    parts = text.split(",")
    if len(parts) != 2:
        raise SystemExit(f"--hole expects X,Y, got {text!r}")
    return float(parts[0]), float(parts[1])


def mesh_2d_cmd(args):
    boundary = read(args.infile, file_format=args.input_format)
    out, report = mesh_domain(
        boundary,
        hsiz=args.hsiz,
        size=args.size,
        holes=[_hole(h) for h in args.hole or []],
        return_report=True,
        **mmg_kwargs(args),
    )
    if not args.quiet:
        print(
            summary("triangulated", report)
            + f"\n  {report['num_subdomains']} subdomain(s), "
            f"{report['num_holes']} hole(s)"
        )
    write(args.outfile, out, file_format=args.output_format)
    return 0


def add_move_args(parser):
    _add_io_args(parser, "moved")
    parser.add_argument(
        "--displacement",
        type=str,
        required=True,
        help="point_data name of the displacement vector",
    )
    parser.add_argument(
        "--moving",
        type=str,
        action="append",
        required=True,
        help="a cell or side region whose boundary moves (repeatable)",
    )
    parser.add_argument(
        "--lagrangian",
        type=int,
        choices=[0, 1, 2],
        default=1,
        help="MMG's lag: 0 moves only, 1 also swaps (default), 2 also inserts",
    )
    add_mmg_args(parser)


def move_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = move(
        mesh,
        args.displacement,
        args.moving,
        lagrangian=args.lagrangian,
        return_report=True,
        **mmg_kwargs(args),
    )
    if not args.quiet:
        print(summary("moved", report))
    write(args.outfile, out, file_format=args.output_format)
    return 0
