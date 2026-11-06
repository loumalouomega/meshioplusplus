import argparse
import json

from .._cut import cut_level_set, cut_quadrature
from .._helpers import read, write
from ._common import add_io_args


def _plane(text):
    try:
        values = [float(v) for v in text.split(",")]
    except ValueError:
        values = []
    if len(values) != 6:
        raise argparse.ArgumentTypeError(
            f"a plane is six comma-separated numbers, origin then normal, not {text!r}"
        )
    return values


def _common_args(parser, outfile_help):
    add_io_args(parser, outfile_help, in_help="mesh to cut")
    parser.add_argument(
        "--level-set",
        "-l",
        action="append",
        default=[],
        metavar="NAME",
        help="a scalar point_data array that is a level set (repeatable)",
    )
    parser.add_argument(
        "--plane",
        action="append",
        type=_plane,
        default=[],
        metavar="OX,OY,OZ,NX,NY,NZ",
        help="a plane level set (x - origin) . normal, named plane0, plane1, ... (repeatable)",
    )
    parser.add_argument(
        "--selection",
        "-s",
        type=str,
        default="",
        help="the part, such as 'phi < 0', 'phi = 0' (the interface) or "
        "'a < 0 and b > 0' (default '<first> < 0')",
    )
    parser.add_argument(
        "--cut-only",
        action="store_true",
        help="only the pieces of the cut cells, not the whole cells inside",
    )
    parser.add_argument(
        "--backend",
        choices=("auto", "lut", "quadrays"),
        default="auto",
        help="lut: straight pieces; quadrays: curved, high-order (default auto: "
        "quadrays for a quadratic level set)",
    )
    parser.add_argument(
        "--triangulate",
        action="store_true",
        help="lut: split the pieces into simplices",
    )
    parser.add_argument(
        "--no-parent-ids",
        action="store_true",
        help="do not write the parent cell of each piece or point",
    )
    parser.add_argument("--json", action="store_true", help="print the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="do not print the summary"
    )


def add_args(parser):
    _common_args(parser, "the part to be written to")
    parser.add_argument(
        "--degree",
        type=int,
        default=2,
        help="quadrays: the Lagrange degree of the pieces (default 2)",
    )


def add_quadrature_args(parser):
    _common_args(parser, "the quadrature points to be written to")
    parser.add_argument(
        "--order",
        type=int,
        default=2,
        help="the polynomial degree the rules integrate exactly (default 2)",
    )


def _planes(args):
    return [(tuple(p[:3]), tuple(p[3:])) for p in args.plane]


def cut_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = cut_level_set(
        mesh,
        args.level_set,
        planes=_planes(args),
        selection=args.selection,
        include_uncut=not args.cut_only,
        backend=args.backend,
        degree=args.degree,
        triangulate=args.triangulate,
        record_parent_ids=not args.no_parent_ids,
        return_report=True,
    )
    if args.json:
        print(json.dumps(report))
    elif not args.quiet:
        print(
            f"'{report['selection']}': {report['num_pieces']} cells "
            f"({report['num_cut_cells']} of {report['num_cells']} cut), measure "
            f"{report['measure']:.10g} (dimension {report['dimension']}); {report['backend']}"
        )
    write(args.outfile, out, file_format=args.output_format)


def cut_quadrature_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = cut_quadrature(
        mesh,
        args.level_set,
        planes=_planes(args),
        selection=args.selection,
        order=args.order,
        include_uncut=not args.cut_only,
        backend=args.backend,
        triangulate=args.triangulate,
        record_parent_ids=not args.no_parent_ids,
        return_report=True,
    )
    if args.json:
        print(json.dumps(report))
    elif not args.quiet:
        print(
            f"'{report['selection']}': {report['num_points']} points in "
            f"{report['num_rules']} rules, weights sum {report['sum_weights']:.10g} "
            f"(dimension {report['dimension']}); {report['backend']}"
        )
    write(args.outfile, out, file_format=args.output_format)
