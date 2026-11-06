"""The ``hull`` and ``delaunay`` verbs: the convex hull of a mesh's points and
their Delaunay triangulation in the ``x, y`` plane (``doc/hull.md``)."""

from .._helpers import read, write
from .._hull import convex_hull, delaunay_2d
from ._common import add_io_args
from ._json import emit_json


def _io_args(parser):
    add_io_args(
        parser, "output mesh", in_help="mesh or point-cloud file", output_format=False
    )
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def add_hull_args(parser):
    _io_args(parser)
    parser.add_argument(
        "--dimension",
        choices=("0", "2", "3"),
        default="0",
        help="2 (the x, y footprint), 3, or 0 to decide from the points",
    )


def add_delaunay_args(parser):
    _io_args(parser)


def hull_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = convex_hull(mesh, dimension=int(args.dimension), return_report=True)
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return 0
    kind = "triangles" if report["dimension"] == 3 else "edges"
    print(
        f"{report['dimension']}-D hull: {report['num_vertices']} vertices, "
        f"{report['num_facets']} {kind}"
    )
    return 0


def delaunay_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = delaunay_2d(mesh, return_report=True)
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return 0
    text = f"Delaunay triangulation: {report['num_triangles']} triangles"
    if report["num_duplicates"]:
        text += f", {report['num_duplicates']} duplicate point(s) left out"
    print(text)
    return 0
