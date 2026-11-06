"""The ``coordinates`` verb: positions, vector and tensor fields in cylindrical
or spherical components about an axis (``doc/coordinates.md``)."""

from .._coordinates import convert_coordinates
from .._helpers import read, write
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(parser, "output mesh", in_help="mesh file", output_format=False)
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument(
        "--system", choices=("cylindrical", "spherical"), default="cylindrical"
    )
    parser.add_argument("--axis", type=str, default="0,0,1", metavar="X,Y,Z")
    parser.add_argument("--origin", type=str, default="0,0,0", metavar="X,Y,Z")
    parser.add_argument("--reference", type=str, default="1,0,0", metavar="X,Y,Z")
    parser.add_argument(
        "--arrays",
        type=str,
        default="",
        metavar="A,B",
        help="vector or tensor point/cell data to rewrite in local components",
    )
    parser.add_argument(
        "--inverse",
        action="store_true",
        help="read the arrays as local components and rewrite them as Cartesian",
    )
    parser.add_argument(
        "--no-record", action="store_true", help="do not attach coords:<system>"
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def _vec3(text, flag):
    try:
        values = [float(v) for v in text.split(",")]
    except ValueError:
        raise SystemExit(f"coordinates: {flag} expects three numbers") from None
    if len(values) != 3:
        raise SystemExit(f"coordinates: {flag} expects three numbers")
    return values


def coordinates_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = convert_coordinates(
        mesh,
        args.system,
        axis=_vec3(args.axis, "--axis"),
        origin=_vec3(args.origin, "--origin"),
        reference=_vec3(args.reference, "--reference"),
        arrays=[a for a in args.arrays.split(",") if a],
        inverse=args.inverse,
        record_coordinates=not args.no_record,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, "system": args.system, **report})
        return 0
    print(
        f"{args.system} components: {report['num_arrays']} array(s) rewritten, "
        f"{report['num_singular']} point(s) on the axis"
    )
    return 0
