"""The ``convert-units`` verb: every declared field unit
(``field_data["meshio:unit:<array>"]``) and the points to a consistent unit
system (``doc/units.md``)."""

from .._helpers import read, write
from .._units import convert_units
from ._hull import _io_args
from ._json import emit_json


def add_args(parser):
    _io_args(parser)
    parser.add_argument(
        "--system",
        choices=("SI", "mm-t-s", "mm-kg-ms", "cm-g-s", "mm-g-ms"),
        default="SI",
        help="target consistent unit system (default SI)",
    )
    parser.add_argument(
        "--lenient",
        action="store_true",
        help="leave arrays without a declared unit unchanged (warn) instead of failing",
    )


def convert_units_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    try:
        out, report = convert_units(
            mesh, args.system, strict=not args.lenient, return_report=True
        )
    except ValueError as exc:
        raise SystemExit(f"convert-units: {exc}") from None
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(
            {
                "output": args.outfile,
                "system": args.system,
                "num_arrays_converted": int(report["num_arrays_converted"]),
            }
        )
        return 0
    for w in report.get("warnings", []):
        print(f"warning: {w}")
    print(f"{args.system}: {report['num_arrays_converted']} array(s) converted")
    return 0
