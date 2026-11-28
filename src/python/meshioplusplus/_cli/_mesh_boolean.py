"""The ``mesh-boolean`` verb: union, difference or intersection of two closed,
oriented triangle surfaces through the registered boolean provider
(``doc/mesh_boolean.md``)."""

from .._helpers import _reader_names, _writer_names, read, write
from .._mesh_boolean import mesh_boolean
from ._json import emit_json


def add_args(parser):
    parser.add_argument("a", type=str, help="first closed surface")
    parser.add_argument("b", type=str, help="second closed surface")
    parser.add_argument("outfile", type=str, help="result surface to be written to")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=_reader_names(),
        help="input file format (both inputs)",
        default=None,
    )
    parser.add_argument(
        "--output-format",
        "-o",
        type=str,
        choices=_writer_names(),
        help="output file format",
        default=None,
    )
    parser.add_argument(
        "--op",
        choices=("union", "difference", "intersection"),
        default="union",
        help="the boolean (default union)",
    )
    parser.add_argument(
        "--backend", default="", help="provider name (default: the registered one)"
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def mesh_boolean_cmd(args):
    a = read(args.a, file_format=args.input_format)
    b = read(args.b, file_format=args.input_format)
    try:
        out, report = mesh_boolean(
            a, b, args.op, backend=args.backend, return_report=True
        )
    except RuntimeError as exc:
        raise SystemExit(f"mesh-boolean: {exc}") from None
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return 0
    print(f"{report['backend']}: {report['num_triangles']} triangle(s)")
    return 0
