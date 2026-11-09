"""The ``measure`` and ``align-datums`` verbs: elements, dimensions and GD&T
tolerances, and the 3-2-1 datum alignment (``doc/measure.md``)."""

from .._helpers import read, write
from .._measure import align_datums, measure
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="the scanned surface")
    parser.add_argument("spec", type=str, help="the measurement spec (JSON)")
    parser.add_argument(
        "outfile",
        type=str,
        nargs="?",
        default=None,
        help="write the mesh with its element table here",
    )
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def measure_cmd(args):
    out, r = measure(read(args.infile, file_format=args.input_format), args.spec)
    if args.outfile is not None:
        write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(r)
        return 1 if r["num_fail"] > 0 else 0
    for e in r["elements"]:
        p = e["point"]
        print(
            f"element {e['name']:<12} ({p[0]:.9g}, {p[1]:.9g}, {p[2]:.9g}) {e['source']}"
        )
    for c in r["dimensions"] + r["tolerances"]:
        verdict = (
            "" if c["tolerance"] is None else ("  PASS" if c["pass"] else "  FAIL")
        )
        print(f"{c['name']:<16} {c['type']:<16} {c['value']:.9g}{verdict}")
    for w in r["warnings"]:
        print(f"warning: {w}")
    return 1 if r["num_fail"] > 0 else 0


def add_align_args(parser):
    parser.add_argument("infile", type=str, help="a measured mesh (measure)")
    parser.add_argument("outfile", type=str, help="write the aligned mesh here")
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument(
        "--primary", type=str, required=True, help="the element setting Z"
    )
    parser.add_argument(
        "--secondary", type=str, required=True, help="the element setting X"
    )
    parser.add_argument(
        "--origin", type=str, default=None, help="a point-like element or a third plane"
    )
    parser.add_argument("--json", action="store_true", help="emit the matrix as JSON")


def align_datums_cmd(args):
    out, matrix = align_datums(
        read(args.infile, file_format=args.input_format),
        args.primary,
        args.secondary,
        args.origin,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"matrix": matrix.tolist()})
        return 0
    for row in matrix:
        print(" ".join(f"{v:.17g}" for v in row))
    return 0
