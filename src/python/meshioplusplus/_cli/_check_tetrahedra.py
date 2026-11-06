import sys

from .._check_tetrahedra import check_tetrahedra
from .._helpers import _reader_names, read
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be checked")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=_reader_names(),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--zero-volume-factor",
        type=float,
        default=1e-6,
        help="a volume under this fraction of the longest edge cubed counts as "
        "zero (default: 1e-6)",
    )
    parser.add_argument("--json", action="store_true", help="emit the result as JSON")


def format_tetrahedra(report):
    """The human-readable summary ``check-tetrahedra`` prints."""
    verdict = "valid" if report["valid"] else "NOT valid"
    lines = [
        f"{verdict} ({report['num_tetrahedra']} tetrahedra examined)",
        f"  inverted:    {report['num_inverted']}",
        f"  degenerate:  {report['num_degenerate']}",
        f"  overlaps:    {report['num_overlaps']}",
        f"  min volume:  {report['min_volume']:.6g}",
        f"  max volume:  {report['max_volume']:.6g}",
    ]
    lines += [f"  warning: {w}" for w in report["warnings"]]
    return "\n".join(lines) + "\n"


def check_tetrahedra_cmd(args):
    # 0 valid, 1 not, 2 could not check -- the `check` verb's convention.
    try:
        mesh = read(args.infile, file_format=args.input_format)
        report = check_tetrahedra(mesh, zero_volume_factor=args.zero_volume_factor)
    except Exception as exc:  # noqa: BLE001 -- reported, exit status 2
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.json:
        emit_json(report)
    else:
        print(format_tetrahedra(report), end="")
    return 0 if report["valid"] else 1
