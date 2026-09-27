import sys

from .._helpers import read, reader_map
from .._quality_gate import check_quality, format_quality_gate
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be checked")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--require",
        "-r",
        type=str,
        action="append",
        default=[],
        metavar="SPEC",
        help="thresholds, e.g. 'scaled_jacobian >= 0.2; aspect_ratio <= 5 @ 1%%' "
        "(repeatable)",
    )
    parser.add_argument(
        "--gate",
        type=str,
        default=None,
        metavar="FILE",
        help="read thresholds from FILE, one clause per line (# comments)",
    )
    parser.add_argument(
        "--max-inverted",
        type=int,
        default=0,
        help="the most inverted cells allowed; negative disables (default: 0)",
    )
    parser.add_argument(
        "--max-degenerate",
        type=int,
        default=0,
        help="the most degenerate cells allowed; negative disables (default: 0)",
    )
    parser.add_argument("--json", action="store_true", help="emit the result as JSON")


def check_cmd(args):
    # 0 pass, 1 fail, 2 could not check (unreadable file, bad specification):
    # a CI job must tell "the mesh is bad" from "the check did not run".
    try:
        require = list(args.require)
        if args.gate:
            with open(args.gate, encoding="utf-8") as f:
                require.append(f.read())
        mesh = read(args.infile, file_format=args.input_format)
        result = check_quality(
            mesh,
            require,
            max_inverted=args.max_inverted,
            max_degenerate=args.max_degenerate,
        )
    except Exception as exc:  # noqa: BLE001 -- reported, exit status 2
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.json:
        emit_json(result)
    else:
        print(format_quality_gate(result), end="")
    return 0 if result["passed"] else 1
