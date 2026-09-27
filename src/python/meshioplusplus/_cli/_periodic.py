import numpy as np

from .._helpers import read, reader_map
from .._periodic import match_periodic_nodes
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file carrying both regions")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    parser.add_argument("--slave", required=True, help="region whose nodes are mapped")
    parser.add_argument("--master", required=True, help="region they map onto")
    parser.add_argument(
        "--translate",
        type=str,
        metavar="DX,DY,DZ",
        help="slave-to-master translation (write --translate=-1,0,0 for a leading minus)",
    )
    parser.add_argument(
        "--rotate",
        type=str,
        metavar="AXIS,DEG",
        help="rotation by DEG degrees about AXIS -- x, y or z, or AX,AY,AZ (so "
        "z,90 or 0,0,1,90) -- applied before --translate",
    )
    parser.add_argument(
        "--origin",
        type=str,
        metavar="X,Y,Z",
        help="point the rotation turns about (default: the origin)",
    )
    parser.add_argument(
        "--matrix",
        type=str,
        metavar="M00,...,M33",
        help="16 comma-separated numbers, a row-major 4x4 affine matrix, overriding "
        "--translate/--rotate",
    )
    parser.add_argument(
        "--atol", type=float, default=1e-8, help="match tolerance (default: 1e-8)"
    )
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="report unmatched slave nodes instead of failing",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="write the pairs as 'slave,master' CSV rows (0-based ids)",
    )
    parser.add_argument("--json", action="store_true", help="emit the pairs as JSON")


def _numbers(text, count, flag):
    if text is None:
        return None
    try:
        values = [float(v) for v in text.split(",")]
    except ValueError:
        raise SystemExit(f"periodic: {flag} expects comma-separated numbers") from None
    if count is not None and len(values) != count:
        raise SystemExit(f"periodic: {flag} expects {count} numbers, got {len(values)}")
    return values


def _rotation(text):
    if text is None:
        return None
    parts = text.split(",")
    if len(parts) == 2 and parts[0].lower() in ("x", "y", "z"):
        return (parts[0].lower(), float(parts[1]))
    values = _numbers(text, 4, "--rotate")
    return (values[:3], values[3])


def periodic_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    pairs, report = match_periodic_nodes(
        mesh,
        args.slave,
        args.master,
        translate=_numbers(args.translate, 3, "--translate"),
        rotate=_rotation(args.rotate),
        matrix=_numbers(args.matrix, 16, "--matrix"),
        origin=_numbers(args.origin, 3, "--origin"),
        atol=args.atol,
        require_complete=not args.allow_incomplete,
        return_report=True,
    )
    if args.output:
        np.savetxt(args.output, pairs, fmt="%d", delimiter=",", header="slave,master")
    if args.json:
        emit_json(
            {
                "num_pairs": len(pairs),
                "pairs": pairs,
                "unmatched": report["unmatched"],
                "num_fixed": report["num_fixed"],
                "max_residual": report["max_residual"],
            }
        )
        return 0
    print(f"periodic pairs {args.slave} -> {args.master}: {len(pairs)}")
    print(f"  unmatched slave nodes: {len(report['unmatched'])}")
    print(f"  fixed points:          {report['num_fixed']}")
    print(f"  max residual:          {report['max_residual']:.3g}")
    if not args.output:
        for s, m in pairs[:20]:
            print(f"  {s} -> {m}")
        if len(pairs) > 20:
            print(f"  ... ({len(pairs) - 20} more; use --output or --json)")
    return 0
