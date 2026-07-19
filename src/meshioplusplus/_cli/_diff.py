from .._diff import _format_report, diff
from .._helpers import read, reader_map


def add_args(parser):
    parser.add_argument("infile_a", type=str, help="first mesh file")
    parser.add_argument("infile_b", type=str, help="second mesh file")
    parser.add_argument(
        "--input-format-a",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input format of the first file",
        default=None,
    )
    parser.add_argument(
        "--input-format-b",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input format of the second file",
        default=None,
    )
    parser.add_argument(
        "--atol",
        type=float,
        default=1e-12,
        help="absolute tolerance in abs_err <= atol + rtol*|expected| (default: 1e-12)",
    )
    parser.add_argument(
        "--rtol",
        type=float,
        default=1e-9,
        help="relative tolerance in abs_err <= atol + rtol*|expected| (default: 1e-9)",
    )
    parser.add_argument(
        "--unordered",
        action="store_true",
        help="match points by spatial proximity (tolerant to a shuffled node order)",
    )
    parser.add_argument(
        "--exact",
        action="store_true",
        help="only a bitwise-identical result passes (tolerated drift exits nonzero)",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="print nothing; communicate equality only via the exit code",
    )


def diff_cmd(args):
    mesh_a = read(args.infile_a, file_format=args.input_format_a)
    mesh_b = read(args.infile_b, file_format=args.input_format_b)

    report = diff(
        mesh_a,
        mesh_b,
        atol=args.atol,
        rtol=args.rtol,
        unordered=args.unordered,
    )

    if not args.quiet:
        print(_format_report(report))

    verdict = report["verdict"]
    equal = verdict == "identical" or (
        verdict == "equal within tolerance" and not args.exact
    )
    return 0 if equal else 1
