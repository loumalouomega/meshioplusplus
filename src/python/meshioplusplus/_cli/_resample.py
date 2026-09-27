from .._sequence import parse_times, resample_sequence, sequence_times
from ._json import emit_json


def add_args(parser):
    parser.add_argument(
        "infile",
        type=str,
        help="the input sequence: a quoted glob ('out_*.vtu') or the first file",
    )
    parser.add_argument(
        "outfile",
        type=str,
        help="a {step}/{index} pattern (one file per target time) or a series file",
    )
    parser.add_argument(
        "--input",
        action="append",
        default=None,
        metavar="FILE",
        help="an EXTRA input file, appended after the positional infile; repeat for each",
    )
    parser.add_argument(
        "--input-format", "-i", type=str, default=None, help="input file format"
    )
    parser.add_argument(
        "--time-from",
        type=str,
        default="auto",
        choices=["auto", "file", "filename", "index"],
        help="where the source step times come from (default: auto)",
    )
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument(
        "--times",
        type=str,
        default=None,
        metavar="SPEC",
        help="the target times: START:STOP:STEP or a comma list T1,T2,...",
    )
    target.add_argument(
        "--times-from",
        type=str,
        default=None,
        metavar="PATTERN",
        help="take the target times from a second sequence (a quoted glob)",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="linear",
        choices=["linear", "nearest", "previous"],
        help="linear blends the bracketing steps; nearest / previous sample and hold",
    )
    parser.add_argument(
        "--clamp",
        action="store_true",
        help="a target outside the source range takes the end step instead of failing",
    )
    parser.add_argument(
        "--blend-points",
        action="store_true",
        help="blend the point coordinates too (a moving mesh)",
    )
    parser.add_argument("--json", action="store_true", help="emit a summary as JSON")


def resample_cmd(args):
    source = [args.infile, *args.input] if args.input else args.infile
    targets = (
        parse_times(args.times)
        if args.times is not None
        else sequence_times(args.times_from, time_from=args.time_from)
    )
    resample_sequence(
        source,
        args.outfile,
        times=targets,
        method=args.method,
        extrapolate="clamp" if args.clamp else "error",
        blend_points=args.blend_points,
        file_format=args.input_format,
        time_from=args.time_from,
    )
    if args.json:
        emit_json({"method": args.method, "num_steps": len(targets), "times": targets})
    else:
        print(f"resample ({args.method}): {len(targets)} step(s) -> {args.outfile}")
    return 0
