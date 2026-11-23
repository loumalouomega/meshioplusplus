"""The ``sequence-stats`` and ``sequence-derivative`` verbs: statistics and time
derivatives over the steps of a sequence (``doc/sequence_stats.md``)."""

from .._helpers import write
from .._sequence import write_sequence
from .._sequence_stats import STATISTICS, sequence_derivative, sequence_statistics
from ._common import split_csv
from ._json import emit_json


def _common(parser):
    parser.add_argument(
        "infile",
        type=str,
        help="the input sequence: a quoted glob ('out_*.vtu') or the first file",
    )
    parser.add_argument("outfile", type=str)
    parser.add_argument(
        "--input-format", "-i", type=str, default=None, help="input file format"
    )
    parser.add_argument(
        "--time-from",
        type=str,
        default="auto",
        choices=["auto", "file", "filename", "index"],
        help="where the step times come from (default: auto)",
    )
    parser.add_argument(
        "--location",
        type=str,
        default="point",
        choices=["point", "cell"],
        help="point_data or cell_data (default: point)",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def add_stats_args(parser):
    _common(parser)
    parser.add_argument(
        "--arrays",
        type=str,
        default=None,
        metavar="NAME[,NAME...]",
        help="the arrays to reduce (default: every numeric array of the first step)",
    )
    parser.add_argument(
        "--stats",
        type=str,
        default="mean,std",
        metavar=",".join(STATISTICS),
        help="comma-separated statistics (default: mean,std)",
    )
    parser.add_argument(
        "--ddof",
        type=int,
        default=0,
        help="delta degrees of freedom of std: 0 population, 1 sample (default 0)",
    )
    parser.add_argument("--output-format", "-o", type=str, default=None)


def sequence_stats_cmd(args):
    out, report = sequence_statistics(
        args.infile,
        split_csv(args.arrays) if args.arrays else None,
        tuple(split_csv(args.stats)),
        location=args.location,
        ddof=args.ddof,
        file_format=args.input_format,
        time_from=args.time_from,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return 0
    print(
        f"{', '.join(split_csv(args.stats))} of {len(report['arrays'])} array(s) over "
        f"{report['num_steps']} step(s), t = {report['first_time']} to {report['last_time']}"
    )
    return 0


def add_derivative_args(parser):
    _common(parser)
    parser.add_argument(
        "--array", type=str, required=True, help="the array to differentiate in time"
    )
    parser.add_argument(
        "--scheme",
        type=str,
        default="central",
        choices=["central", "forward"],
        help="central (second order inside, one-sided at the ends) or forward",
    )
    parser.add_argument(
        "--output-name",
        type=str,
        default=None,
        help="the derivative array's name (default: <array>:dt)",
    )


def sequence_derivative_cmd(args):
    steps = sequence_derivative(
        args.infile,
        args.array,
        location=args.location,
        scheme=args.scheme,
        output=args.output_name,
        file_format=args.input_format,
        time_from=args.time_from,
    )
    written = write_sequence(args.outfile, steps)
    if args.json:
        emit_json({"output": args.outfile, "num_steps": len(written), "files": written})
        return 0
    print(f"wrote the time derivative of '{args.array}' for {len(written)} step(s)")
    return 0
