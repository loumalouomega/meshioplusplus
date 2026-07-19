from .._helpers import read, reader_map, write
from .._quality import attach_quality, compute_quality


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="write the metrics into this file as cell_data",
    )


def _format(value: float) -> str:
    return "  n/a " if value != value else f"{value:8.4f}"  # NaN check


def quality(args):
    mesh = read(args.infile, file_format=args.input_format)
    report = compute_quality(mesh)

    print(f"Mesh quality report for {args.infile}")
    print(
        f"  cells: {report['num_cells']}   "
        f"inverted: {report['num_inverted']}   "
        f"degenerate: {report['num_degenerate']}"
    )
    print(f"  {'metric':<24}{'min':>10}{'mean':>10}{'max':>10}{'count':>10}")
    for name, s in report["metrics"].items():
        if s["count"] == 0:
            continue
        print(
            f"  {name:<24}{_format(s['min']):>10}{_format(s['mean']):>10}"
            f"{_format(s['max']):>10}{s['count']:>10}"
        )

    if args.output is not None:
        write(args.output, attach_quality(mesh))

    return 0
