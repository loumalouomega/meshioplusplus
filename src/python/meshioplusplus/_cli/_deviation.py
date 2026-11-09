from .._deviation import deviation
from .._helpers import _reader_names, _writer_names, read, write
from ._common import split_csv
from ._json import emit_json


def add_args(parser):
    parser.add_argument(
        "scan", type=str, help="scan (point cloud or surface) to measure"
    )
    parser.add_argument(
        "reference",
        type=str,
        help="reference surface, or with --patch a mesh fit-primitives wrote",
    )
    parser.add_argument(
        "outfile",
        type=str,
        nargs="?",
        default=None,
        help="optional: write the scan with deviation:distance here",
    )
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=_reader_names(),
        help="input file format (the scan)",
        default=None,
    )
    parser.add_argument(
        "--reference-format",
        type=str,
        choices=_reader_names(),
        help="input file format of the reference",
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
        "--patch",
        type=int,
        default=None,
        metavar="ID",
        help="measure against the reference's fitted primitive of this patch label",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.0,
        metavar="T",
        help="report the share of points within +-T",
    )
    parser.add_argument(
        "--min-share",
        type=float,
        default=None,
        metavar="S",
        help="exit with status 1 when the share within --tolerance is below S",
    )
    parser.add_argument(
        "--bins", type=int, default=20, help="histogram bins (default: 20)"
    )
    parser.add_argument(
        "--range",
        type=str,
        default=None,
        metavar="LO,HI",
        help="histogram range (default: +-max |d|)",
    )
    parser.add_argument(
        "--max-normal-angle",
        type=float,
        default=90.0,
        metavar="DEG",
        help="far-side rule: skip reference triangles facing away by more (default: 90; "
        "180 turns it off)",
    )
    parser.add_argument(
        "--reference-region",
        type=str,
        default="",
        help="named cell region of the reference",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def deviation_cmd(args):
    scan = read(args.scan, file_format=args.input_format)
    reference = read(args.reference, file_format=args.reference_format)
    rng = None
    if args.range is not None:
        parts = split_csv(args.range)
        if len(parts) != 2:
            raise SystemExit("meshio++: deviation: --range takes LO,HI")
        rng = (float(parts[0]), float(parts[1]))
    out, r = deviation(
        scan,
        reference,
        patch=args.patch,
        tolerance=args.tolerance,
        bins=args.bins,
        range=rng,
        max_normal_angle=args.max_normal_angle,
        reference_region=args.reference_region,
        attach=args.outfile is not None,
        return_report=True,
    )
    if args.outfile is not None:
        write(args.outfile, out, file_format=args.output_format)
    failed = (
        args.min_share is not None
        and args.tolerance > 0.0
        and not r["share_within"] >= args.min_share
    )
    if args.json:
        obj = dict(r)
        if args.min_share is not None:
            obj["min_share"] = args.min_share
            obj["passed"] = not failed
        emit_json(obj)
    else:
        print(f"{r['num_measured']} of {r['num_points']} points measured")
        if r["num_measured"] > 0:
            p = r["percentiles"]
            print(
                f"  min {r['min']:.9g}, max {r['max']:.9g}, mean {r['mean']:.9g}, "
                f"rms {r['rms']:.9g}, sigma {r['sigma']:.9g}"
            )
            print(
                f"  p1 {p['p1']:.6g}, p5 {p['p5']:.6g}, p50 {p['p50']:.6g}, "
                f"p95 {p['p95']:.6g}, p99 {p['p99']:.6g}; suggested range "
                f"+-{r['suggested_range']:.6g}"
            )
        if args.tolerance > 0.0:
            print(f"  {100.0 * r['share_within']:.2f}% within +-{args.tolerance:g}")
        if args.min_share is not None:
            print(f"  {'FAIL' if failed else 'pass'}: share >= {args.min_share:g}")
        for w in r["warnings"]:
            print(f"  warning: {w}")
    return 1 if failed else 0
