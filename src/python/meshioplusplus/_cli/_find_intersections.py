import sys

from .._find_intersections import find_intersections
from .._helpers import _writer_map, read, reader_map, write
from ._json import emit_json

_JSON_KEYS = (
    "num_pairs",
    "num_triangle_pairs",
    "num_coplanar",
    "num_segments",
    "num_cells_a",
    "num_cells_b",
    "truncated",
    "clean",
)


def add_args(parser):
    parser.add_argument("infile", type=str, help="surface mesh file to be checked")
    parser.add_argument(
        "outfile",
        type=str,
        nargs="?",
        default=None,
        help="optional: write the mesh with the intersections:count cell data here",
    )
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--output-format",
        "-o",
        type=str,
        choices=sorted(list(_writer_map.keys())),
        help="output file format",
        default=None,
    )
    parser.add_argument(
        "--other",
        type=str,
        default=None,
        metavar="FILE",
        help="a second surface: report where the two meet instead of self-intersections",
    )
    parser.add_argument(
        "--other-format",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input format of --other",
        default=None,
    )
    parser.add_argument(
        "--segments",
        type=str,
        default=None,
        metavar="FILE",
        help="write the intersection segments as a line mesh here",
    )
    parser.add_argument(
        "--region", type=str, default="", help="named cell region to test"
    )
    parser.add_argument(
        "--region-other",
        type=str,
        default="",
        help="named cell region of --other to test",
    )
    parser.add_argument(
        "--max-pairs",
        type=int,
        default=0,
        metavar="N",
        help="list at most N pairs (default: 0, all); the counts still cover every pair",
    )
    parser.add_argument("--json", action="store_true", help="emit the result as JSON")


def find_intersections_cmd(args):
    # 0 clean, 1 intersections found, 2 could not check -- the `check` verb's convention.
    try:
        mesh = read(args.infile, file_format=args.input_format)
        other = None
        if args.other is not None:
            other = read(args.other, file_format=args.other_format)
        out, report = find_intersections(
            mesh,
            other,
            region=args.region,
            region_other=args.region_other,
            attach_flags=args.outfile is not None,
            segments=args.segments is not None,
            max_pairs=args.max_pairs,
            return_report=True,
        )
        if args.outfile is not None:
            write(args.outfile, out, file_format=args.output_format)
        if args.segments is not None:
            write(args.segments, report["segments"])
    except Exception as exc:  # noqa: BLE001 -- reported, exit status 2
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.json:
        obj = {k: report[k] for k in _JSON_KEYS}
        obj["pairs"] = report["pairs"]
        obj["warnings"] = report["warnings"]
        emit_json(obj)
    else:
        what = "intersections" if other is None else "contacts between the two surfaces"
        if report["clean"]:
            print(f"no {what}")
        else:
            print(
                f"{report['num_pairs']} intersecting cell pairs "
                f"({report['num_triangle_pairs']} triangle pairs, "
                f"{report['num_coplanar']} coplanar), {report['num_segments']} segments"
            )
        for w in report["warnings"]:
            print(f"  warning: {w}")
    return 0 if report["clean"] else 1
