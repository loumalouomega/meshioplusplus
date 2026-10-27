import sys

from .._check_conformity import check_conformity
from .._helpers import _writer_map, read, reader_map, write
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be checked")
    parser.add_argument(
        "outfile",
        type=str,
        nargs="?",
        default=None,
        help="optional: write the mesh with the conformity:hanging point data here",
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
        "--tolerance",
        type=float,
        default=1e-9,
        help="distances at or below this fraction of the facets' bounding-box "
        "diagonal count as zero (default: 1e-9)",
    )
    parser.add_argument("--json", action="store_true", help="emit the result as JSON")


def format_conformity(report):
    """The human-readable summary ``check-conformity`` prints."""
    verdict = "conforming" if report["conforming"] else "NOT conforming"
    lines = [
        f"{verdict} ({report['dimension']}-D cells examined)",
        f"  hanging nodes:       {report['num_hanging_nodes']}",
        f"  T-junctions:         {report['num_t_junctions']}",
        f"  partial overlaps:    {report['num_partial_overlaps']}",
        f"  unwelded facets:     {report['num_unwelded_facets']}",
        f"  boundary facets:     {report['num_boundary_facets']}",
        f"  non-manifold facets: {report['num_non_manifold_facets']}",
        f"  interface mode:      {report['interface_mode']}",
    ]
    lines += [f"  warning: {w}" for w in report["warnings"]]
    return "\n".join(lines) + "\n"


def check_conformity_cmd(args):
    # 0 conforming, 1 not, 2 could not check -- the `check` verb's convention.
    try:
        mesh = read(args.infile, file_format=args.input_format)
        out, report = check_conformity(
            mesh,
            tolerance=args.tolerance,
            attach_flags=args.outfile is not None,
            return_report=True,
        )
        if args.outfile is not None:
            write(args.outfile, out, file_format=args.output_format)
    except Exception as exc:  # noqa: BLE001 -- reported, exit status 2
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.json:
        emit_json(report)
    else:
        print(format_conformity(report), end="")
    return 0 if report["conforming"] else 1
