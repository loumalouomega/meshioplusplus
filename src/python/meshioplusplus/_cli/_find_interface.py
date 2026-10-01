from .._helpers import _writer_map, read, reader_map, write
from .._interfaces import find_interface
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh containing part A")
    parser.add_argument("outfile", type=str, help="matched interface facet mesh")
    parser.add_argument("--region-a", required=True, help="Cell region for part A")
    parser.add_argument("--region-b", required=True, help="Cell region for part B")
    parser.add_argument(
        "--mesh-b", default=None, help="optional separate mesh containing part B"
    )
    parser.add_argument(
        "--input-format", "-i", choices=sorted(reader_map), default=None
    )
    parser.add_argument("--format-b", choices=sorted(reader_map), default=None)
    parser.add_argument(
        "--output-format", "-o", choices=sorted(_writer_map), default=None
    )
    parser.add_argument(
        "--mode", choices=("conforming", "proximity"), default="conforming"
    )
    parser.add_argument("--master", choices=("a", "b"), default="a")
    parser.add_argument("--gap-tolerance", type=float, default=0.0)
    parser.add_argument("--angle-tolerance", type=float, default=30.0)
    parser.add_argument("--overlap-tolerance", type=float, default=0.0)
    parser.add_argument("--json", action="store_true", help="emit a JSON report")


def find_interface_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    mesh_b = read(args.mesh_b, file_format=args.format_b) if args.mesh_b else None
    out, report = find_interface(
        mesh,
        args.region_a,
        args.region_b,
        mesh_b=mesh_b,
        mode=args.mode,
        master=args.master,
        gap_tolerance=args.gap_tolerance,
        angle_tolerance=args.angle_tolerance,
        overlap_tolerance=args.overlap_tolerance,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    result = {
        key: report[key]
        for key in ("num_pairs", "area", "max_gap", "unmatched_a", "unmatched_b")
    }
    result["side_a"] = report["side_a"].entries.tolist()
    result["side_b"] = report["side_b"].entries.tolist()
    if args.json:
        emit_json(result)
    else:
        print(
            f"find interface: {result['num_pairs']} facet pair(s), area {result['area']:g}, "
            f"max gap {result['max_gap']:g} written"
        )
    return 0
