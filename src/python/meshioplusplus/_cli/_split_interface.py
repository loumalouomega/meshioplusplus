import json

from .._helpers import _reader_names, _writer_names, read, write
from .._interfaces import split_interface
from .._regions import Region
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh containing the interface")
    parser.add_argument("outfile", type=str, help="split mesh to write")
    parser.add_argument(
        "--side-region", default=None, help="name of a Side region on the input mesh"
    )
    parser.add_argument(
        "--side-entries",
        default=None,
        help='JSON list of [global_cell, local_facet] pairs, e.g. "[[4, 2], [7, 0]]"',
    )
    parser.add_argument("--input-format", "-i", choices=_reader_names(), default=None)
    parser.add_argument("--output-format", "-o", choices=_writer_names(), default=None)
    parser.add_argument("--add-cohesive", action="store_true")
    parser.add_argument("--json", action="store_true", help="emit a JSON report")


def split_interface_cmd(args):
    if (args.side_region is None) == (args.side_entries is None):
        raise ValueError("supply exactly one of --side-region and --side-entries")
    mesh = read(args.infile, file_format=args.input_format)
    side = args.side_region
    if args.side_entries is not None:
        side = Region("cli:interface", "side", json.loads(args.side_entries))
    out, report = split_interface(
        mesh, side, add_cohesive=args.add_cohesive, return_report=True
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(report)
    else:
        print(
            f"split interface: duplicated {report['num_duplicated_points']} point(s), "
            f"added {report['num_cohesive_cells']} cohesive cell(s)"
        )
    return 0
