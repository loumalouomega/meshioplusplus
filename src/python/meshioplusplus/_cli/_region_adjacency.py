from .._helpers import _writer_map, read, reader_map, write
from .._interfaces import region_adjacency
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh to analyze")
    parser.add_argument("outfile", type=str, help="shared-facet mesh to write")
    parser.add_argument(
        "--input-format", "-i", choices=sorted(reader_map), default=None
    )
    parser.add_argument(
        "--output-format", "-o", choices=sorted(_writer_map), default=None
    )
    parser.add_argument(
        "--regions",
        nargs="+",
        default=None,
        help="Cell region names (default: every Cell region, or each block if none exist)",
    )
    parser.add_argument("--json", action="store_true", help="emit a JSON summary")
    parser.add_argument("--quiet", "-q", action="store_true", help="suppress summary")


def region_adjacency_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    regions = None
    if args.regions:
        regions = [name for value in args.regions for name in value.split(",") if name]
    out = region_adjacency(mesh, regions)
    write(args.outfile, out, file_format=args.output_format)
    nfacets = sum(len(block.data) for block in out.cells)
    area = sum(float(values.sum()) for values in out.cell_data["interface:measure"])
    if args.json:
        emit_json({"facets": nfacets, "measure": area})
    elif not args.quiet:
        print(f"region adjacency: {nfacets} shared facets, measure {area:g} written")
    return 0
