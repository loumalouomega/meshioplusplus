from .._helpers import _writer_map, read, reader_map, write
from .._refine import refine


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="refined mesh to be written to")
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
        "--levels",
        type=int,
        default=1,
        help=(
            "how many times to subdivide; each level splits a triangle/quad "
            "into 4 and a tetra/wedge/hexahedron into 8 (default: 1)"
        ),
    )
    parser.add_argument(
        "--record-parent-ids",
        action="store_true",
        help="attach refine:parent_cell cell_data of the original cell indices",
    )
    # At most one selector; with none, every cell is refined.
    parser.add_argument(
        "--cells",
        type=str,
        default=None,
        help=(
            "comma-separated global (block-major) cell indices to refine; "
            "everything else is only split as far as conformity demands"
        ),
    )
    parser.add_argument(
        "--region",
        type=str,
        default=None,
        help=(
            "name of a region to refine (a cell region selects its cells, a "
            "point region every cell touching it; a side region is an error)"
        ),
    )
    parser.add_argument(
        "--where",
        type=str,
        default=None,
        metavar="EXPR",
        help=(
            "refine the cells satisfying a threshold on a scalar cell_data "
            'array, e.g. --where "quality:scaled_jacobian < 0.3"'
        ),
    )
    parser.add_argument(
        "--closure",
        type=str,
        choices=("redgreen", "propagate", "balanced"),
        default="redgreen",
        help=(
            "how to resolve hanging nodes: redgreen keeps the extra refinement "
            "local (default), propagate works for every cell type but reaches "
            "the whole edge-connected component, balanced KEEPS the hanging "
            "nodes and only enforces 2:1 balance (output is not conforming)"
        ),
    )
    parser.add_argument(
        "--record-levels",
        action="store_true",
        help="attach refine:level cell_data of each cell's refinement depth",
    )
    parser.add_argument(
        "--record-hierarchy",
        action="store_true",
        help=(
            "attach refine:cell_id/refine:parent_id cell_data -- the persistent "
            "parent/child hierarchy a multigrid caller resolves across the "
            "sequence of meshes it keeps; also forces refine:entity to be "
            "attached even when the closure leaves no hanging node"
        ),
    )


def refine_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    cells = None
    if args.cells is not None:
        cells = [int(x) for x in args.cells.split(",") if x.strip()]
    out = refine(
        mesh,
        levels=args.levels,
        record_parent_ids=args.record_parent_ids,
        cells=cells,
        region=args.region,
        where=args.where,
        closure=args.closure,
        record_levels=args.record_levels,
        record_hierarchy=args.record_hierarchy,
    )
    write(args.outfile, out, file_format=args.output_format)
    return 0
