"""The ``graph`` verb: write a mesh's node or cell graph as a METIS or Scotch
graph file (``doc/adjacency.md``)."""

from .._adjacency import write_graph
from .._helpers import read
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(
        parser,
        "graph file: .graph or .metis for METIS, .grf for Scotch",
        in_help="mesh file",
        output_format=False,
    )
    parser.add_argument(
        "--kind",
        choices=("node", "clique", "cell"),
        default="node",
        help="node: nodes joined by a cell edge (default); clique: nodes sharing a "
        "cell, the FEM sparsity pattern; cell: cells sharing a facet, the dual graph",
    )
    parser.add_argument(
        "--format",
        choices=("metis", "scotch"),
        default=None,
        help="graph-file dialect (default: from the extension)",
    )
    parser.add_argument("--json", action="store_true", help="emit the summary as JSON")


def graph_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    report = write_graph(args.outfile, mesh, args.kind, file_format=args.format)
    if args.json:
        emit_json({"output": args.outfile, "kind": args.kind, **report})
        return
    print(
        f"{args.kind} graph -> {args.outfile} ({report['format']}): "
        f"{report['num_vertices']} vertices, {report['num_edges']} edges, "
        f"max degree {report['max_degree']}"
    )
