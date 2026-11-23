"""The ``mesh-operator`` verb: the mesh's P1 operators as a sparse matrix in an
``.npz`` file (``doc/mesh_operator.md``)."""

from .._helpers import read
from .._mesh_operator import OPERATORS, mesh_operator
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(
        parser,
        "the .npz file the matrix is written to (indptr, indices, values, shape)",
        in_help="triangle surface mesh",
        output_format=False,
    )
    parser.add_argument(
        "--kind",
        type=str,
        default="laplacian",
        choices=OPERATORS,
        help="the operator (default: laplacian)",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def mesh_operator_cmd(args):
    import numpy as np

    mesh = read(args.infile, file_format=args.input_format)
    (indptr, indices, values, shape), report = mesh_operator(
        mesh, args.kind, sparse=False, return_report=True
    )
    np.savez(
        args.outfile,
        indptr=indptr,
        indices=indices,
        values=values,
        shape=np.asarray(shape),
    )
    info = {
        "output": args.outfile,
        "kind": args.kind,
        "shape": list(shape),
        "nnz": int(len(values)),
        "num_degenerate": report["num_degenerate"],
    }
    if args.json:
        emit_json(info)
        return 0
    print(
        f"wrote the {args.kind} operator ({shape[0]} x {shape[1]}, "
        f"{info['nnz']} entries) to {args.outfile}"
    )
    return 0
