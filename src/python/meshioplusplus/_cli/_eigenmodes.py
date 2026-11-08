"""The ``eigenmodes`` verb: the lowest eigenpairs of a surface's
Laplace-Beltrami operator (``doc/eigenmodes.md``)."""

from .._eigenmodes import laplacian_eigenmodes
from .._helpers import read, write
from ._hull import _io_args
from ._json import emit_json


def add_args(parser):
    _io_args(parser)
    parser.add_argument(
        "--num-modes", "-k", type=int, default=10, help="eigenpairs wanted"
    )
    parser.add_argument("--tolerance", type=float, default=1e-8)
    parser.add_argument(
        "--max-iterations", type=int, default=0, help="block iterations (0: 300)"
    )


def eigenmodes_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = laplacian_eigenmodes(
        mesh,
        args.num_modes,
        tolerance=args.tolerance,
        max_iterations=args.max_iterations,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return
    print("eigenvalues: " + " ".join(f"{v:.9g}" for v in report["eigenvalues"]))
    tail = "" if report["converged"] else ", not converged"
    print(f"{report['iterations']} iteration(s){tail}")
