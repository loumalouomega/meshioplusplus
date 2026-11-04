"""The ``register`` verb: the rigid motion laying SOURCE onto TARGET, by
iterative closest points (``doc/registration.md``)."""

from .._helpers import _reader_names, read, write
from .._registration import register_rigid
from .._transform import transform
from ._json import emit_json


def add_args(parser):
    parser.add_argument("source", type=str, help="the mesh or point cloud to move")
    parser.add_argument("target", type=str, help="the mesh to lay it onto")
    for side in ("source", "target"):
        parser.add_argument(
            f"--input-format-{side}",
            type=str,
            choices=_reader_names(),
            default=None,
            help=f"input format of the {side}",
        )
    parser.add_argument(
        "--output", type=str, default="", help="write the moved source here"
    )
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument(
        "--method",
        choices=("point_to_point", "point_to_plane"),
        default="point_to_point",
    )
    parser.add_argument("--max-iterations", type=int, default=50)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument(
        "--max-distance",
        type=float,
        default=0.0,
        help="drop pairs farther apart than this (0: keep all)",
    )
    parser.add_argument(
        "--rotate-data",
        action="store_true",
        help="rotate vector and tensor data of the written source too",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def register_cmd(args):
    source = read(args.source, file_format=args.input_format_source)
    target = read(args.target, file_format=args.input_format_target)
    matrix, report = register_rigid(
        source,
        target,
        method=args.method,
        max_iterations=args.max_iterations,
        tolerance=args.tolerance,
        max_distance=args.max_distance,
    )
    if args.output:
        moved = transform(source, matrix=matrix, rotate_vector_data=args.rotate_data)
        write(args.output, moved, file_format=args.output_format)
    if args.json:
        doc = {"output": args.output} if args.output else {}
        doc["matrix"] = [[float(v) for v in row] for row in matrix]
        doc.update(report)
        emit_json(doc)
        return 0
    for row in matrix:
        print(" ".join(f"{v:14.9g}" for v in row))
    tail = "" if report["converged"] else ", not converged"
    print(
        f"rms {report['rms']:.9g} (from {report['initial_rms']:.9g}) after "
        f"{report['iterations']} iteration(s){tail}"
    )
    return 0
