"""The ``register`` verb: the rigid motion laying SOURCE onto TARGET, by
iterative closest points (``doc/registration.md``)."""

from .._helpers import _reader_names, read, write
from .._registration import (
    auto_align,
    find_symmetry_plane,
    register_landmarks,
    register_rigid,
)
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
        "--initial",
        choices=("identity", "auto"),
        default="identity",
        help="the start: the identity, or 'auto', the best of the 24 principal-frame "
        "candidates",
    )
    parser.add_argument(
        "--symmetric-check",
        action="store_true",
        help="with --initial auto, warn when the part looks symmetric",
    )
    parser.add_argument(
        "--source-region",
        type=str,
        default="",
        help="cell region of the source to pair",
    )
    parser.add_argument(
        "--target-region",
        type=str,
        default="",
        help="cell region of the target to match",
    )
    parser.add_argument(
        "--adaptive-trim",
        type=float,
        default=0.0,
        metavar="K",
        help="drop pairs beyond K median pair distances at each iteration (0: off)",
    )
    parser.add_argument(
        "--rotate-data",
        action="store_true",
        help="rotate vector and tensor data of the written source too",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def _triples(text, what):
    rows = []
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        parts = [float(v) for v in item.split(",")]
        if len(parts) != 3:
            raise SystemExit(
                f"meshio++: register-landmarks: {what} takes x,y,z;x,y,z;..."
            )
        rows.append(parts)
    return rows


def add_landmark_args(parser):
    parser.add_argument(
        "source_points", type=str, help="picked source points 'x,y,z;x,y,z;...'"
    )
    parser.add_argument(
        "target_points", type=str, help="their target points 'x,y,z;x,y,z;...'"
    )
    parser.add_argument(
        "--weights", type=str, default=None, help="one weight per pair 'w,w,...'"
    )
    parser.add_argument("--mesh", type=str, default="", help="a mesh to move")
    parser.add_argument(
        "--output", type=str, default="", help="write the moved --mesh here"
    )
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def register_landmarks_cmd(args):
    weights = None
    if args.weights is not None:
        weights = [float(v) for v in args.weights.split(",") if v.strip()]
    matrix, report = register_landmarks(
        _triples(args.source_points, "the source"),
        _triples(args.target_points, "the target"),
        weights,
    )
    if args.mesh and args.output:
        write(
            args.output,
            transform(read(args.mesh), matrix=matrix),
            file_format=args.output_format,
        )
    if args.json:
        doc = {"output": args.output} if args.mesh and args.output else {}
        doc["matrix"] = [[float(v) for v in row] for row in matrix]
        doc.update(report)
        emit_json(doc)
        return
    for row in matrix:
        print(" ".join(f"{v:14.9g}" for v in row))
    print(f"rms {report['rms']:.9g} over {len(report['residuals'])} pair(s)")


def add_auto_align_args(parser):
    parser.add_argument("infile", type=str, help="the mesh to pose")
    parser.add_argument("outfile", type=str, help="the posed mesh")
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument(
        "--angle-tolerance",
        type=float,
        default=10.0,
        metavar="DEG",
        help="degrees within which a normal agrees with a direction (default: 10)",
    )
    parser.add_argument(
        "--subdivisions",
        type=int,
        default=3,
        help="the Gauss-sphere icosahedron's splits (default: 3, 642 bins)",
    )
    parser.add_argument("--region", type=str, default="", help="the voting cell region")
    parser.add_argument(
        "--rotate-data", action="store_true", help="rotate vector and tensor data too"
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def auto_align_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = auto_align(
        mesh,
        angle_tolerance=args.angle_tolerance,
        subdivisions=args.subdivisions,
        region=args.region,
        rotate_vector_data=args.rotate_data,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        doc = dict(report)
        doc["matrix"] = [[float(v) for v in row] for row in report["matrix"]]
        emit_json(doc)
        return
    print(
        f"posed on its {'largest face' if report['mode'] == 'face' else 'axis'}: "
        f"{report['share']:.4g} of the area agrees (face {report['face_share']:.4g}, "
        f"axis {report['axis_share']:.4g})"
    )
    for w in report["warnings"]:
        print(f"warning: {w}")


def add_symmetry_args(parser):
    parser.add_argument("infile", type=str, help="the mesh or point cloud")
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument(
        "--seed-origin", type=str, default=None, help="a seed plane's point 'x,y,z'"
    )
    parser.add_argument(
        "--seed-normal", type=str, default=None, help="the seed plane's normal 'x,y,z'"
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.0,
        help="a mirrored sample this close counts as symmetric (0: 0.005 of the "
        "diagonal)",
    )
    parser.add_argument(
        "--min-share",
        type=float,
        default=0.9,
        help="refuse a plane with a smaller share of symmetric samples (default: 0.9)",
    )
    parser.add_argument("--max-samples", type=int, default=4000)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def symmetry_plane_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    seed = None
    if (args.seed_origin is None) != (args.seed_normal is None):
        raise SystemExit(
            "meshio++: symmetry-plane: --seed-origin and --seed-normal go together"
        )
    if args.seed_origin is not None:
        seed = (
            [float(v) for v in args.seed_origin.split(",")],
            [float(v) for v in args.seed_normal.split(",")],
        )
    r = find_symmetry_plane(
        mesh,
        seed=seed,
        tolerance=args.tolerance,
        min_share=args.min_share,
        max_samples=args.max_samples,
    )
    if args.json:
        emit_json(r)
        return
    o, n = r["origin"], r["normal"]
    print(
        f"plane through ({o[0]:.9g}, {o[1]:.9g}, {o[2]:.9g}) normal "
        f"({n[0]:.9g}, {n[1]:.9g}, {n[2]:.9g})"
    )
    print(
        f"rms {r['rms']:.9g}, {r['share']:.4g} of the samples within {r['tolerance']:.6g}"
    )


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
        initial="auto" if args.initial == "auto" else None,
        symmetric_check=args.symmetric_check,
        source_region=args.source_region,
        target_region=args.target_region,
        adaptive_trim=args.adaptive_trim,
    )
    if args.output:
        moved = transform(source, matrix=matrix, rotate_vector_data=args.rotate_data)
        write(args.output, moved, file_format=args.output_format)
    if args.json:
        doc = {"output": args.output} if args.output else {}
        doc["matrix"] = [[float(v) for v in row] for row in matrix]
        doc.update(report)
        emit_json(doc)
        return
    for row in matrix:
        print(" ".join(f"{v:14.9g}" for v in row))
    tail = "" if report["converged"] else ", not converged"
    print(
        f"rms {report['rms']:.9g} (from {report['initial_rms']:.9g}) after "
        f"{report['iterations']} iteration(s){tail}"
    )
    for w in report["warnings"]:
        print(f"warning: {w}")
