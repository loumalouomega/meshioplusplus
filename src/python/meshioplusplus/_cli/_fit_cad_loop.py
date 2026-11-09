"""The reverse-engineering loop verbs (roadmap 7.7): ``face-support``,
``fit-cad-iterate`` and ``constrain-primitives`` (``doc/fit_cad_loop.md``,
``doc/constrain_primitives.md``)."""

from .._fit_cad_loop import constrain_primitives, face_support, fit_cad_iterate
from .._helpers import read, write
from ._json import emit_json


def add_face_support_args(parser):
    parser.add_argument("solid", type=str, help="a fit_cad re-mesh")
    parser.add_argument("scan", type=str)
    parser.add_argument("outfile", type=str, nargs="?", default=None)
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument("--tolerance", type=float, default=0.0)
    parser.add_argument("--measured-share", type=float, default=0.9)
    parser.add_argument("--inferred-share", type=float, default=0.1)
    parser.add_argument("--hole-data", type=str, default="repair:hole")
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def face_support_cmd(args):
    out, r = face_support(
        read(args.solid, file_format=args.input_format),
        read(args.scan, file_format=args.input_format),
        tolerance=args.tolerance,
        measured_share=args.measured_share,
        inferred_share=args.inferred_share,
        hole_data=args.hole_data,
        return_report=True,
    )
    if args.outfile:
        write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(r)
        return 0
    print(
        f"{r['num_measured']} measured, {r['num_partial']} partial, "
        f"{r['num_inferred']} inferred face(s)"
    )
    return 0


def add_iterate_args(parser):
    parser.add_argument("scan", type=str)
    parser.add_argument("outfile", type=str)
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument(
        "--step", type=str, default=None, help="write the best solid as STEP"
    )
    parser.add_argument(
        "--part-class",
        choices=("any", "analytic", "cast", "freeform"),
        default="any",
    )
    parser.add_argument("--target-within", type=float, default=0.0)
    parser.add_argument("--target-explained", type=float, default=0.0)
    parser.add_argument("--segment-tolerance", type=float, default=0.0)
    parser.add_argument("--fit-tolerance", type=float, default=0.0)
    parser.add_argument("--max-rounds", type=int, default=4)
    parser.add_argument("--worst-faces", type=int, default=3)
    parser.add_argument("--refine-factor", type=float, default=0.25)
    parser.add_argument("--backend", choices=("auto", "occt", "native"), default="auto")
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def fit_cad_iterate_cmd(args):
    out, r = fit_cad_iterate(
        read(args.scan, file_format=args.input_format),
        args.step,
        part_class=args.part_class,
        target_within=args.target_within,
        target_explained=args.target_explained,
        segment_tolerance=args.segment_tolerance,
        fit_tolerance=args.fit_tolerance,
        max_rounds=args.max_rounds,
        worst_faces=args.worst_faces,
        refine_factor=args.refine_factor,
        backend=args.backend,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    r.pop("skin")
    if args.json:
        emit_json(r)
        return 0
    print(
        f"{len(r['rounds'])} round(s), best {r['best_round']} ({r['stop_reason']}); "
        f"target {'met' if r['target_met'] else 'not met'}"
    )
    return 0


def _relation(text):
    parts = text.split(":")
    if len(parts) != 3:
        raise SystemExit(
            f"constrain-primitives: a relation is PATCH_A:PATCH_B:RELATION, not {text!r}"
        )
    return int(parts[0]), int(parts[1]), parts[2]


def add_constrain_args(parser):
    parser.add_argument("infile", type=str)
    parser.add_argument("outfile", type=str)
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument("--face-tolerance", type=float, default=0.0)
    parser.add_argument("--angle-tolerance", type=float, default=-1.0)
    parser.add_argument("--distance-tolerance", type=float, default=-1.0)
    parser.add_argument("--min-fraction", type=float, default=-1.0)
    parser.add_argument("--force", action="append", default=[], metavar="A:B:REL")
    parser.add_argument("--forbid", action="append", default=[], metavar="A:B:REL")
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def constrain_primitives_cmd(args):
    out, r = constrain_primitives(
        read(args.infile, file_format=args.input_format),
        face_tolerance=args.face_tolerance,
        angle_tolerance=args.angle_tolerance,
        distance_tolerance=args.distance_tolerance,
        min_fraction=args.min_fraction,
        force=[_relation(t) for t in args.force],
        forbid=[_relation(t) for t in args.forbid],
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        r["unmet"] = [list(u) for u in r["unmet"]]
        emit_json(r)
        return 0
    print(
        f"{r['num_constrained']} primitive(s) constrained in {r['num_groups']} "
        f"group(s); {len(r['relations'])} relation(s), {len(r['unmet'])} forced "
        "relation(s) unmet"
    )
    return 0


def add_profile_args(parser):
    parser.add_argument("sketch", type=str, help="a fit_curves sketch")
    parser.add_argument("outfile", type=str)
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument(
        "--features",
        type=str,
        required=True,
        help="';'-separated: extrude:D[:symmetric], revolve:OX,OY,OZ:DX,DY,DZ[:ANGLE], "
        "fillet:R, chamfer:D, shell:W",
    )
    parser.add_argument("--step", type=str, default=None)
    parser.add_argument("--scan", type=str, default=None)
    parser.add_argument("--provider", type=str, default="")
    parser.add_argument("--deflection", type=float, default=0.0)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def profile_solid_cmd(args):
    from .._profile_solid import profile_solid

    scan = read(args.scan, file_format=args.input_format) if args.scan else None
    out, r = profile_solid(
        read(args.sketch, file_format=args.input_format),
        args.features,
        args.step,
        provider=args.provider,
        deflection=args.deflection,
        scan=scan,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(r)
        return 0
    note = "" if r["valid"] else " (not a valid solid)"
    print(f"{r['num_faces']} face(s), volume {r['volume']!r}{note} [{r['provider']}]")
    for w in r["warnings"]:
        print(f"warning: {w}")
    return 0
