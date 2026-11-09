"""The ``step-elements`` verb: fitted primitives as exact STEP geometry
(``doc/formats/step.md``)."""

from .._helpers import read
from .._step_elements import write_step_elements
from ._json import emit_json


def add_args(parser):
    parser.add_argument(
        "infile", type=str, help="a fitted mesh (segment, then fit-primitives)"
    )
    parser.add_argument("outfile", type=str, help="the .step file to write")
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument(
        "--patches",
        type=str,
        default=None,
        help="comma-separated patch labels to write (default: every fitted patch)",
    )
    parser.add_argument(
        "--segments",
        type=int,
        default=64,
        help="facets round a cylinder, cone or torus in the writer's tessellation",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def step_elements_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    patches = (
        None
        if args.patches is None
        else [int(float(v)) for v in args.patches.split(",") if v.strip()]
    )
    r = write_step_elements(args.outfile, mesh, patches=patches, segments=args.segments)
    if args.json:
        emit_json(r)
        return 0
    print(
        f"{r['num_planes']} plane(s), {r['num_spheres']} sphere(s), "
        f"{r['num_cylinders']} cylinder(s), {r['num_cones']} cone(s), "
        f"{r['num_tori']} torus/tori; {r['num_skipped']} skipped"
    )
    for w in r["warnings"]:
        print(f"warning: {w}")
    return 0
