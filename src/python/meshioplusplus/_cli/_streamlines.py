"""The ``streamlines`` verb: lines along a point vector field, as a mesh of
segments (``doc/streamlines.md``)."""

from .._helpers import read, write
from .._streamlines import streamlines
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(
        parser,
        "output mesh of line segments",
        in_help="mesh carrying the vector point array",
        output_format=False,
    )
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument(
        "--array", type=str, required=True, help="the 2- or 3-component point array"
    )
    parser.add_argument(
        "--seeds",
        type=int,
        default=40,
        help="about this many seeds, ranked at equal measure (default 40)",
    )
    parser.add_argument(
        "--seed-points",
        type=str,
        default="",
        metavar="X,Y,Z[;X,Y,Z...]",
        help="explicit seeds, semicolon-separated triples (instead of --seeds)",
    )
    parser.add_argument(
        "--direction",
        choices=("both", "forward", "backward"),
        default="both",
        help="which way from each seed (default both)",
    )
    parser.add_argument(
        "--integrator",
        choices=("rk4", "rk45"),
        default="rk4",
        help="fixed-step RK4, or adaptive Dormand-Prince 5(4) (default rk4)",
    )
    parser.add_argument(
        "--step",
        type=float,
        default=0.0,
        help="the step as a length at the fastest speed (default 0.005 of the diagonal)",
    )
    parser.add_argument(
        "--max-length",
        type=float,
        default=0.5,
        help="the longest path each way, a fraction of the diagonal (default 0.5)",
    )
    parser.add_argument("--max-steps", type=int, default=0)
    parser.add_argument(
        "--terminal-speed",
        type=float,
        default=0.0,
        help="a line stops where the speed falls to this (default 1e-9 of the fastest)",
    )
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def _seed_points(text):
    pts = []
    for triple in text.split(";"):
        if not triple.strip():
            continue
        try:
            v = [float(c) for c in triple.split(",")]
        except ValueError:
            raise SystemExit("streamlines: --seed-points expects x,y,z;x,y,z") from None
        if len(v) == 2:
            v.append(0.0)
        if len(v) != 3:
            raise SystemExit("streamlines: --seed-points expects x,y,z;x,y,z")
        pts.append(v)
    return pts or None


def streamlines_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = streamlines(
        mesh,
        args.array,
        seeds=args.seeds,
        seed_points=_seed_points(args.seed_points),
        direction=args.direction,
        integrator=args.integrator,
        step=args.step,
        max_length=args.max_length,
        max_steps=args.max_steps,
        terminal_speed=args.terminal_speed,
        tolerance=args.tolerance,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(
            {
                "output": args.outfile,
                "num_lines": report["num_lines"],
                "num_points": len(out.points),
                "num_seeds_outside": report["num_seeds_outside"],
                "dim": report["dim"],
            }
        )
        return 0
    text = f"traced {report['num_lines']} streamline(s), {len(out.points)} point(s)"
    if report["num_seeds_outside"]:
        text += f", {report['num_seeds_outside']} seed(s) outside the mesh"
    print(text)
    return 0
