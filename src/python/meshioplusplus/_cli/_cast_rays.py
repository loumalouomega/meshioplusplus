"""The ``cast-rays`` verb: the nearest surface hit of a batch of rays
(``doc/ray_queries.md``)."""

from .._cast_rays import _directions_of, _hits_mesh, _xyz, cast_rays
from .._helpers import write
from ._common import add_transfer_io_args, read_source_target
from ._json import emit_json


def add_args(parser):
    add_transfer_io_args(
        parser,
        "vertex mesh of the hit points (ray:index, ray:distance, ray:cell, ray:normal)",
    )
    parser.add_argument(
        "--direction",
        type=str,
        default=None,
        metavar="X,Y,Z",
        help="one direction for every ray (default: the point array 'directions' of the rays file)",
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=0.0,
        help="ignore hits farther than this (default: no limit)",
    )
    parser.add_argument(
        "--region", type=str, default=None, help="restrict to this cell region"
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def cast_rays_cmd(args):
    surface, rays = read_source_target(args)
    direction = None
    if args.direction is not None:
        try:
            direction = [float(v) for v in args.direction.split(",")]
        except ValueError:
            raise SystemExit("cast-rays: --direction expects three numbers") from None
        if len(direction) != 3:
            raise SystemExit("cast-rays: --direction expects three numbers")
    try:
        directions = _directions_of(rays, direction)
    except ValueError as exc:
        raise SystemExit(f"cast-rays: {exc}") from None
    res = cast_rays(
        surface,
        _xyz(rays.points),
        directions,
        max_distance=args.max_distance,
        region=args.region,
    )
    write(args.outfile, _hits_mesh(res), file_format=args.output_format)
    n = len(res["hit"])
    if args.json:
        emit_json({"output": args.outfile, "num_rays": n, "num_hits": res["num_hits"]})
        return 0
    print(f"{res['num_hits']} of {n} rays hit the surface")
    return 0
