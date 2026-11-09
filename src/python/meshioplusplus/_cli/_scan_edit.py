"""The scan-editing verbs: ``remove-islands``, ``flip-normals``, ``clip``,
``offset-surface``, ``thicken`` and ``hollow`` (roadmap 7.8;
``doc/remove_islands.md``, ``doc/flip_normals.md``, ``doc/clip.md``,
``doc/offset_surface.md``)."""

from .._clip import clip
from .._flip_normals import flip_normals
from .._helpers import read, write
from .._offset_surface import hollow, offset_surface, thicken
from .._remove_islands import remove_islands
from ._json import emit_json
from ._sweep_common import vec3


def _io(parser):
    parser.add_argument("infile", type=str)
    parser.add_argument("outfile", type=str)
    parser.add_argument("--input-format", "-i", type=str, default=None)
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def add_islands_args(parser):
    _io(parser)
    parser.add_argument("--connectivity", choices=("edge", "point"), default="edge")
    parser.add_argument("--min-diagonal", type=float, default=0.01)
    parser.add_argument("--min-area", type=float, default=0.0)
    parser.add_argument("--min-cells", type=int, default=0)
    parser.add_argument(
        "--absolute",
        action="store_true",
        help="the thresholds are lengths, not fractions",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="label instead of removing"
    )


def remove_islands_cmd(args):
    out, r = remove_islands(
        read(args.infile, file_format=args.input_format),
        connectivity=args.connectivity,
        min_diagonal=args.min_diagonal,
        min_area=args.min_area,
        min_cells=args.min_cells,
        relative=not args.absolute,
        dry_run=args.dry_run,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(
            {
                k: r[k]
                for k in (
                    "num_components",
                    "num_removed",
                    "num_cells_removed",
                    "num_points_removed",
                )
            }
        )
        return 0
    print(
        f"{r['num_components']} component(s); {r['num_removed']} removed "
        f"({r['num_cells_removed']} cell(s), {r['num_points_removed']} point(s))"
        + (" [dry run]" if args.dry_run else "")
    )
    return 0


def add_flip_args(parser):
    _io(parser)
    parser.add_argument("--region", type=str, default=None)


def flip_normals_cmd(args):
    out, r = flip_normals(
        read(args.infile, file_format=args.input_format),
        args.region,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(r)
        return 0
    print(f"{r['num_flipped']} cell(s) flipped")
    return 0


def add_clip_args(parser):
    _io(parser)
    parser.add_argument("--origin", type=str, default="0,0,0")
    parser.add_argument("--normal", type=str, default="0,0,1")
    parser.add_argument("--keep", choices=("below", "above"), default="below")
    parser.add_argument("--no-cap", action="store_true", help="leave the cut open")


def clip_cmd(args):
    out, r = clip(
        read(args.infile, file_format=args.input_format),
        vec3(args.origin, "clip", "origin"),
        vec3(args.normal, "clip", "normal"),
        keep=args.keep,
        cap=not args.no_cap,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(r)
        return 0
    print(
        f"{r['num_kept']} kept, {r['num_split']} split; {r['num_cap_triangles']} cap "
        f"triangle(s) over {r['num_loops']} loop(s)"
    )
    for w in r["warnings"]:
        print(f"warning: {w}")
    return 0


def _offset_args(parser, amount):
    _io(parser)
    parser.add_argument(f"--{amount}", type=float, required=True)
    parser.add_argument(
        "--smooth", type=int, default=0, help="direction smoothing passes"
    )


def add_offset_args(parser):
    _offset_args(parser, "distance")


def add_thicken_args(parser):
    _offset_args(parser, "thickness")
    parser.add_argument(
        "--side", choices=("inward", "outward", "both"), default="inward"
    )


def add_hollow_args(parser):
    _offset_args(parser, "wall")


def _offset_report(args, out, r):
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(r)
        return 0
    print(
        f"max move {r['max_move']!r}; {r['num_folded']} folded face(s), "
        f"{r['num_rim_quads']} rim quad(s)"
    )
    for w in r["warnings"]:
        print(f"warning: {w}")
    return 0


def offset_surface_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, r = offset_surface(
        mesh, args.distance, smooth_iterations=args.smooth, return_report=True
    )
    return _offset_report(args, out, r)


def thicken_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, r = thicken(
        mesh,
        args.thickness,
        side=args.side,
        smooth_iterations=args.smooth,
        return_report=True,
    )
    return _offset_report(args, out, r)


def hollow_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, r = hollow(mesh, args.wall, smooth_iterations=args.smooth, return_report=True)
    return _offset_report(args, out, r)
