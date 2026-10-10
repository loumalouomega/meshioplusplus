"""``meshioplusplus snapshot``: the software rasterizer from the command line.

The flags are the native CLI's (``src/cpp/cli/main.cpp``, ``render_flag_specs``)
so the two verbs take the same command lines.
"""

import shutil
import sys

from .._helpers import read, reader_map
from .._render import _core_module, _environment_color_depth
from .._render import parse_color as _rgba
from .._render import parse_cutaway, render_text, snapshot


def _numbers(text, flag):
    out = []
    for item in text.split(","):
        if not item.strip():
            out.append(None)
            continue
        try:
            out.append(float(item))
        except ValueError:
            raise SystemExit(
                f"meshio++: --{flag} expects numbers separated by commas, not '{text}'"
            ) from None
    return out


def _parse_clip(text):
    clip = _numbers(text, "clip")
    if len(clip) != 2:
        raise SystemExit(
            "meshio++: --clip expects LOW,HIGH percentiles (either may be empty)"
        )
    return tuple(clip)


def _parse_levels(text):
    levels = _numbers(text, "iso-levels")
    if any(v is None for v in levels):
        raise SystemExit("meshio++: --iso-levels expects numbers separated by commas")
    return levels


def add_render_args(parser):
    """The render flags shared by every software-rendering verb."""
    camera = parser.add_argument_group("camera")
    camera.add_argument(
        "--view",
        choices=["iso", "+x", "-x", "+y", "-y", "+z", "-z"],
        default=None,
        help="named view: the camera sits on that side (overrides --azimuth/--elevation)",
    )
    camera.add_argument("--azimuth", type=float, default=None, metavar="DEG")
    camera.add_argument("--elevation", type=float, default=None, metavar="DEG")
    camera.add_argument("--roll", type=float, default=None, metavar="DEG")
    camera.add_argument("--perspective", action="store_true", help="a pinhole camera")
    camera.add_argument("--fov", type=float, default=None, metavar="DEG")
    camera.add_argument("--zoom", type=float, default=None)
    camera.add_argument("--pan-x", type=float, default=None)
    camera.add_argument("--pan-y", type=float, default=None)
    look = parser.add_argument_group("appearance")
    look.add_argument("--shading", choices=["none", "flat", "smooth"], default="flat")
    look.add_argument("--one-sided", action="store_true", help="leave back faces unlit")
    look.add_argument("--ambient", type=float, default=None)
    look.add_argument("--split-angle", type=float, default=None, metavar="DEG")
    look.add_argument("--edges", choices=["none", "all", "feature"], default="none")
    look.add_argument("--feature-angle", type=float, default=None, metavar="DEG")
    for flag in ("edge-color", "fill", "line-color", "background"):
        look.add_argument(f"--{flag}", type=_rgba, default=None, metavar="#RRGGBB[AA]")
    look.add_argument("--point-radius", type=float, default=None)
    look.add_argument("--supersample", type=int, choices=[1, 2, 4], default=None)
    look.add_argument("--axes", action="store_true", help="draw the world axes")
    look.add_argument("--scale-bar", action="store_true", help="draw a scale bar")
    color = parser.add_argument_group("colouring")
    color.add_argument(
        "--color-by", type=str, default=None, help="data array to colour by"
    )
    color.add_argument("--component", type=int, default=None)
    color.add_argument("--cmap", type=str, default=None)
    color.add_argument("--vmin", type=float, default=None)
    color.add_argument("--vmax", type=float, default=None)
    color.add_argument("--nan-color", type=_rgba, default=None, metavar="#RRGGBB[AA]")
    color.add_argument("--colorbar", action="store_true")
    field = parser.add_argument_group("field rendering")
    field.add_argument(
        "--reduce", choices=["mises", "hydrostatic", "principal"], default=None
    )
    field.add_argument(
        "--expr", type=str, default=None, help="colour by a data_calc expression"
    )
    field.add_argument(
        "--clip", type=str, default=None, metavar="LOW,HIGH", help="percentiles"
    )
    field.add_argument(
        "--symmetric", action="store_true", help="range symmetric about zero"
    )
    field.add_argument("--scale", choices=["linear", "log", "symlog"], default=None)
    field.add_argument("--scale-threshold", type=float, default=None)
    field.add_argument(
        "--categorical", action="store_true", help="integer data as categories"
    )
    field.add_argument(
        "--color-regions", action="store_true", help="colour by cell region"
    )
    field.add_argument("--category-edges", action="store_true")
    field.add_argument("--isolines", type=int, default=None, metavar="N")
    field.add_argument("--iso-levels", type=str, default=None, metavar="A,B,...")
    field.add_argument("--iso-color", type=_rgba, default=None, metavar="#RRGGBB[AA]")
    field.add_argument(
        "--vectors", type=str, default=None, help="vector point array for arrows"
    )
    field.add_argument("--vector-count", type=int, default=None)
    field.add_argument("--vector-length", type=float, default=None)
    field.add_argument(
        "--vector-color", type=_rgba, default=None, metavar="#RRGGBB[AA]"
    )
    field.add_argument(
        "--streamlines", type=str, default=None, help="vector point array to follow"
    )
    field.add_argument("--stream-seeds", type=int, default=None)
    field.add_argument("--stream-length", type=float, default=None)
    field.add_argument(
        "--stream-color", type=_rgba, default=None, metavar="#RRGGBB[AA]"
    )
    field.add_argument(
        "--warp", type=str, default=None, help="displacement point array"
    )
    field.add_argument("--warp-scale", type=float, default=None)
    field.add_argument("--warp-outline", action="store_true")
    field.add_argument(
        "--outline-color", type=_rgba, default=None, metavar="#RRGGBB[AA]"
    )
    field.add_argument(
        "--diagnostic",
        choices=[
            "none",
            "quality",
            "inverted",
            "degenerate",
            "orientation",
            "free-edges",
            "edge-length",
        ],
        default=None,
    )
    field.add_argument("--quality-metric", type=str, default=None)
    cut = parser.add_argument_group("cut-aways")
    cut.add_argument(
        "--cutaway",
        action="append",
        default=None,
        metavar="PLANE",
        help="clip away a half-space (twice at most): PX,PY,PZ,NX,NY,NZ (a point and "
        "the normal of the side kept) or AXIS:OFFSET with AXIS one of +x -x +y -y +z -z",
    )
    cut.add_argument("--cutaway-tint", type=_rgba, default=None, metavar="#RRGGBB[AA]")


def render_options(args):
    """The keyword arguments of :mod:`meshioplusplus._render` from the flags."""
    if args.view is not None and (
        args.azimuth is not None or args.elevation is not None
    ):
        raise SystemExit(
            "meshio++: --view and --azimuth/--elevation are mutually exclusive"
        )
    if args.fov is not None and not args.perspective:
        raise SystemExit("meshio++: --fov requires --perspective")
    diagnostic = args.diagnostic or "none"
    if args.quality_metric is not None and diagnostic != "quality":
        raise SystemExit("meshio++: --quality-metric requires --diagnostic quality")
    mapped = (
        args.color_by is not None
        or args.expr is not None
        or diagnostic in ("quality", "edge-length")
    )
    if not mapped:
        what = "--color-by, --expr or a quality or edge-length diagnostic"
        for flag in (
            "component",
            "cmap",
            "vmin",
            "vmax",
            "nan_color",
            "clip",
            "scale",
            "scale_threshold",
        ):
            if getattr(args, flag) is not None:
                raise SystemExit(
                    f"meshio++: --{flag.replace('_', '-')} requires {what}"
                )
        if args.colorbar or args.symmetric:
            raise SystemExit(f"meshio++: --colorbar and --symmetric require {what}")
    options = {
        "view": args.view,
        "azimuth": args.azimuth,
        "elevation": args.elevation,
        "roll": args.roll,
        "projection": "perspective" if args.perspective else None,
        "fov": args.fov,
        "zoom": args.zoom,
        "shading": args.shading,
        "two_sided": not args.one_sided,
        "ambient": args.ambient,
        "split_angle": args.split_angle,
        "edges": args.edges,
        "feature_angle": args.feature_angle,
        "edge_color": args.edge_color,
        "fill_color": args.fill,
        "line_color": args.line_color,
        "background": args.background,
        "point_radius": args.point_radius,
        "supersample": args.supersample,
        "axes": args.axes,
        "scale_bar": args.scale_bar,
        "color_by": args.color_by,
        "component": args.component,
        "cmap": args.cmap,
        "vmin": args.vmin,
        "vmax": args.vmax,
        "nan_color": args.nan_color,
        "colorbar": args.colorbar if mapped else None,
        "reduce": args.reduce,
        "expr": args.expr,
        "clip": _parse_clip(args.clip) if args.clip is not None else None,
        "symmetric": args.symmetric or None,
        "scale": args.scale,
        "scale_threshold": args.scale_threshold,
        "categorical": args.categorical or None,
        "color_regions": args.color_regions or None,
        "category_edges": args.category_edges or None,
        "isolines": args.isolines,
        "iso_levels": (
            _parse_levels(args.iso_levels) if args.iso_levels is not None else None
        ),
        "iso_color": args.iso_color,
        "vectors": args.vectors,
        "vector_count": args.vector_count,
        "vector_length": args.vector_length,
        "vector_color": args.vector_color,
        "streamlines": args.streamlines,
        "stream_seeds": args.stream_seeds,
        "stream_length": args.stream_length,
        "stream_color": args.stream_color,
        "warp": args.warp,
        "warp_scale": args.warp_scale,
        "warp_outline": args.warp_outline or None,
        "outline_color": args.outline_color,
        "diagnostic": None if diagnostic == "none" else diagnostic.replace("-", "_"),
        "quality_metric": args.quality_metric,
        "cutaway": args.cutaway,
        "cutaway_tint": args.cutaway_tint,
    }
    if args.cutaway_tint is not None and not args.cutaway:
        raise SystemExit("meshio++: --cutaway-tint requires --cutaway")
    if args.cutaway:
        try:
            options["cutaway"] = [list(parse_cutaway(p)) for p in args.cutaway]
        except ValueError as exc:
            raise SystemExit(f"meshio++: --cutaway: {exc}") from None
    if args.pan_x is not None or args.pan_y is not None:
        options["pan"] = (args.pan_x or 0.0, args.pan_y or 0.0)
    return {k: v for k, v in options.items() if v is not None}


def _enable_terminal_output():
    """On Windows, the native CLI's ``enable_terminal_output``: UTF-8 output and
    virtual-terminal processing on the console. False when the console refuses
    (before Windows 10 1511), so the caller can drop the colour escapes."""
    if sys.platform != "win32":
        return True
    try:
        import ctypes

        kernel32 = ctypes.windll.kernel32
        kernel32.SetConsoleOutputCP(65001)
        handle = kernel32.GetStdHandle(-11)  # STD_OUTPUT_HANDLE
        mode = ctypes.c_uint32()
        if not kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            return True  # not a console: a pipe or a file needs nothing
        return bool(kernel32.SetConsoleMode(handle, mode.value | 0x0004))
    except (AttributeError, OSError):
        return True


def _write_text(text):
    # Bytes, not str: a legacy console encoding (cp1252) cannot encode the
    # block characters, and the native CLI writes UTF-8 regardless.
    buffer = getattr(sys.stdout, "buffer", None)
    if buffer is None:
        sys.stdout.write(text)
    else:
        sys.stdout.flush()
        buffer.write(text.encode("utf-8"))
        buffer.flush()
    sys.stdout.flush()


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument(
        "outfile",
        type=str,
        help="output: .png, .txt, .ansi, .html or .cast; - draws in this terminal",
    )
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--width", type=int, default=800, help="PNG width (default: 800)"
    )
    parser.add_argument(
        "--height", type=int, default=600, help="PNG height (default: 600)"
    )
    parser.add_argument("--cols", type=int, default=None, help="text width in cells")
    parser.add_argument("--rows", type=int, default=None, help="text height in cells")
    parser.add_argument(
        "--encoding",
        choices=[
            "halfblock",
            "quadrant",
            "sextant",
            "braille",
            "ascii",
            "kitty",
            "iterm2",
            "sixel",
        ],
        default="halfblock",
    )
    parser.add_argument(
        "--color-depth",
        choices=["auto", "truecolor", "24bit", "256", "16", "mono"],
        default=None,
        help="auto (from NO_COLOR/COLORTERM/TERM) for -, truecolor for files",
    )
    parser.add_argument("--cell-aspect", type=float, default=None)
    parser.add_argument(
        "--tmux",
        action="store_true",
        help="wrap graphics protocols for tmux passthrough",
    )
    parser.add_argument(
        "--no-notes", action="store_true", help="no text under the picture"
    )
    parser.add_argument("--png-compress", type=int, default=None, metavar="0-9")
    parser.add_argument("--cast-frames", type=int, default=None)
    parser.add_argument("--cast-fps", type=float, default=None)
    parser.add_argument("--cast-degrees", type=float, default=None)
    add_render_args(parser)


def snapshot_cmd(args):
    import os

    options = render_options(args)
    to_stdout = args.outfile == "-"
    depth = args.color_depth or ("auto" if to_stdout else "truecolor")
    if depth == "auto":
        depth = _environment_color_depth(_core_module("snapshot"))
    elif depth == "24bit":
        depth = "truecolor"
    graphics = args.encoding in ("kitty", "iterm2", "sixel")
    if graphics and os.environ.get("TMUX") and not args.tmux:
        raise SystemExit(
            f"meshio++: snapshot: inside tmux the {args.encoding} protocol needs --tmux "
            "(and `set -g allow-passthrough on`, tmux 3.3+); a cell encoding works without it"
        )
    text = {"tmux": args.tmux or None, "notes": False if args.no_notes else None}
    if args.cell_aspect is not None:
        text["cell_aspect"] = args.cell_aspect
    text = {k: v for k, v in text.items() if v is not None}

    mesh = read(args.infile, file_format=args.input_format)
    if to_stdout:
        if args.cols is None or args.rows is None:
            size = shutil.get_terminal_size((100, 41))
            note_rows = 1 if args.color_by else 0
            cols = size.columns if args.cols is None else args.cols
            rows = (
                max(1, size.lines - 1 - note_rows) if args.rows is None else args.rows
            )
        else:
            cols, rows = args.cols, args.rows
        if not _enable_terminal_output():
            depth = "mono"
        _write_text(
            render_text(
                mesh,
                cols=cols,
                rows=rows,
                encoding=args.encoding,
                color_depth=depth,
                **text,
                **options,
            )
        )
        return 0
    snap = {
        "png_compress": args.png_compress,
        "cast_frames": args.cast_frames,
        "cast_fps": args.cast_fps,
        "cast_degrees": args.cast_degrees,
    }
    snapshot(
        mesh,
        args.outfile,
        width=args.width,
        height=args.height,
        cols=100 if args.cols is None else args.cols,
        rows=40 if args.rows is None else args.rows,
        encoding=args.encoding,
        color_depth=depth,
        **text,
        **{k: v for k, v in snap.items() if v is not None},
        **options,
    )
    return 0
