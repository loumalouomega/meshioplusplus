"""``meshioplusplus tui``: the interactive terminal viewer from the command line.

The flags are the native CLI's (``src/cpp/cli/main.cpp``, ``cmd_tui``): the
render flags of ``snapshot`` plus the encoding ones, so the two verbs take the
same command lines.
"""

import sys

from .._helpers import read, reader_map
from .._tui import tui
from ._snapshot import _write_text, add_render_args, render_options


def add_args(parser):
    parser.add_argument(
        "infile",
        type=str,
        nargs="+",
        help="mesh file to be read from; several files, or a quoted glob "
        "(`'out_*.vtu'`), make a time series",
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
        default="auto",
    )
    parser.add_argument("--cell-aspect", type=float, default=None)
    parser.add_argument(
        "--tmux",
        action="store_true",
        help="wrap graphics protocols for tmux passthrough",
    )
    parser.add_argument(
        "--replay",
        type=str,
        default=None,
        metavar="FILE",
        help="play a recorded input stream on a --cols x --rows screen and print "
        "what the viewer writes (for tests and documentation figures)",
    )
    parser.add_argument(
        "--compare", type=str, default=None, help="a second mesh drawn beside the first"
    )
    parser.add_argument(
        "--diff",
        action="store_true",
        help="with --compare: draw |B - A| of the --color-by point array",
    )
    parser.add_argument(
        "--separate-ranges",
        action="store_true",
        help="with --compare: each mesh its own colour range",
    )
    parser.add_argument(
        "--follow",
        action="store_true",
        help="watch the series (or the file) and show a new step once it is complete",
    )
    parser.add_argument(
        "--follow-interval",
        type=int,
        default=None,
        metavar="MS",
        help="with --follow: how often to look (default 500)",
    )
    parser.add_argument(
        "--settle",
        type=int,
        default=None,
        metavar="MS",
        help="with --follow: wait until the file has not changed for this long (default 300)",
    )
    parser.add_argument(
        "--fps",
        type=float,
        default=None,
        help="steps per second when playing (default 4)",
    )
    parser.add_argument(
        "--session",
        type=str,
        default=None,
        metavar="FILE",
        help="a session file: read at the start if it exists, written at the end",
    )
    parser.add_argument("--cols", type=int, default=100, help="replay screen width")
    parser.add_argument("--rows", type=int, default=40, help="replay screen height")
    add_render_args(parser)


def tui_cmd(args):
    options = render_options(args)
    depth = None if args.color_depth == "auto" else args.color_depth
    if depth == "24bit":
        depth = "truecolor"
    if args.cell_aspect is not None:
        options["cell_aspect"] = args.cell_aspect
    if args.tmux:
        options["tmux"] = True
    if args.diff and not args.compare:
        raise SystemExit("meshio++: --diff requires --compare FILE")
    if args.follow_interval is not None and not args.follow:
        raise SystemExit("meshio++: --follow-interval requires --follow")
    infiles = args.infile
    glob = len(infiles) == 1 and any(ch in infiles[0] for ch in "*?")
    series = None
    mesh = None
    if glob:
        series = infiles[0]
    elif len(infiles) > 1 or args.follow:
        series = infiles
    else:
        mesh = read(infiles[0], file_format=args.input_format)
    replay = None
    if args.replay is not None:
        with open(args.replay, "rb") as handle:
            replay = handle.read()
    extra = {}
    if args.settle is not None:
        extra["settle"] = args.settle
    if args.fps is not None:
        extra["fps"] = args.fps
    result = tui(
        mesh,
        encoding=args.encoding,
        color_depth=depth,
        title=infiles[0].replace("\\", "/").rsplit("/", 1)[-1],
        replay=replay,
        cols=args.cols,
        rows=args.rows,
        compare=args.compare,
        diff=args.diff,
        shared_range=not args.separate_ranges,
        series=series,
        input_format=args.input_format,
        follow=(args.follow_interval or 500) if args.follow else False,
        session=args.session,
        **extra,
        **options,
    )
    if replay is not None:
        buffer = getattr(sys.stdout, "buffer", None)
        if buffer is None:
            _write_text(result["output"].decode("utf-8"))
        else:
            sys.stdout.flush()
            buffer.write(result["output"])
            buffer.flush()
    return result["exit"]
