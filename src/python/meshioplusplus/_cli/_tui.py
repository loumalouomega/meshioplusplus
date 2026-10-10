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
    parser.add_argument("infile", type=str, help="mesh file to be read from")
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
    mesh = read(args.infile, file_format=args.input_format)
    replay = None
    if args.replay is not None:
        with open(args.replay, "rb") as handle:
            replay = handle.read()
    result = tui(
        mesh,
        encoding=args.encoding,
        color_depth=depth,
        title=args.infile.replace("\\", "/").rsplit("/", 1)[-1],
        replay=replay,
        cols=args.cols,
        rows=args.rows,
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
