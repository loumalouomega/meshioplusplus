"""``meshioplusplus tui``: the interactive terminal viewer from the command line.

The flags are the native CLI's (``src/cpp/cli/main.cpp``, ``cmd_tui``): the
render flags of ``snapshot`` plus the encoding ones, so the two verbs take the
same command lines.
"""

import sys

from .._helpers import read, reader_map
from .._tui import tui
from ._snapshot import _write_text, add_render_args, render_options


def _music_key(text):
    """A note name (C, C#, Db, ...) or a semitone above C, as the native CLI reads it."""
    keys = {
        "C": 0, "C#": 1, "DB": 1, "D": 2, "D#": 3, "EB": 3, "E": 4, "F": 5, "F#": 6,
        "GB": 6, "G": 7, "G#": 8, "AB": 8, "A": 9, "A#": 10, "BB": 10, "B": 11,
    }  # fmt: skip
    if text.upper() in keys:
        return keys[text.upper()]
    try:
        value = int(text)
    except ValueError:
        value = -1
    if not 0 <= value <= 11:
        raise SystemExit(
            f"meshio++: --music-key expects a note name (C, C#, Db, ..., B) or 0 to 11, not '{text}'"
        )
    return value


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
    parser.add_argument(
        "--music",
        action="store_true",
        help="play a generated synthwave loop on an external player (silent in CI "
        "and over SSH); m mutes it, < > change the volume",
    )
    parser.add_argument(
        "--music-out",
        type=str,
        default=None,
        metavar="FILE.wav",
        help="write the loop to a WAV file for your own player",
    )
    parser.add_argument("--volume", type=float, default=None, help="0.05 to 0.9 (0.3)")
    parser.add_argument("--tempo", type=float, default=None, help="60 to 140 BPM (100)")
    parser.add_argument("--music-seed", type=int, default=None)
    parser.add_argument(
        "--music-key", type=str, default=None, help="C, C#, Db, ... B, or 0 to 11 (A)"
    )
    parser.add_argument(
        "--music-over-ssh",
        action="store_true",
        help="play even in an SSH session (the sound comes out of the machine "
        "the process runs on)",
    )
    parser.add_argument(
        "--pulse", action="store_true", help="with --theme: step the grid on every beat"
    )
    parser.add_argument(
        "--reduced-motion",
        action="store_true",
        help="no beat pulse (also: REDUCED_MOTION in the environment)",
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
    for flag, key in (
        ("volume", "volume"),
        ("tempo", "tempo"),
        ("music_seed", "music_seed"),
    ):
        if getattr(args, flag) is not None:
            extra[key] = getattr(args, flag)
    if args.music_key is not None:
        extra["music_key"] = _music_key(args.music_key)
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
        music=args.music,
        music_out=args.music_out,
        music_over_ssh=args.music_over_ssh,
        pulse=args.pulse,
        reduced_motion=args.reduced_motion,
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
