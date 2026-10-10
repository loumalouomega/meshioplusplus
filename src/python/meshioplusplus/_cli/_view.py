import sys
import warnings

from .._helpers import read, reader_map
from .._viewer import has_viewer, screenshot, view


def _common_args(parser):
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
        "--kind",
        type=str,
        choices=["auto", "surface", "volume", "curve", "points"],
        default="auto",
        help=(
            "what to draw; 'surface' draws only the boundary of a volume mesh, "
            "which is much lighter for a large solid (default: auto)"
        ),
    )
    parser.add_argument(
        "--color-by",
        type=str,
        default=None,
        help="name of a data array to colour by on load",
    )
    parser.add_argument(
        "--name", type=str, default="mesh", help="name shown in the viewer"
    )


def add_args(parser):
    _common_args(parser)
    parser.add_argument(
        "--backend",
        type=str,
        choices=["auto", "polyscope", "browser", "terminal"],
        default="auto",
        help=(
            "auto uses the polyscope desktop window when it is installed and a "
            "display is available; else the browser when a display lets it open; "
            "else the interactive terminal viewer (`tui`) when stdin and stdout "
            "are terminals; else an HTML file (default: auto)"
        ),
    )


def view_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    view(
        mesh,
        kind=args.kind,
        backend=args.backend,
        color_by=args.color_by,
        name=args.name,
    )
    return 0


def add_screenshot_args(parser):
    _common_args(parser)
    parser.add_argument("outfile", type=str, help="PNG file to write")
    parser.add_argument(
        "--size",
        type=int,
        nargs=2,
        metavar=("WIDTH", "HEIGHT"),
        default=(1280, 960),
        help="image size in pixels (default: 1280 960)",
    )
    parser.add_argument(
        "--transparent",
        action="store_true",
        help="transparent rather than opaque background",
    )


def screenshot_cmd(args):
    if not has_viewer():
        # No polyscope: `screenshot` falls back to the software rasterizer
        # (what `snapshot` draws). Say so once, on stderr, with the way to get
        # polyscope's rendering instead.
        print(
            "meshio++: screenshot: polyscope is not installed (`pip install "
            "meshioplusplus[viewer]`); drawing with the software rasterizer",
            file=sys.stderr,
        )
    mesh = read(args.infile, file_format=args.input_format)
    with warnings.catch_warnings():
        # The fallback's warning repeats the line printed above.
        warnings.simplefilter("ignore", RuntimeWarning)
        screenshot(
            mesh,
            args.outfile,
            kind=args.kind,
            color_by=args.color_by,
            name=args.name,
            size=tuple(args.size),
            transparent=args.transparent,
        )
    return 0
