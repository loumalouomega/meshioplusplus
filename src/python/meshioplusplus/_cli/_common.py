"""Argument blocks shared by the CLI verbs.

Most verbs read one mesh and many write one: ``INFILE [OUTFILE]`` and the
``--input-format``/``--output-format`` overrides are declared here once, in the
order every verb lists them, so a verb's ``--help`` is unchanged by using it.
"""

from .._helpers import _reader_names, _writer_names


def add_io_args(
    parser,
    out_help=None,
    in_help="mesh file to be read from",
    output_format=None,
):
    """Add ``infile``, ``outfile`` (when ``out_help`` is given) and the format flags.

    ``--output-format`` is added when there is an ``outfile``, unless
    ``output_format`` says otherwise.
    """
    if output_format is None:
        output_format = out_help is not None
    parser.add_argument("infile", type=str, help=in_help)
    if out_help is not None:
        parser.add_argument("outfile", type=str, help=out_help)
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=_reader_names(),
        help="input file format",
        default=None,
    )
    if output_format:
        parser.add_argument(
            "--output-format",
            "-o",
            type=str,
            choices=_writer_names(),
            help="output file format",
            default=None,
        )


# -- value lists ---------------------------------------------------------------


def split_csv(text):
    """``"a, b,,c"`` -> ``["a", "b", "c"]``: stripped, empty entries dropped.

    ``None`` and ``""`` give ``[]``; a caller that tells "not given" apart from
    "given empty" tests for ``None`` first.
    """
    if not text:
        return []
    return [part for part in (p.strip() for p in text.split(",")) if part]


def parse_ints(text, verb, flag, counts=(3,), form="'nx,ny,nz'"):
    """``"4,5,6"`` -> ``[4, 5, 6]``, with ``len`` one of ``counts``.

    ``None`` gives ``None``. A wrong count raises
    ``ValueError("<verb>: <flag> expects <form>")``; a non-integer entry raises
    the ``int`` conversion's own ``ValueError``.
    """
    if text is None:
        return None
    values = [int(x) for x in text.split(",")]
    if len(values) not in counts:
        raise ValueError(f"{verb}: {flag} expects {form}")
    return values


def parse_bounds(text, verb):
    """``"xlo,ylo,zlo,xhi,yhi,zhi"`` -> six floats; ``None`` gives ``None``."""
    if text is None:
        return None
    bounds = [float(x) for x in text.split(",")]
    if len(bounds) != 6:
        raise ValueError(f"{verb}: --bounds expects 'xlo,ylo,zlo,xhi,yhi,zhi'")
    return bounds


# -- option groups -------------------------------------------------------------

_PADDING_RELATIVE_HELP = "grow the box by this fraction of its diagonal on every side"


def add_box_args(
    parser,
    bounds_of="mesh's",
    bounds_help=None,
    padding_relative=0.0,
    padding_relative_help=_PADDING_RELATIVE_HELP,
):
    """Add ``--bounds``, ``--padding`` and ``--padding-relative``, in that order.

    The box of the lattice verbs (``sdf``, ``voxelize``, ``octree``,
    ``octree-mesh``, ``remesh-volume``, ``grid-sample``). ``--bounds``'s help
    names whose bounding box is the default (``bounds_of``) unless
    ``bounds_help`` replaces it whole.
    """
    if bounds_help is None:
        bounds_help = (
            f"'xlo,ylo,zlo,xhi,yhi,zhi'; the {bounds_of} bounding box by default "
            "(negatives need --bounds=)"
        )
    parser.add_argument("--bounds", type=str, default=None, help=bounds_help)
    parser.add_argument(
        "--padding", type=float, default=0.0, help="grow the box by this on every side"
    )
    parser.add_argument(
        "--padding-relative",
        type=float,
        default=padding_relative,
        help=padding_relative_help,
    )


def add_transfer_io_args(parser, out_help):
    """Add ``source``, ``target``, ``outfile`` and the format flags.

    The two-mesh verbs (``interpolate``, ``conservative-interpolate``) read
    both meshes with the one ``--input-format``.
    """
    parser.add_argument("source", type=str, help="mesh file whose data is sampled")
    parser.add_argument("target", type=str, help="mesh file receiving the samples")
    parser.add_argument("outfile", type=str, help=out_help)
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=_reader_names(),
        help="input file format (both source and target)",
        default=None,
    )
    parser.add_argument(
        "--output-format",
        "-o",
        type=str,
        choices=_writer_names(),
        help="output file format",
        default=None,
    )


def read_source_target(args):
    """The ``(source, target)`` meshes :func:`add_transfer_io_args` names."""
    from .._helpers import read

    source = read(args.source, file_format=args.input_format)
    target = read(args.target, file_format=args.input_format)
    return source, target


def add_on_conflict_arg(parser):
    """Add the transfer verbs' ``--on-conflict`` (``error``/``overwrite``/``suffix``)."""
    parser.add_argument(
        "--on-conflict",
        type=str,
        choices=["error", "overwrite", "suffix"],
        default="error",
        help=(
            "what to do when a transferred name already exists on the target: "
            "error (default), overwrite, or suffix (writes to NAME_interp)"
        ),
    )


def add_transfer_quiet_arg(parser):
    """Add the transfer verbs' ``--quiet``/``-q``."""
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="do not print the transfer summary",
    )


def add_decimate_args(parser, unit, help_texts, boundary_flag):
    """The option block of ``decimate`` and ``decimate-volume``.

    ``unit`` is the count flag's noun (``"faces"`` or ``"cells"``, giving
    ``--target-faces``/``--target-cells``); ``help_texts`` maps ``in``,
    ``ratio``, ``count``, ``max_error``, ``placement``, ``features`` and
    ``angle`` to each verb's own wording; ``boundary_flag`` is the
    ``(flags, kwargs)`` of the boundary switch, whose polarity differs between
    the two verbs. The order is the verbs' ``--help`` order.
    """
    add_io_args(parser, "decimated mesh to be written to", in_help=help_texts["in"])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ratio", type=float, default=None, help=help_texts["ratio"])
    group.add_argument(
        f"--target-{unit}", type=int, default=None, help=help_texts["count"]
    )
    group.add_argument(
        "--max-error", type=float, default=None, help=help_texts["max_error"]
    )
    parser.add_argument(
        "--placement",
        choices=["optimal", "midpoint", "endpoint"],
        default="optimal",
        help=help_texts["placement"],
    )
    flags, kwargs = boundary_flag
    parser.add_argument(*flags, **kwargs)
    parser.add_argument(
        "--no-preserve-features",
        dest="preserve_features",
        action="store_false",
        help=help_texts["features"],
    )
    parser.add_argument(
        "--feature-angle", type=float, default=30.0, help=help_texts["angle"]
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="do not print the collapse summary",
    )


def decimate_summary(out, report, noun, removed_key):
    """The one-line summary both decimate verbs print."""
    kept = sum(len(cb.data) for cb in out.cells)
    return (
        f"decimated to {kept} {noun} "
        f"({report[removed_key]} removed, "
        f"{report['points_removed']} points removed, "
        f"{report['collapses_rejected']} collapses rejected, "
        f"max error {report['max_error_applied']:.6g})"
    )
