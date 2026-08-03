import pathlib

import numpy as np

from .._colormap import NAMES as CMAP_NAMES
from .._helpers import _writer_map, read, reader_map, write

# Colouring reaches the writer as keyword arguments, and only the SVG and TikZ
# writers accept them -- every other writer would raise on the unexpected kwarg,
# so the flags are validated against the resolved output format up front.
_COLOR_FORMATS = ("svg", "tikz")


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
        "--output-format",
        "-o",
        type=str,
        choices=sorted(list(_writer_map.keys())),
        help="output file format",
        default=None,
    )
    parser.add_argument(
        "--ascii",
        "-a",
        action="store_true",
        help="write in ASCII format variant (where applicable, default: binary)",
    )
    parser.add_argument("outfile", type=str, help="mesh file to be written to")
    parser.add_argument(
        "--float-format",
        "-f",
        type=str,
        help="float format used in output ASCII files (default: .16e)",
    )
    parser.add_argument(
        "--sets-to-int-data",
        "-s",
        action="store_true",
        help="if possible, convert sets to integer data (useful if the output type does not support sets)",
    )
    parser.add_argument(
        "--int-data-to-sets",
        "-d",
        action="store_true",
        help="if possible, convert integer data to sets (useful if the output type does not support integer data)",
    )
    parser.add_argument(
        "--points-only",
        action="store_true",
        help=(
            "read geometry only, dropping every data array. Readers with a "
            "native selective path skip the arrays outright rather than "
            "reading and discarding them."
        ),
    )
    parser.add_argument(
        "--arrays",
        type=str,
        default=None,
        help=(
            "comma-separated data-array names to keep; everything else is "
            "skipped. Pass an empty string to keep none. Names absent from the "
            "file are ignored."
        ),
    )
    parser.add_argument(
        "--time-step",
        type=int,
        default=0,
        metavar="N",
        help=(
            "which step of a multi-step file to read: 0 (default) is the "
            "first, negative counts from the end (-1 = last). Out of range is "
            "an error naming the available count. Honoured by formats "
            "carrying a time series (currently exodus); a negative value "
            "needs the --time-step=-1 form."
        ),
    )
    seq = parser.add_argument_group(
        "transient sequences (multi-file datasets)",
        "Treat a set of files, or the steps inside one file, as one ordered "
        "dataset. See doc/sequences.md.",
    )
    seq.add_argument(
        "--input",
        action="append",
        default=None,
        metavar="FILE",
        help=(
            "an EXTRA input file, appended after the positional infile; repeat "
            "for each. Use this when your shell has already expanded a glob: "
            "`convert a.vtu --input b.vtu --input c.vtu out.xdmf`. Quoting the "
            "pattern instead (`convert 'out_*.vtu' out.xdmf`) lets the CLI "
            "expand it itself and is usually what you want."
        ),
    )
    seq.add_argument(
        "--times",
        type=str,
        default=None,
        metavar="T1,T2,...",
        help="comma-separated time value per step; the count must match.",
    )
    seq.add_argument(
        "--time-from",
        type=str,
        choices=("auto", "file", "filename", "index"),
        default="auto",
        help=(
            "where each step's time comes from: auto (the documented "
            "precedence), file, filename (its last digit run) or index."
        ),
    )
    seq.add_argument(
        "--sequence",
        action="store_true",
        help="force sequence handling even when nothing in the paths suggests it.",
    )
    seq.add_argument(
        "--no-sequence",
        action="store_true",
        help=(
            "force single-file handling -- the escape hatch for a filename that "
            "genuinely contains '*' or '{step}', which is legal on POSIX."
        ),
    )
    color = parser.add_argument_group(
        "data-driven colouring (svg/tikz output only)",
        "Colour the drawn faces by a data array instead of a flat fill.",
    )
    color.add_argument(
        "--color-by",
        type=str,
        default=None,
        metavar="NAME",
        help=(
            "point_data or cell_data array to colour by. Point data uses the "
            "mean of a face's corner values; cell data its owning cell's value."
        ),
    )
    color.add_argument(
        "--component",
        type=int,
        default=None,
        metavar="I",
        help="component of a multi-component array (default: its magnitude)",
    )
    color.add_argument(
        "--cmap",
        type=str,
        default="viridis",
        choices=CMAP_NAMES,
        help="colormap to map the value through (default: viridis)",
    )
    color.add_argument(
        "--vmin",
        type=float,
        default=None,
        metavar="V",
        help="low end of the colour range (default: the drawn faces' minimum)",
    )
    color.add_argument(
        "--vmax",
        type=float,
        default=None,
        metavar="V",
        help="high end of the colour range (default: the drawn faces' maximum)",
    )
    color.add_argument(
        "--nan-color",
        type=str,
        default=None,
        metavar="C",
        help="colour for faces whose value is NaN or infinite",
    )
    color.add_argument(
        "--colorbar",
        action="store_true",
        help="append a gradient bar with min/max labels",
    )


def _color_kwargs(args):
    """Validate the colouring flags and turn them into writer kwargs.

    Returns an empty dict when no colouring was requested. The modifier flags
    are meaningless without ``--color-by``, and the whole group is meaningless
    for any writer but SVG/TikZ, so both are errors rather than silent no-ops.
    """
    modifiers = {
        "--component": args.component is not None,
        "--cmap": args.cmap != "viridis",
        "--vmin": args.vmin is not None,
        "--vmax": args.vmax is not None,
        "--nan-color": args.nan_color is not None,
        "--colorbar": args.colorbar,
    }
    if args.color_by is None:
        used = [flag for flag, given in modifiers.items() if given]
        if used:
            raise ValueError(f"{', '.join(used)} require(s) --color-by")
        return {}

    fmt = args.output_format or pathlib.Path(args.outfile).suffix.lstrip(".").lower()
    if fmt not in _COLOR_FORMATS:
        raise ValueError(
            f"--color-by is only supported for {'/'.join(_COLOR_FORMATS)} output, "
            f"not '{fmt}'"
        )
    if args.ascii:
        raise ValueError(f"--ascii has no meaning for {fmt} output")

    kwargs = {
        "color_by": args.color_by,
        "component": args.component,
        "cmap": args.cmap,
        "vmin": args.vmin,
        "vmax": args.vmax,
        "colorbar": args.colorbar,
    }
    # Leave nan_color unset so each writer keeps its own format-native default.
    if args.nan_color is not None:
        kwargs["nan_color"] = args.nan_color
    return kwargs


def _wants_sequence(args):
    """Whether this invocation describes a transient sequence.

    Inferred from the paths -- extra ``--input``s, a glob in the input, or a
    ``{step}``/``{index}`` token in the output -- with ``--sequence`` and
    ``--no-sequence`` as explicit overrides. ``--no-sequence`` exists because a
    filename may legitimately contain ``*`` on POSIX.
    """
    from .._sequence import pattern_has_token

    if args.no_sequence and args.sequence:
        raise ValueError("--sequence and --no-sequence are mutually exclusive")
    if args.no_sequence:
        return False
    if args.sequence:
        return True
    if args.input:
        return True
    if "*" in args.infile or "?" in args.infile:
        return True
    if pattern_has_token(args.outfile):
        return True
    # A multi-step input aimed at a single-step output must not quietly become
    # step 0 -- route it here so the sequence driver refuses by name, pointing
    # at '{step}' and at --time-step. An explicit --time-step IS a deliberate
    # single-step selection, so it opts out. The check is skipped entirely for
    # the formats that cannot carry time (see _sequence._TIME_CAPABLE_READERS),
    # which is what keeps it free for almost every convert.
    from .._sequence import num_steps

    if not args.time_step:
        return num_steps(args.infile, args.input_format) > 1
    return False


def _convert_sequence(args, arrays, write_kwargs):
    """The transient half of `convert`: fan-in, fan-out or N->N."""
    from .._sequence import (
        pattern_has_token,
        read_sequence,
        sequence_entries,
        write_sequence,
    )

    # A pre-expanded argv (`--input a --input b`) and a quoted pattern
    # ('out_*.vtu') must reach exactly the same code, so both become the same
    # `source` here.
    source = [args.infile, *args.input] if args.input else args.infile
    times = None
    if args.times is not None:
        times = [float(t) for t in args.times.split(",") if t]

    entries = sequence_entries(
        source,
        file_format=args.input_format,
        times=times,
        time_from=args.time_from,
    )
    steps = read_sequence(
        source,
        file_format=args.input_format,
        times=times,
        time_from=args.time_from,
        points_only=args.points_only,
        arrays=arrays,
    )
    written = write_sequence(
        args.outfile, steps, file_format=args.output_format, **write_kwargs
    )
    what = "fan-out" if pattern_has_token(args.outfile) else "fan-in"
    if not pattern_has_token(args.outfile) and len(written) == len(entries):
        what = "sequence"
    print(f"{what}: {len(entries)} step(s) -> {len(written)} file(s)")
    return 0


def convert(args):
    if args.points_only and args.arrays is not None:
        raise ValueError("--points-only and --arrays are mutually exclusive")

    # Validate before reading: a bad flag combination should fail immediately,
    # not after loading a large mesh.
    color_kwargs = _color_kwargs(args)

    arrays = None
    if args.arrays is not None:
        arrays = [name for name in args.arrays.split(",") if name]

    if _wants_sequence(args):
        if args.sets_to_int_data or args.int_data_to_sets:
            raise ValueError(
                "-s/-d are not available for a sequence; convert one step at a "
                "time, or use the `pipeline` verb's Operations chain"
            )
        write_kwargs = dict(color_kwargs)
        if args.float_format is not None:
            write_kwargs["float_fmt"] = args.float_format
        if args.ascii:
            write_kwargs["binary"] = False
        return _convert_sequence(args, arrays, write_kwargs)

    # read mesh data
    mesh = read(
        args.infile,
        file_format=args.input_format,
        points_only=args.points_only,
        arrays=arrays,
        time_step=args.time_step,
    )

    # Some converters (like VTK) require `points` to be contiguous.
    mesh.points = np.ascontiguousarray(mesh.points)

    if (args.points_only or arrays is not None) and (
        args.sets_to_int_data or args.int_data_to_sets
    ):
        raise ValueError(
            "--points-only/--arrays cannot be combined with -s/-d: the set "
            "conversions operate on the very data arrays that were skipped"
        )

    if args.sets_to_int_data:
        mesh.point_sets_to_data()
        mesh.cell_sets_to_data()

    if args.int_data_to_sets:
        # Snapshot the keys: *_data_to_sets(key) deletes that key from the dict,
        # so iterating the live view raises "dictionary changed size".
        for key in list(mesh.point_data):
            mesh.point_data_to_sets(key)
        for key in list(mesh.cell_data):
            mesh.cell_data_to_sets(key)

    # write it out
    kwargs = {"file_format": args.output_format}
    if args.float_format is not None:
        kwargs["float_fmt"] = args.float_format
    if args.ascii:
        kwargs["binary"] = False
    kwargs.update(color_kwargs)

    write(args.outfile, mesh, **kwargs)
