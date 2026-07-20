"""The nested ``meshioplusplus data <verb>`` command group.

The nine verbs are thin front-ends over the five data operations:
``rename``/``drop``/``keep`` call ``data_manage``, ``to-cell``/``to-point``
call the two averaging entry points, ``clamp``/``normalize`` call
``data_condition`` with a preset mode, and ``calc``/``info`` map one-to-one.

This is the repository's first two-level subcommand. It needs no change to
``_main.py``'s dispatch because ``set_defaults(func=...)`` on the *inner*
parser overrides the outer one, so the existing ``args.func(args)`` call
already reaches the right handler.
"""

import json

from .._data_average import cell_data_to_point_data, point_data_to_cell_data
from .._data_calc import data_calc
from .._data_condition import data_condition
from .._data_info import data_info
from .._data_manage import data_drop, data_keep, data_rename
from .._helpers import _writer_map, read, reader_map, write

_LOCATION_FLAGS = ("point", "cell", "field")


def _add_io_args(parser, with_output=True):
    """The shared ``INFILE [OUTFILE]`` + format-override arguments."""
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    if with_output:
        parser.add_argument("outfile", type=str, help="mesh file to be written to")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    if with_output:
        parser.add_argument(
            "--output-format",
            "-o",
            type=str,
            choices=sorted(list(_writer_map.keys())),
            help="output file format",
            default=None,
        )


def _split_names(value):
    """``"a,b"`` -> ``["a", "b"]``, ignoring empty entries."""
    if not value:
        return []
    return [part for part in (p.strip() for p in value.split(",")) if part]


def _split_rename(value):
    """``"OLD:NEW"`` -> ``("OLD", "NEW")``, splitting on the LAST colon.

    Data names routinely contain colons (``gmsh:physical``), so splitting on the
    last one lets ``--point gmsh:physical:tag`` rename ``gmsh:physical`` to
    ``tag``. The same rule is implemented in the native CLI and documented in
    ``doc/cli.md``.
    """
    if ":" not in value:
        raise ValueError(
            f"data rename: expected OLD:NEW, got '{value}' (no ':' separator)"
        )
    old, _, new = value.rpartition(":")
    if not old or not new:
        raise ValueError(f"data rename: expected OLD:NEW, got '{value}'")
    return old, new


def _split_assignment(value):
    """``"name = expr"`` -> ``("name", "expr")``, splitting on the FIRST ``=``."""
    if "=" not in value:
        raise ValueError(
            f"data calc: expected 'NAME = EXPRESSION', got '{value}' (no '=')"
        )
    name, _, expr = value.partition("=")
    name = name.strip()
    expr = expr.strip()
    if not name or not expr:
        raise ValueError(f"data calc: expected 'NAME = EXPRESSION', got '{value}'")
    return name, expr


# --- data info -------------------------------------------------------------


def add_info_args(parser):
    _add_io_args(parser, with_output=False)
    parser.add_argument(
        "--json", action="store_true", help="emit the summary as JSON"
    )


def info_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    arrays = data_info(mesh)
    if args.json:
        print(json.dumps(arrays, indent=2))
        return 0
    print("<meshio++ data summary>")
    if not arrays:
        print("  (the mesh carries no data arrays)")
        return 0
    header = f"  {'location':<11} {'name':<20} {'dtype':<6} {'comp':>4} " \
             f"{'entries':>8} {'min':>12} {'max':>12} {'mean':>12} {'nan':>5} {'inf':>5}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for a in arrays:
        print(
            f"  {a['location']:<11} {a['name']:<20} {a['dtype']:<6} "
            f"{a['num_components']:>4} {a['num_entries']:>8} "
            f"{a['min']:>12.6g} {a['max']:>12.6g} {a['mean']:>12.6g} "
            f"{a['num_nan']:>5} {a['num_inf']:>5}"
        )
        if a["inconsistent_blocks"]:
            print(f"  {'':<11} {'':<20} (warning: cell blocks disagree in components)")
    return 0


# --- data rename / drop / keep --------------------------------------------


def add_rename_args(parser):
    _add_io_args(parser)
    for loc in _LOCATION_FLAGS:
        parser.add_argument(
            f"--{loc}",
            action="append",
            default=[],
            metavar="OLD:NEW",
            help=f"rename a {loc}_data array (repeatable; split on the last ':')",
        )


def rename_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    total = 0
    for loc in _LOCATION_FLAGS:
        for spec in getattr(args, loc):
            old, new = _split_rename(spec)
            mesh = data_rename(mesh, loc, old, new)
            print(f"renamed {loc}_data '{old}' -> '{new}'")
            total += 1
    if total == 0:
        print("data rename: nothing to do (pass --point/--cell/--field OLD:NEW)")
    write(args.outfile, mesh, file_format=args.output_format)
    return 0


def add_drop_args(parser):
    _add_io_args(parser)
    for loc in _LOCATION_FLAGS:
        parser.add_argument(
            f"--{loc}",
            type=str,
            default=None,
            metavar="A,B",
            help=f"comma-separated {loc}_data array names to drop",
        )
    parser.add_argument(
        "--ignore-missing",
        action="store_true",
        help="skip names that do not exist instead of failing",
    )


def drop_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    for loc in _LOCATION_FLAGS:
        names = _split_names(getattr(args, loc))
        if names:
            mesh = data_drop(mesh, loc, names, ignore_missing=args.ignore_missing)
            print(f"dropped {loc}_data: {', '.join(names)}")
    write(args.outfile, mesh, file_format=args.output_format)
    return 0


def add_keep_args(parser):
    _add_io_args(parser)
    for loc in _LOCATION_FLAGS:
        parser.add_argument(
            f"--{loc}",
            type=str,
            default=None,
            metavar="A,B",
            help=f"comma-separated {loc}_data array names to keep "
            "(the location is left untouched if this is omitted)",
        )
    parser.add_argument(
        "--ignore-missing",
        action="store_true",
        help="skip names that do not exist instead of failing",
    )


def keep_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    for loc in _LOCATION_FLAGS:
        value = getattr(args, loc)
        # A location the user did not name at all is left alone; naming it with
        # an empty list means "keep nothing there".
        if value is None:
            continue
        names = _split_names(value)
        mesh = data_keep(mesh, loc, names, ignore_missing=args.ignore_missing)
        print(f"kept {loc}_data: {', '.join(names) if names else '(nothing)'}")
    write(args.outfile, mesh, file_format=args.output_format)
    return 0


# --- data to-cell / to-point ----------------------------------------------


def add_to_cell_args(parser):
    _add_io_args(parser)
    parser.add_argument(
        "--keys",
        type=str,
        default=None,
        metavar="A,B",
        help="comma-separated point_data names (default: all)",
    )
    parser.add_argument(
        "--target-suffix",
        type=str,
        default="",
        help="append this to each output name (default: same name)",
    )


def to_cell_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    keys = _split_names(args.keys)
    out = point_data_to_cell_data(mesh, keys=keys, suffix=args.target_suffix)
    moved = keys if keys else sorted(mesh.point_data.keys())
    print(f"averaged {len(moved)} point_data array(s) onto the cells")
    write(args.outfile, out, file_format=args.output_format)
    return 0


def add_to_point_args(parser):
    _add_io_args(parser)
    parser.add_argument(
        "--keys",
        type=str,
        default=None,
        metavar="A,B",
        help="comma-separated cell_data names (default: all)",
    )
    parser.add_argument(
        "--weighted",
        action="store_true",
        help="weight by cell measure (area/volume) instead of counting cells equally",
    )
    parser.add_argument(
        "--target-suffix",
        type=str,
        default="",
        help="append this to each output name (default: same name)",
    )


def to_point_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    keys = _split_names(args.keys)
    out = cell_data_to_point_data(
        mesh, keys=keys, weighted=args.weighted, suffix=args.target_suffix
    )
    moved = keys if keys else sorted(mesh.cell_data.keys())
    how = "measure-weighted" if args.weighted else "unweighted"
    print(f"averaged {len(moved)} cell_data array(s) onto the points ({how})")
    write(args.outfile, out, file_format=args.output_format)
    return 0


# --- data calc -------------------------------------------------------------


def add_calc_args(parser):
    _add_io_args(parser)
    for loc in _LOCATION_FLAGS:
        parser.add_argument(
            f"--{loc}",
            action="append",
            default=[],
            metavar="NAME = EXPR",
            help=f"derive a {loc}_data array (repeatable)",
        )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="allow replacing an array that already exists",
    )


def calc_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    total = 0
    for loc in _LOCATION_FLAGS:
        for spec in getattr(args, loc):
            name, expr = _split_assignment(spec)
            mesh = data_calc(
                mesh, expr, location=loc, output=name, overwrite=args.overwrite
            )
            print(f"computed {loc}_data '{name}' = {expr}")
            total += 1
    if total == 0:
        print("data calc: nothing to do (pass --point/--cell/--field 'NAME = EXPR')")
    write(args.outfile, mesh, file_format=args.output_format)
    return 0


# --- data clamp / normalize ------------------------------------------------


def _add_condition_common(parser):
    _add_io_args(parser)
    for loc in _LOCATION_FLAGS:
        parser.add_argument(
            f"--{loc}",
            type=str,
            default=None,
            metavar="A,B",
            help=f"comma-separated {loc}_data names (default: all at that location)",
        )
    parser.add_argument(
        "--magnitude",
        action="store_true",
        help="condition by row magnitude instead of per component",
    )
    parser.add_argument(
        "--nan",
        type=str,
        choices=["ignore", "replace", "fail"],
        default="ignore",
        help="what reaches the output for non-finite values (default: ignore)",
    )
    parser.add_argument(
        "--nan-value",
        type=float,
        default=0.0,
        help="replacement used with --nan replace",
    )
    parser.add_argument(
        "--suffix",
        type=str,
        default="",
        help="store the result as NAME+SUFFIX instead of replacing in place",
    )


def _selected_locations(args):
    """The (location, names) pairs the user actually asked for."""
    out = []
    for loc in _LOCATION_FLAGS:
        value = getattr(args, loc)
        if value is not None:
            out.append((loc, _split_names(value)))
    return out


def add_clamp_args(parser):
    _add_condition_common(parser)
    parser.add_argument(
        "--min", dest="lo", type=float, required=True, help="lower bound"
    )
    parser.add_argument(
        "--max", dest="hi", type=float, required=True, help="upper bound"
    )


def clamp_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    selected = _selected_locations(args)
    if not selected:
        print("data clamp: nothing to do (pass --point/--cell/--field NAME)")
    for loc, names in selected:
        mesh = data_condition(
            mesh,
            location=loc,
            keys=names,
            mode="clamp",
            scope="magnitude" if args.magnitude else "component",
            lo=args.lo,
            hi=args.hi,
            nan_policy=args.nan,
            nan_replacement=args.nan_value,
            suffix=args.suffix,
        )
        label = ", ".join(names) if names else "(all)"
        print(f"clamped {loc}_data {label} to [{args.lo:g}, {args.hi:g}]")
    write(args.outfile, mesh, file_format=args.output_format)
    return 0


def add_normalize_args(parser):
    _add_condition_common(parser)
    parser.add_argument(
        "--to",
        type=str,
        default="0,1",
        metavar="LO,HI",
        help="target range (default: 0,1)",
    )
    parser.add_argument(
        "--zero-mean",
        action="store_true",
        help="standardize to zero mean and unit standard deviation instead",
    )


def normalize_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    parts = _split_names(args.to)
    if not args.zero_mean and len(parts) != 2:
        raise ValueError(f"data normalize: --to expects LO,HI, got '{args.to}'")
    lo, hi = (0.0, 1.0) if args.zero_mean else (float(parts[0]), float(parts[1]))
    mode = "standardize" if args.zero_mean else "normalize"

    selected = _selected_locations(args)
    if not selected:
        print("data normalize: nothing to do (pass --point/--cell/--field NAME)")
    for loc, names in selected:
        mesh = data_condition(
            mesh,
            location=loc,
            keys=names,
            mode=mode,
            scope="magnitude" if args.magnitude else "component",
            lo=lo,
            hi=hi,
            nan_policy=args.nan,
            nan_replacement=args.nan_value,
            suffix=args.suffix,
        )
        label = ", ".join(names) if names else "(all)"
        target = "zero mean / unit std" if args.zero_mean else f"[{lo:g}, {hi:g}]"
        print(f"normalized {loc}_data {label} to {target}")
    write(args.outfile, mesh, file_format=args.output_format)
    return 0


# --- group wiring ----------------------------------------------------------

_VERBS = (
    ("info", "Summarize every data array (dtype, shape, min/max/mean, NaN/inf)",
     add_info_args, info_cmd),
    ("rename", "Rename data arrays", add_rename_args, rename_cmd),
    ("drop", "Drop data arrays by name", add_drop_args, drop_cmd),
    ("keep", "Keep only the named data arrays", add_keep_args, keep_cmd),
    ("to-cell", "Average point_data onto the cells", add_to_cell_args, to_cell_cmd),
    ("to-point", "Average cell_data onto the points", add_to_point_args, to_point_cmd),
    ("calc", "Derive a new array from an elementwise expression",
     add_calc_args, calc_cmd),
    ("clamp", "Clamp values into [min, max]", add_clamp_args, clamp_cmd),
    ("normalize", "Rescale values to a target range (or zero mean / unit std)",
     add_normalize_args, normalize_cmd),
)


def add_args(parser):
    sub = parser.add_subparsers(title="data subcommands", dest="data_command")
    sub.required = True
    for name, help_text, add, func in _VERBS:
        p = sub.add_parser(name, help=help_text)
        add(p)
        p.set_defaults(func=func)


def data_cmd(args):
    """Reached only if argparse somehow admits a bare ``data`` invocation."""
    print("error: 'data' requires a subcommand: " + ", ".join(v[0] for v in _VERBS))
    return 2
