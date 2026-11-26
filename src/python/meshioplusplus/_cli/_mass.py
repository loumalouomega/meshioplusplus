"""The ``mass`` verb: mass, centre of mass and inertia of a body
(``doc/mass_properties.md``)."""

from .._helpers import read
from .._mass_properties import mass_properties
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(parser)
    parser.add_argument(
        "--density", type=float, default=1.0, help="density of every cell (default 1)"
    )
    parser.add_argument(
        "--density-array",
        default=None,
        metavar="NAME",
        help="scalar cell_data array of per-cell densities (replaces --density)",
    )
    parser.add_argument(
        "--region", default=None, help="restrict the total to this cell region"
    )
    parser.add_argument(
        "--thickness",
        type=float,
        default=0.0,
        help="treat surface cells as a shell of this thickness (open surfaces allowed)",
    )
    parser.add_argument(
        "--about",
        choices=("centroid", "origin"),
        default="centroid",
        help="reference point of the inertia tensor",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def _fmt(v):
    return "[" + ", ".join(f"{x:.6g}" for x in v) + "]"


def _print_entry(e, kind, indent="  "):
    label = "area" if kind == "shell" else "volume"
    print(f"{indent}cells:              {e['num_cells']}")
    print(f"{indent}{label + ':':<19}{e['measure']:.6g}")
    print(f"{indent}mass:               {e['mass']:.6g}")
    print(f"{indent}centre of mass:     {_fmt(e['centre'])}")
    print(f"{indent}inertia xx yy zz xy yz zx: {_fmt(e['inertia'])}")
    print(f"{indent}principal moments:  {_fmt(e['principal_moments'])}")
    print(f"{indent}gyration radii:     {_fmt(e['gyration_radii'])}")


def mass_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    report = mass_properties(
        mesh,
        density=args.density,
        density_array=args.density_array,
        region=args.region,
        thickness=args.thickness,
        about=args.about,
    )
    if args.json:
        emit_json(report)
        return
    print(f"<meshio++ mass properties: {report['kind']}>")
    _print_entry(report["total"], report["kind"])
    for entry in report["regions"]:
        print(f"  region '{entry['name']}':")
        _print_entry(entry, report["kind"], indent="    ")
    for w in report["warnings"]:
        print(f"  warning: {w}")
