from .. import _primitives
from .._helpers import _writer_map, write

#: shape name on the command line -> (function, value options, choice options)
_SHAPES = {
    "rectangle": (
        _primitives.rectangle,
        {"dims": "nx,ny", "size": "lx,ly", "origin": "x,y[,z]"},
        {
            "cell_type": ("quad", ["quad", "triangle"]),
            "split": ("zigzag", ["zigzag", "up", "down", "center"]),
        },
    ),
    "box": (
        _primitives.box,
        {"dims": "nx,ny,nz", "size": "lx,ly,lz", "origin": "x,y,z"},
        {"cell_type": ("hexahedron", ["hexahedron", "tetra", "wedge"])},
    ),
    "disk": (
        _primitives.disk,
        {"radius": None, "n": None, "center": "x,y,z"},
        {"cell_type": ("quad", ["quad", "triangle"])},
    ),
    "uv-sphere": (
        _primitives.uv_sphere,
        {"radius": None, "nu": None, "nv": None, "center": "x,y,z"},
        {"cell_type": ("triangle", ["triangle", "quad"])},
    ),
    "icosa-sphere": (
        _primitives.icosa_sphere,
        {"radius": None, "subdivisions": None, "center": "x,y,z"},
        {},
    ),
    "cylinder": (
        _primitives.cylinder,
        {"radius": None, "height": None, "nxy": None, "nz": None, "center": "x,y,z"},
        {"cell_type": ("hexahedron", ["hexahedron", "wedge"])},
    ),
    "reference-cells": (_primitives.reference_cells, {}, {}),
}

_INTEGERS = {"n", "nu", "nv", "nxy", "nz", "subdivisions"}
_FLOATS = {"radius", "height"}
_VECTOR_KINDS = {"dims": int, "size": float, "origin": float, "center": float}


def add_args(parser):
    sub = parser.add_subparsers(dest="shape", metavar="SHAPE", required=True)
    for name, (_, values, choices) in _SHAPES.items():
        p = sub.add_parser(name, help=_doc_line(name))
        p.add_argument("outfile", type=str, help="mesh to be written to")
        p.add_argument(
            "--output-format",
            "-o",
            type=str,
            choices=sorted(list(_writer_map.keys())),
            help="output file format",
            default=None,
        )
        for key, metavar in values.items():
            p.add_argument(
                "--" + key.replace("_", "-"),
                type=str,
                default=None,
                metavar=(metavar or key.upper()),
                help=f"{key} (the library default when omitted; negatives need "
                f"--{key.replace('_', '-')}=)",
            )
        for key, (default, options) in choices.items():
            p.add_argument(
                "--" + key.replace("_", "-"),
                type=str,
                default=default,
                choices=options,
                help=f"cell type (default {default})",
            )
        if name == "reference-cells":
            p.add_argument(
                "--types",
                type=str,
                default=None,
                help="comma-separated cell types (every supported type when omitted)",
            )
            p.add_argument(
                "--pairs",
                action="store_true",
                help="two cells sharing a facet instead of one (linear types only)",
            )
        else:
            p.add_argument(
                "--order",
                type=int,
                default=1,
                choices=[1, 2],
                help="2 for quadratic cells with mid-edge nodes on the exact geometry",
            )
            p.add_argument(
                "--max-cells",
                type=int,
                default=_primitives.DEFAULT_MAX_CELLS,
                help="refuse above this many cells (0 lifts the limit)",
            )
        p.add_argument(
            "--quiet", "-q", action="store_true", help="suppress the summary"
        )


def _doc_line(name):
    fn = _SHAPES[name][0]
    return (fn.__doc__ or "").strip().splitlines()[0]


def _parse(key, text):
    if key in _VECTOR_KINDS:
        return [_VECTOR_KINDS[key](x) for x in text.split(",")]
    if key in _INTEGERS:
        return int(text)
    return float(text)


def generate_cmd(args):
    fn, values, choices = _SHAPES[args.shape]
    kwargs = {}
    for key in values:
        text = getattr(args, key)
        if text is not None:
            kwargs[key] = _parse(key, text)
    for key in choices:
        kwargs[key] = getattr(args, key)
    if args.shape == "reference-cells":
        kwargs["types"] = args.types
        kwargs["pairs"] = args.pairs
    else:
        kwargs["order"] = args.order
        kwargs["max_cells"] = args.max_cells
    mesh = fn(**kwargs)

    if not args.quiet:
        print(f"generated {args.shape}")
        print(f"  points: {len(mesh.points)}")
        for block in mesh.cells:
            print(f"  {block.type}: {len(block.data)}")
        for region in mesh.regions:
            print(f"  region {region.name}: {len(region.entries)}")

    write(args.outfile, mesh, file_format=args.output_format)
    return 0
