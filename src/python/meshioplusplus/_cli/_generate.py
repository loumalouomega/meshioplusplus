from .. import _primitives
from .._helpers import _writer_map, write

QUAD_OR_TRIANGLE = {"cell_type": ("quad", ["quad", "triangle"])}
HEX_OR_TETRA = {"cell_type": ("hexahedron", ["hexahedron", "tetra"])}

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
    "reference-triangle": (
        _primitives.reference_triangle,
        {"n": None, "size": "lx,ly", "origin": "x,y[,z]"},
        {},
    ),
    "ngon": (
        _primitives.ngon,
        {"sides": None, "radius": None, "rings": None, "center": "x,y,z"},
        {},
    ),
    "annulus": (
        _primitives.annulus,
        {
            "inner_radius": None,
            "outer_radius": None,
            "ntheta": None,
            "nr": None,
            "angle": None,
            "center": "x,y,z",
        },
        {"cell_type": ("quad", ["quad", "triangle"])},
    ),
    "lshape": (
        _primitives.lshape,
        {"n": None, "size": None, "center": "x,y,z"},
        {"cell_type": ("quad", ["quad", "triangle"])},
    ),
    "plate-with-hole": (
        _primitives.plate_with_hole,
        {
            "radius": None,
            "size": None,
            "nt": None,
            "nr": None,
            "center": "x,y,z",
        },
        {"cell_type": ("quad", ["quad", "triangle"])},
    ),
    "honeycomb": (
        _primitives.honeycomb,
        {"dims": "nx,ny", "edge": None, "origin": "x,y[,z]"},
        {},
    ),
    "reference-cells": (_primitives.reference_cells, {}, {}),
    "line": (
        _primitives.line,
        {"n": None, "start": "x,y,z", "end": "x,y,z"},
        {},
    ),
    "circle": (
        _primitives.circle,
        {"radius": None, "n": None, "center": "x,y,z"},
        {},
    ),
    "helix": (
        _primitives.helix,
        {
            "radius": None,
            "pitch": None,
            "turns": None,
            "n": None,
            "center": "x,y,z",
        },
        {},
    ),
    "sphere-points": (
        _primitives.sphere_points,
        {"n": None, "radius": None, "center": "x,y,z"},
        {},
    ),
    "box-points": (
        _primitives.box_points,
        {"n": None, "size": "lx,ly,lz", "origin": "x,y,z", "seed": None},
        {},
    ),
    "tetra-sphere": (
        _primitives.tetra_sphere,
        {"radius": None, "subdivisions": None, "center": "x,y,z"},
        {},
    ),
    "octa-sphere": (
        _primitives.octa_sphere,
        {"radius": None, "subdivisions": None, "center": "x,y,z"},
        {},
    ),
    "cube-sphere": (
        _primitives.cube_sphere,
        {"radius": None, "n": None, "center": "x,y,z"},
        QUAD_OR_TRIANGLE,
    ),
    "geo-sphere": (
        _primitives.geo_sphere,
        {"radius": None, "n": None, "center": "x,y,z"},
        {},
    ),
    "torus": (
        _primitives.torus,
        {
            "major_radius": None,
            "minor_radius": None,
            "nu": None,
            "nv": None,
            "center": "x,y,z",
        },
        QUAD_OR_TRIANGLE,
    ),
    "moebius": (
        _primitives.moebius,
        {
            "radius": None,
            "width": None,
            "nu": None,
            "nv": None,
            "num_twists": None,
            "center": "x,y,z",
        },
        QUAD_OR_TRIANGLE,
    ),
    "tube": (
        _primitives.tube,
        {"radius": None, "height": None, "ntheta": None, "nz": None, "center": "x,y,z"},
        QUAD_OR_TRIANGLE,
    ),
    "cone": (
        _primitives.cone,
        {"radius": None, "height": None, "ntheta": None, "nz": None, "center": "x,y,z"},
        QUAD_OR_TRIANGLE,
    ),
    "capsule": (
        _primitives.capsule,
        {
            "radius": None,
            "height": None,
            "ntheta": None,
            "nz": None,
            "ncap": None,
            "center": "x,y,z",
        },
        QUAD_OR_TRIANGLE,
    ),
    "box-surface": (
        _primitives.box_surface,
        {"dims": "nx,ny,nz", "size": "lx,ly,lz", "origin": "x,y,z"},
        QUAD_OR_TRIANGLE,
    ),
    "platonic": (
        _primitives.platonic,
        {"radius": None, "center": "x,y,z"},
        {
            "name": (
                "icosahedron",
                ["tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron"],
            )
        },
    ),
    "ball": (
        _primitives.ball,
        {"radius": None, "n": None, "center": "x,y,z"},
        HEX_OR_TETRA,
    ),
    "spherical-shell": (
        _primitives.spherical_shell,
        {
            "inner_radius": None,
            "outer_radius": None,
            "n": None,
            "nr": None,
            "ratio": None,
            "center": "x,y,z",
        },
        {},
    ),
    "pipe": (
        _primitives.pipe,
        {
            "inner_radius": None,
            "outer_radius": None,
            "height": None,
            "ntheta": None,
            "nr": None,
            "nz": None,
            "ratio": None,
            "center": "x,y,z",
        },
        {},
    ),
    "fichera": (
        _primitives.fichera,
        {"n": None, "size": None, "center": "x,y,z"},
        HEX_OR_TETRA,
    ),
    "solid-torus": (
        _primitives.solid_torus,
        {
            "major_radius": None,
            "minor_radius": None,
            "n": None,
            "nu": None,
            "center": "x,y,z",
        },
        {},
    ),
}

_INTEGERS = {
    "n",
    "nu",
    "nv",
    "nxy",
    "nz",
    "subdivisions",
    "sides",
    "rings",
    "ntheta",
    "nr",
    "nt",
    "ncap",
    "num_twists",
    "seed",
}
_FLOATS = {
    "radius",
    "height",
    "inner_radius",
    "outer_radius",
    "angle",
    "edge",
    "pitch",
    "turns",
    "major_radius",
    "minor_radius",
    "width",
    "ratio",
}
_VECTOR_KINDS = {
    "dims": int,
    "size": float,
    "origin": float,
    "center": float,
    "start": float,
    "end": float,
}
#: shapes whose ``size`` is one number (a side), not a vector
_SCALAR_SIZE = {"lshape", "plate-with-hole", "fichera"}
#: shapes with no ``order`` option (a hexagon, a point or a pentagon face has no
#: quadratic form)
_NO_ORDER = {"reference-cells", "honeycomb", "sphere-points", "box-points", "platonic"}


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
                help=f"{key.replace('_', ' ')} (default {default})",
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
            if name not in _NO_ORDER:
                p.add_argument(
                    "--order",
                    type=int,
                    default=1,
                    choices=list(range(1, _primitives.MAX_ORDER + 1)),
                    help="2 for quadratic cells with mid-edge nodes on the exact "
                    "geometry, 3 to 10 for the VTK_LAGRANGE_* family",
                )
                p.add_argument(
                    "--complete",
                    action="store_true",
                    help="with --order 2, quad9 / hexahedron27 / wedge18 instead of "
                    "the serendipity types",
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


def _parse(key, text, shape=None):
    if key == "size" and shape in _SCALAR_SIZE:
        return float(text)
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
            kwargs[key] = _parse(key, text, args.shape)
    for key in choices:
        kwargs[key] = getattr(args, key)
    if args.shape == "reference-cells":
        kwargs["types"] = args.types
        kwargs["pairs"] = args.pairs
    else:
        if args.shape not in _NO_ORDER:
            kwargs["order"] = args.order
            kwargs["complete"] = args.complete
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
