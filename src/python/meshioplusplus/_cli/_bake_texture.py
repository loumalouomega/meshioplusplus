"""The ``bake-texture`` verb: a point or cell array baked into a texture through
``uv`` (``doc/bake_texture.md``)."""

import numpy as np

from .._bake_texture import bake_texture
from .._helpers import read, write
from ._common import add_io_args
from ._json import emit_json


def add_args(parser):
    add_io_args(
        parser,
        "output mesh carrying the texture (glTF embeds it, OBJ writes .mtl and .png)",
        in_help="surface mesh with point_data uv",
        output_format=False,
    )
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument("--array", required=True, help="the point or cell data array")
    parser.add_argument("--width", type=int, default=512, help="pixels (default 512)")
    parser.add_argument("--height", type=int, default=512, help="pixels (default 512)")
    parser.add_argument("--cmap", default="viridis", help="colormap (default viridis)")
    parser.add_argument("--vmin", type=float, default=None, help="low end of the range")
    parser.add_argument(
        "--vmax", type=float, default=None, help="high end of the range"
    )
    parser.add_argument(
        "--component", type=int, default=None, help="component of a vector array"
    )
    parser.add_argument("--image", default=None, help="also write the texture as a PNG")
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def bake_texture_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, _image, report = bake_texture(
        mesh,
        args.array,
        size=(args.width, args.height),
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
        component=args.component,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.image:
        with open(args.image, "wb") as f:
            f.write(np.asarray(out.field_data["texture:png"], dtype=np.uint8).tobytes())
    if args.json:
        emit_json(
            {
                "output": args.outfile,
                "vmin": report["vmin"],
                "vmax": report["vmax"],
                "covered_pixels": report["covered_pixels"],
            }
        )
        return 0
    print(
        f"{report['covered_pixels']} pixel(s) covered, range "
        f"[{report['vmin']:.6g}, {report['vmax']:.6g}]"
    )
    return 0
