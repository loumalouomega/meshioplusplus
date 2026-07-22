from .._helpers import _writer_map, read, reader_map, write
from .._transform import transform


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="transformed mesh to be written to")
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
    group = parser.add_argument_group("transform (choose one)")
    group.add_argument("--translate", type=str, help="translation 'x,y,z'")
    group.add_argument(
        "--scale", type=str, help="scale 'sx,sy,sz' or a single scalar 's'"
    )
    group.add_argument(
        "--rotate", type=str, help="rotation 'axis,deg' (axis: x|y|z or nx,ny,nz)"
    )
    group.add_argument(
        "--matrix", type=str, help="row-major 4x4 affine matrix 'm00,...,m33'"
    )
    group.add_argument(
        "--scale-units", type=float, help="uniform unit-scale factor (e.g. 0.001)"
    )
    parser.add_argument(
        "--rotate-data",
        action="store_true",
        help="also rotate vector/tensor point_data by the transform",
    )


def _floats(text):
    return [float(x) for x in text.split(",")]


def _parse_rotate(text):
    parts = text.split(",")
    if len(parts) == 2 and parts[0].strip().lower() in ("x", "y", "z"):
        return (parts[0].strip().lower(), float(parts[1]))
    if len(parts) == 4:  # nx,ny,nz,deg
        return ([float(parts[0]), float(parts[1]), float(parts[2])], float(parts[3]))
    raise ValueError(
        "transform: --rotate expects 'axis,deg' (axis x|y|z) or 'nx,ny,nz,deg'"
    )


def transform_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    kwargs = {}
    if args.translate is not None:
        kwargs["translate"] = _floats(args.translate)
    if args.scale is not None:
        vals = _floats(args.scale)
        kwargs["scale"] = vals[0] if len(vals) == 1 else vals
    if args.rotate is not None:
        kwargs["rotate"] = _parse_rotate(args.rotate)
    if args.matrix is not None:
        kwargs["matrix"] = _floats(args.matrix)
    if args.scale_units is not None:
        kwargs["scale_units"] = args.scale_units
    if not kwargs:
        raise ValueError(
            "transform: give one of --translate/--scale/--rotate/--matrix/--scale-units"
        )

    out = transform(mesh, rotate_vector_data=args.rotate_data, **kwargs)
    write(args.outfile, out, file_format=args.output_format)
    return 0
