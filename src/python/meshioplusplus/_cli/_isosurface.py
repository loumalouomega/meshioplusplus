from .._helpers import _writer_map, read, reader_map, write
from .._isosurface import isosurface


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="contour mesh to be written to")
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
        "--array",
        type=str,
        required=True,
        help="name of the point_data array to contour "
        "(cell_data has no level set -- convert it with 'data to-point' first)",
    )
    parser.add_argument(
        "--values",
        type=str,
        required=True,
        help="isovalues 'v1,v2,...' (negatives need --values=-1,0)",
    )
    parser.add_argument(
        "--component",
        type=int,
        default=None,
        help="component of a multi-component array; the row magnitude by default",
    )
    parser.add_argument(
        "--record-parent-ids",
        action="store_true",
        help="attach the originating input cell id as iso:parent_cell cell_data",
    )


def isosurface_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    values = [float(x) for x in args.values.split(",")]
    if not values:
        raise ValueError("isosurface: --values expects at least one isovalue")

    out = isosurface(
        mesh,
        args.array,
        values,
        component=args.component,
        record_parent_ids=args.record_parent_ids,
    )

    write(args.outfile, out, file_format=args.output_format)
    return 0
