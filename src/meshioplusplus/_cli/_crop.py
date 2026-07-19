from .._crop import crop
from .._helpers import _writer_map, read, reader_map, write


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="cropped mesh to be written to")
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
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--bbox", type=str, help="bounding box 'xmin,ymin,zmin,xmax,ymax,zmax'"
    )
    group.add_argument(
        "--plane", type=str, help="half-space 'px,py,pz,nx,ny,nz' (point + normal)"
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["all", "any"],
        default="all",
        help="keep cell if ALL (default) or ANY node is inside",
    )
    parser.add_argument(
        "--record-ids",
        action="store_true",
        help="attach original point/cell ids as data arrays",
    )


def _floats(text):
    return [float(x) for x in text.split(",")]


def crop_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    if args.bbox is not None:
        vals = _floats(args.bbox)
        if len(vals) != 6:
            raise ValueError("crop: --bbox expects 'xmin,ymin,zmin,xmax,ymax,zmax'")
        out = crop(mesh, bbox=vals, mode=args.mode, record_ids=args.record_ids)
    else:
        vals = _floats(args.plane)
        if len(vals) != 6:
            raise ValueError("crop: --plane expects 'px,py,pz,nx,ny,nz'")
        out = crop(
            mesh,
            plane=(vals[:3], vals[3:]),
            mode=args.mode,
            record_ids=args.record_ids,
        )

    write(args.outfile, out, file_format=args.output_format)
    return 0
