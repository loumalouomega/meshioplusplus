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
    group.add_argument(
        "--where",
        type=str,
        help="a cell_data predicate, e.g. 'sdf:distance < 0' or "
        "'quality:scaled_jacobian < 0.3'. A non-finite cell value never matches. "
        "Compose it with `distance_to_surface --location center` for "
        "inside/outside a surface",
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["all", "any"],
        default="all",
        help="keep cell if ALL (default) or ANY node is inside. Applies to "
        "--bbox/--plane only: a --where predicate is already one value per cell "
        "and has nothing to reduce",
    )
    parser.add_argument(
        "--record-ids",
        action="store_true",
        help="attach original point/cell ids as data arrays",
    )


def _floats(text):
    return [float(x) for x in text.split(",")]


#: Longest first, so "<=" is not read as "<" followed by a stray "=".
_COMPARES = ("<=", ">=", "==", "!=", "<", ">")


def _parse_where(text):
    """Split ``NAME OP VALUE`` on the comparison operator.

    Scanned longest-operator-first from the RIGHT-hand end of the name, because
    array names routinely contain characters that look like nothing else here
    (`sdf:distance`, `quality:scaled_jacobian`) but never contain an operator.
    """
    for op in _COMPARES:
        idx = text.find(op)
        if idx > 0:
            name = text[:idx].strip()
            value = text[idx + len(op) :].strip()
            if name and value:
                return name, op, float(value)
    raise ValueError(
        "crop: --where expects 'NAME OP VALUE', with OP one of " + ", ".join(_COMPARES)
    )


def crop_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    if args.bbox is not None:
        vals = _floats(args.bbox)
        if len(vals) != 6:
            raise ValueError("crop: --bbox expects 'xmin,ymin,zmin,xmax,ymax,zmax'")
        out = crop(mesh, bbox=vals, mode=args.mode, record_ids=args.record_ids)
    elif args.plane is not None:
        vals = _floats(args.plane)
        if len(vals) != 6:
            raise ValueError("crop: --plane expects 'px,py,pz,nx,ny,nz'")
        out = crop(
            mesh,
            plane=(vals[:3], vals[3:]),
            mode=args.mode,
            record_ids=args.record_ids,
        )
    else:
        out = crop(mesh, where=_parse_where(args.where), record_ids=args.record_ids)

    write(args.outfile, out, file_format=args.output_format)
    return 0
