from .._agglomerate import agglomerate
from .._helpers import _writer_map, read, reader_map, write
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="coarsened mesh to be written to")
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
        "--target-group-size",
        type=int,
        default=8,
        help="approximate member cells per output group (default: 8)",
    )
    parser.add_argument(
        "--merge-coplanar-faces",
        action="store_true",
        help="fuse the coplanar faces two groups (or a group and the boundary) "
        "share into one polygon",
    )
    parser.add_argument(
        "--coplanar-angle",
        type=float,
        default=1.0,
        metavar="DEG",
        help="largest normal deviation still coplanar (default: 1)",
    )
    parser.add_argument(
        "--min-sphericity",
        type=float,
        default=0.0,
        metavar="S",
        help="refuse an absorption that would drop a group's sphericity below S "
        "(default: 0, off)",
    )
    parser.add_argument("--json", action="store_true", help="emit a summary as JSON")


def agglomerate_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = agglomerate(
        mesh,
        target_group_size=args.target_group_size,
        merge_coplanar_faces=args.merge_coplanar_faces,
        coplanar_angle=args.coplanar_angle,
        min_sphericity=args.min_sphericity,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(
            {
                "cells_in": sum(len(c) for c in mesh.cells),
                "cells_out": sum(len(c) for c in out.cells),
                "num_faces_merged": report["num_faces_merged"],
                "num_rejected": report["num_rejected"],
            }
        )
    return 0
