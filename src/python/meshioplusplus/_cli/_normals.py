from .._helpers import _writer_map, read, reader_map, write
from .._normals import compute_normals


def add_args(parser):
    parser.add_argument("infile", type=str, help="surface mesh to be read from")
    parser.add_argument(
        "outfile", type=str, help="mesh with the normals attached to be written to"
    )
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
        "--cell", action="store_true", help="also attach the cell normals"
    )
    parser.add_argument(
        "--no-point", action="store_true", help="do not attach the point normals"
    )
    parser.add_argument(
        "--weight",
        type=str,
        choices=["angle", "area"],
        default="angle",
        help="how incident faces are weighted into a point normal (default: angle)",
    )
    parser.add_argument(
        "--split-angle",
        type=float,
        default=None,
        metavar="DEG",
        help="duplicate points where the surface creases by more than DEG degrees "
        "(0-180), so every point carries one normal (default: no split)",
    )
    parser.add_argument(
        "--record-parent-ids",
        action="store_true",
        help="also attach normals:parent_point, each output point's input point",
    )
    parser.add_argument(
        "--region",
        type=str,
        default="",
        help="restrict to this named cell region (default: every surface cell)",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def normals_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    out, report = compute_normals(
        mesh,
        point_normals=not args.no_point,
        cell_normals=args.cell,
        weight=args.weight,
        split_angle=args.split_angle,
        region=args.region,
        record_parent_ids=args.record_parent_ids,
        return_report=True,
    )

    if not args.quiet:
        q = report["quality"]
        print(f"normals ({args.weight})")
        print(f"  points added by the split: {report['num_added_points']}")
        print(f"  points split:              {report['num_split_points']}")
        print(f"  isolated points (NaN):     {report['num_isolated']}")
        print(f"  undefined points (NaN):    {report['num_undefined']}")
        print(f"  degenerate triangles:      {report['num_degenerate']}")
        print(
            "  surface: "
            + (
                "watertight"
                if q["watertight"]
                else (
                    f"{q['boundary_edges']} boundary / "
                    f"{q['non_manifold_edges']} non-manifold / "
                    f"{q['inconsistent_pairs']} inconsistent edge(s)"
                )
            )
        )

    write(args.outfile, out, file_format=args.output_format)
    return 0
