from .._helpers import _writer_map, read, reader_map, write
from .._merge import merge


def add_args(parser):
    parser.add_argument(
        "files",
        nargs="+",
        type=str,
        help="two or more input meshes followed by the output file",
    )
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format (applied to every input)",
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
        "--weld",
        action="store_true",
        help="merge coincident nodes within --atol",
    )
    parser.add_argument(
        "--atol",
        type=float,
        default=1e-8,
        help="coincidence tolerance for --weld (default: 1e-8)",
    )
    parser.add_argument(
        "--data-policy",
        type=str,
        choices=["intersection", "fill"],
        default="intersection",
        help="how to treat data keys missing from some inputs (default: intersection)",
    )
    parser.add_argument(
        "--drop-duplicate-cells",
        action="store_true",
        help="with --weld, drop cells that become identical after welding",
    )
    parser.add_argument(
        "--no-source-tag",
        action="store_true",
        help="do not add the per-cell source_mesh_id tag",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="do not print the merge summary",
    )


def merge_cmd(args):
    if len(args.files) < 2:
        raise SystemExit(
            "merge: need at least one input mesh and an output file "
            "(e.g. `meshioplusplus merge a.vtu b.vtu out.vtu`)"
        )
    *infiles, outfile = args.files
    meshes = [read(f, file_format=args.input_format) for f in infiles]

    total_in_points = sum(len(m.points) for m in meshes)
    total_in_cells = sum(len(cb.data) for m in meshes for cb in m.cells)

    out = merge(
        meshes,
        weld=args.weld,
        atol=args.atol,
        source_tag=not args.no_source_tag,
        data_policy=args.data_policy,
        drop_duplicate_cells=args.drop_duplicate_cells,
    )

    out_points = len(out.points)
    out_cells = sum(len(cb.data) for cb in out.cells)
    if not args.quiet:
        print(f"merged {len(meshes)} meshes")
        print(f"  points in:  {total_in_points}")
        print(f"  cells in:   {total_in_cells}")
        if args.weld:
            print(f"  points welded: {total_in_points - out_points}")
        print(f"  points out: {out_points}")
        print(f"  cells out:  {out_cells}")

    write(outfile, out, file_format=args.output_format)
    return 0
