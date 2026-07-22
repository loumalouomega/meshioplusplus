from .._clean import clean
from .._helpers import _writer_map, read, reader_map, write


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="cleaned mesh to be written to")
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
        "--weld", action="store_true", help="fuse coincident points within --atol"
    )
    parser.add_argument(
        "--atol", type=float, default=1e-8, help="weld tolerance (default 1e-8)"
    )
    parser.add_argument(
        "--remove-orphans", action="store_true", help="drop unused points"
    )
    parser.add_argument(
        "--drop-degenerate", action="store_true", help="drop degenerate cells"
    )
    parser.add_argument(
        "--drop-duplicates", action="store_true", help="drop exact-duplicate cells"
    )


def clean_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    # With no step flags given, run the sensible default set (prune + de-dup, no
    # weld). If any step flag is given, honor exactly the ones requested.
    any_flag = (
        args.weld or args.remove_orphans or args.drop_degenerate or args.drop_duplicates
    )
    if any_flag:
        opts = dict(
            weld=args.weld,
            atol=args.atol,
            remove_orphans=args.remove_orphans,
            drop_degenerate=args.drop_degenerate,
            drop_duplicate_cells=args.drop_duplicates,
        )
    else:
        opts = dict(
            weld=False,
            atol=args.atol,
            remove_orphans=True,
            drop_degenerate=True,
            drop_duplicate_cells=True,
        )

    out, report = clean(mesh, return_report=True, **opts)
    print("cleaned mesh")
    print(f"  points welded:            {report['points_welded']}")
    print(f"  points removed (orphan):  {report['points_removed_orphan']}")
    print(f"  cells dropped (degenerate): {report['cells_dropped_degenerate']}")
    print(f"  cells dropped (duplicate):  {report['cells_dropped_duplicate']}")
    write(args.outfile, out, file_format=args.output_format)
    return 0
