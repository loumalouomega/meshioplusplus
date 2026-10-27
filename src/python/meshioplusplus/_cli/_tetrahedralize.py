import json

from .._helpers import _writer_map, read, reader_map, write
from .._tetrahedralize import tetrahedralize


def add_args(parser):
    parser.add_argument("infile", type=str, help="closed triangle surface")
    parser.add_argument("outfile", type=str, help="tetrahedra to be written to")
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
        "--max-size",
        type=float,
        default=0.0,
        help="the largest edge asked of the tetrahedra (default 0: none)",
    )
    parser.add_argument(
        "--max-volume",
        type=float,
        default=0.0,
        help="the largest volume asked of a tetrahedron (default 0: none)",
    )
    parser.add_argument(
        "--grading",
        type=float,
        default=0.3,
        help="how fast the size grows away from the surface, in [0, 1] (default 0.3)",
    )
    parser.add_argument(
        "--optimize-steps",
        type=int,
        default=3,
        help="passes of the backend's volume optimisation, 0 to 100 (default 3)",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=50_000_000,
        help="refuse a result of more tetrahedra (default 50000000)",
    )
    parser.add_argument(
        "--face-array",
        type=str,
        default="cad:face",
        help="the cell data naming each triangle's CAD face (default cad:face)",
    )
    parser.add_argument(
        "--backend",
        choices=("auto", "netgen", "octree"),
        default="auto",
        help="the tetrahedra's backend: netgen keeps the surface's triangles, octree "
        "(native) fits them; auto is netgen when built in, else octree (default auto)",
    )
    parser.add_argument("--json", action="store_true", help="print the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="do not print the summary"
    )


def tetrahedralize_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = tetrahedralize(
        mesh,
        backend=args.backend,
        max_size=args.max_size,
        max_volume=args.max_volume,
        grading=args.grading,
        optimize_steps=args.optimize_steps,
        max_cells=args.max_cells,
        face_array=args.face_array,
        return_report=True,
    )
    if args.json:
        print(json.dumps(report))
    elif not args.quiet:
        print(
            f"{report['num_tetrahedra']} tetrahedra in {report['num_regions']} "
            f"region(s), {report['num_new_points']} new points "
            f"({report['num_points']} in all); dihedral "
            f"{report['min_dihedral_deg']:.2f} to {report['max_dihedral_deg']:.2f} deg, "
            f"volume {report['volume']:.6g}; {report['backend']}, boundary "
            f"{'preserved' if report['boundary_preserved'] else 'NOT preserved'} "
            f"({report['num_boundary_triangles']} triangles)"
        )
    write(args.outfile, out, file_format=args.output_format)
