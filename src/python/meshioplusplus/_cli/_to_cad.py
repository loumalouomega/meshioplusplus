from .._fit_cad import fit_cad
from .._fit_primitives import TYPES, fit_primitives
from .._helpers import read, reader_map, write
from .._segment import segment
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="closed surface mesh to be read from")
    parser.add_argument("outfile", type=str, help="STEP file the solid is written to")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--curvature-tolerance",
        type=float,
        default=0.0,
        help="segment's curvature test, for a skin without segment:patch "
        "(default: 0, off)",
    )
    parser.add_argument(
        "--min-patch-cells",
        type=int,
        default=1,
        help="segment's smallest patch, merged below it, for a skin without "
        "segment:patch (default: 1; a noisy scan wants tens)",
    )
    parser.add_argument(
        "--fit-tolerance",
        type=float,
        default=1e-3,
        help="fit_primitives' relative tolerance, for a skin without the "
        "fit:patch_* table (default: 1e-3)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=None,
        help="largest distance a face may lie from the skin "
        "(default: 3x the largest accepted patch tolerance)",
    )
    parser.add_argument(
        "--no-exact-edges",
        action="store_true",
        help="approximate every edge between analytic faces with a B-spline",
    )
    parser.add_argument(
        "--types",
        type=str,
        default=",".join(TYPES[1:]),
        help="comma-separated fitted types kept as analytic faces (default: all)",
    )
    parser.add_argument(
        "--mesh-deflection",
        type=float,
        default=None,
        help="chordal deflection of the re-meshed solid (default: tolerance / 4)",
    )
    parser.add_argument(
        "--max-demotion-rounds",
        type=int,
        default=8,
        help="rounds of demoting failed faces before faceting all (default: 8)",
    )
    parser.add_argument(
        "--mesh-out", type=str, default=None, help="also write the re-meshed solid here"
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def to_cad_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    if "segment:patch" not in mesh.cell_data:
        mesh = segment(
            mesh,
            curvature_tolerance=args.curvature_tolerance,
            min_patch_cells=args.min_patch_cells,
        )
    if "fit:patch_id" not in mesh.field_data:
        mesh = fit_primitives(mesh, tolerance=args.fit_tolerance)
    out, report = fit_cad(
        mesh,
        args.outfile,
        tolerance=args.tolerance,
        exact_edges=not args.no_exact_edges,
        types=[t.strip() for t in args.types.split(",") if t.strip()],
        mesh_deflection=args.mesh_deflection,
        max_demotion_rounds=args.max_demotion_rounds,
        return_report=True,
    )
    if args.mesh_out:
        write(args.mesh_out, out)
    if args.json:
        emit_json(report)
    elif not args.quiet:
        print(
            f"to-cad: {report['num_faces']} face(s) in {report['num_solids']} "
            f"solid(s), {'valid' if report['valid'] else 'NOT valid'}"
        )
        for t in TYPES[1:]:
            print(f"  {t + ':':10s} {report['num_' + t]}")
        print(f"  {'facets:':10s} {report['num_faceted_faces']}")
        print(f"  {'demoted:':10s} {report['num_demoted']} patch(es)")
        print(
            f"  volume ratio {report['volume_ratio']:.6g}, max deviation "
            f"{report['max_deviation']:.6g} (tolerance {report['tolerance']:.6g})"
        )
    return 0 if report["valid"] else 2
