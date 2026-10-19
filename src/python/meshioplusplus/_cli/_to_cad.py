import numpy as np

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
        "--segment-method",
        type=str,
        choices=["primitive", "angle"],
        default="primitive",
        help="how a skin without segment:patch is segmented: where the surface "
        "its facets lie on changes, or by angles alone (default: primitive)",
    )
    parser.add_argument(
        "--segment-tolerance",
        type=float,
        default=1e-3,
        help="the primitive segmentation's tolerance, relative to the diagonal "
        "(default: 1e-3)",
    )
    parser.add_argument(
        "--curvature-tolerance",
        type=float,
        default=0.0,
        help="segment's curvature test, for a skin without segment:patch, angle "
        "method (default: 0, off)",
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
        default=None,
        help="fit_primitives' relative tolerance, for a skin without the "
        "fit:patch_* table (default: the primitive segmentation's own, as a "
        "length; 1e-3 after the angle method)",
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
        "--no-freeform",
        action="store_true",
        help="facet the patches no primitive explains instead of fitting B-spline faces",
    )
    parser.add_argument(
        "--freeform-degree",
        type=int,
        default=3,
        help="highest degree of a freeform B-spline face, 3 to 8 (default: 3)",
    )
    parser.add_argument(
        "--freeform-max-control",
        type=int,
        default=64,
        help="most grid points per direction a freeform fit samples (default: 64)",
    )
    parser.add_argument(
        "--no-freeform-split",
        action="store_true",
        help="facet a freeform patch that folds over every flattening instead of "
        "splitting it into pieces",
    )
    parser.add_argument(
        "--no-fillets",
        action="store_true",
        help="do not try a patch between two fitted faces as a rolling-ball fillet "
        "or blend",
    )
    parser.add_argument(
        "--no-threads",
        action="store_true",
        help="do not try a patch as a helical thread face",
    )
    parser.add_argument(
        "--no-merge-faces",
        action="store_true",
        help="keep faces on one surface apart once the solid is built",
    )
    parser.add_argument(
        "--backend",
        type=str,
        choices=["auto", "occt", "native"],
        default="auto",
        help="the engine: OpenCASCADE, or the kernel-free native backend that "
        "writes the fitted patches as analytic STEP faces (default: auto, "
        "OpenCASCADE when the core carries it)",
    )
    parser.add_argument(
        "--mesh-out", type=str, default=None, help="also write the re-meshed solid here"
    )
    parser.add_argument(
        "--cad-json",
        type=str,
        default=None,
        help="also write the solid as Kratos CAD JSON (NURBS faces, loops, edges) here",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def to_cad_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    fit_kwargs = {}
    if args.fit_tolerance is not None:
        fit_kwargs["tolerance"] = args.fit_tolerance
    if "segment:patch" not in mesh.cell_data:
        primitive = args.segment_method == "primitive"
        mesh = segment(
            mesh,
            curvature_tolerance=args.curvature_tolerance,
            min_patch_cells=args.min_patch_cells,
            method=args.segment_method,
            tolerance=args.segment_tolerance,
        )
        if primitive and args.fit_tolerance is None:
            # Fit at the tolerance the patches were grown at, as a length.
            pts = np.asarray(mesh.points, dtype=np.float64)
            diag = float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))
            if diag > 0.0:
                fit_kwargs = {
                    "tolerance": args.segment_tolerance * diag,
                    "relative": False,
                }
    if "fit:patch_id" not in mesh.field_data:
        mesh = fit_primitives(mesh, **fit_kwargs)
    out, report = fit_cad(
        mesh,
        args.outfile,
        tolerance=args.tolerance,
        exact_edges=not args.no_exact_edges,
        types=[t.strip() for t in args.types.split(",") if t.strip()],
        mesh_deflection=args.mesh_deflection,
        max_demotion_rounds=args.max_demotion_rounds,
        freeform=not args.no_freeform,
        freeform_degree=args.freeform_degree,
        freeform_max_control=args.freeform_max_control,
        freeform_split=not args.no_freeform_split,
        fillets=not args.no_fillets,
        threads=not args.no_threads,
        return_report=True,
        merge_faces=not args.no_merge_faces,
        backend=args.backend,
    )
    if args.mesh_out:
        write(args.mesh_out, out)
    if args.cad_json:
        write(args.cad_json, out, file_format="cad_json")
    if args.json:
        emit_json(report)
    elif not args.quiet:
        print(
            f"to-cad ({report['backend']}): {report['num_faces']} face(s) in "
            f"{report['num_solids']} solid(s), "
            f"{'valid' if report['valid'] else 'NOT valid'}"
        )
        for t in TYPES[1:]:
            print(f"  {t + ':':10s} {report['num_' + t]}")
        split = report["num_split_patches"]
        print(
            f"  {'bspline:':10s} {report['num_bspline']}"
            + (f" ({split} split piece(s))" if split > 0 else "")
        )
        variable, blend = report["num_variable_fillet"], report["num_blend"]
        conic = report["num_conic_blend"]
        print(
            f"  {'fillet:':10s} {report['num_fillet']}"
            + (
                f" ({variable} of varying radius, {blend} blend(s))"
                if variable > 0 or blend > 0
                else ""
            )
            + (f" ({conic} conic)" if conic > 0 else "")
        )
        rounded = report["num_rounded_thread"]
        tapered = report["num_tapered_thread"]
        varying = report["num_variable_pitch_thread"]
        print(
            f"  {'helical:':10s} {report['num_helical']}"
            + (f" ({rounded} rounded)" if rounded > 0 else "")
            + (
                f" ({tapered} tapered, {varying} of varying pitch)"
                if tapered > 0 or varying > 0
                else ""
            )
        )
        print(f"  {'chamfer:':10s} {report['num_chamfer']}")
        print(f"  {'vertex:':10s} {report['num_vertex_blend']} vertex blend(s)")
        print(f"  {'facets:':10s} {report['num_faceted_faces']}")
        print(f"  {'demoted:':10s} {report['num_demoted']} patch(es)")
        print(
            f"  volume ratio {report['volume_ratio']:.6g}, max deviation "
            f"{report['max_deviation']:.6g} (tolerance {report['tolerance']:.6g})"
        )
    return 0 if report["valid"] else 2
