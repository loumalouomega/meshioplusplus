import numpy as np

from .._fit_cad import DEFAULT_RANSAC_SEED, fit_cad
from .._fit_primitives import TYPES, fit_primitives
from .._helpers import read, reader_map, write
from .._segment import segment
from ._json import emit_json
from ._refine_patches import parse_seed


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
        "--refine",
        action="store_true",
        help="first refine the patches by variational shape approximation "
        "(refine_patches), so more of the skin becomes analytic faces",
    )
    parser.add_argument(
        "--max-refine-iterations",
        type=int,
        default=10,
        metavar="N",
        help="--refine: the most iterations (default: 10)",
    )
    parser.add_argument(
        "--refine-epsilon",
        type=float,
        default=1e-3,
        metavar="E",
        help="--refine: stop when an iteration adds less than this to the area "
        "share primitives explain (default: 1e-3)",
    )
    parser.add_argument(
        "--normal-weight",
        type=float,
        default=1.0,
        metavar="W",
        help="--refine: the weight of the normal term in a cell's cost (default: 1)",
    )
    parser.add_argument(
        "--no-refine-merge",
        action="store_true",
        help="--refine: do not merge neighbouring patches one primitive explains",
    )
    parser.add_argument(
        "--no-refine-split",
        action="store_true",
        help="--refine: do not split a primitive out of a freeform patch",
    )
    parser.add_argument(
        "--ransac-seed",
        type=parse_seed,
        default=DEFAULT_RANSAC_SEED,
        metavar="S",
        help="--refine: the seed of the RANSAC split (default: 0x6d657368696f2b2b)",
    )
    parser.add_argument(
        "--refine-tolerance",
        type=float,
        default=1e-3,
        metavar="T",
        help="--refine: how far a point may lie from its patch's primitive, "
        "relative to the diagonal (default: 1e-3)",
    )
    parser.add_argument(
        "--freeform-general",
        action="store_true",
        help="try a patch the plane and cylinder freeform fit leaves unexplained "
        "as a general freeform B-spline face (disc parametrisation, iterative fit)",
    )
    parser.add_argument(
        "--freeform-tolerance",
        type=float,
        default=0.0,
        metavar="T",
        help="--freeform-general: the fit's target at the vertices (default: 0, "
        "a quarter of the tolerance)",
    )
    parser.add_argument(
        "--freeform-smoothing",
        type=float,
        default=1e-12,
        metavar="W",
        help="--freeform-general: the thin-plate weight (default: 1e-12)",
    )
    parser.add_argument(
        "--max-freeform-poles",
        type=int,
        default=40,
        metavar="N",
        help="--freeform-general: the most poles per direction, 4 to 1024 "
        "(default: 40)",
    )
    parser.add_argument(
        "--max-freeform-iterations",
        type=int,
        default=30,
        metavar="N",
        help="--freeform-general: the iteration budget (default: 30)",
    )
    parser.add_argument(
        "--max-freeform-splits",
        type=int,
        default=2,
        metavar="N",
        help="--freeform-general: how many times a patch the fit cannot explain "
        "is halved (default: 2)",
    )
    parser.add_argument(
        "--sew-tolerance",
        type=float,
        default=0.0,
        metavar="T",
        help="the tolerance the faces are sewn at (default: 0, the tolerance)",
    )
    parser.add_argument(
        "--extend-surfaces",
        action="store_true",
        help="extend freeform B-spline surfaces past their patches before the "
        "edges are built",
    )
    parser.add_argument(
        "--snap-junctions",
        action="store_true",
        help="place the vertices where analytic faces meet by a least-squares "
        "solve on their tangent planes",
    )
    parser.add_argument(
        "--face-deviation",
        action="store_true",
        help="measure each face's two-sided distance to the skin "
        "(fit_cad:face_hausdorff, fit_cad:face_rms, hausdorff_max/_rms)",
    )
    parser.add_argument(
        "--freeform-mode",
        choices=["legacy", "adaptive", "network", "auto"],
        default="legacy",
        help="how patches no primitive or feature fit explains are treated "
        "(default legacy: every output as before)",
    )
    parser.add_argument(
        "--freeform-parameterization",
        choices=["auto", "projection", "lscm", "harmonic"],
        default="auto",
        help="the freeform parametrisation outside the legacy mode",
    )
    parser.add_argument(
        "--g1-tolerance",
        type=float,
        default=1.0,
        help="largest normal jump across an edge between network faces, degrees",
    )
    parser.add_argument(
        "--extraordinary-vertex",
        choices=["g0", "g1", "g1x", "g1xx"],
        default="g1x",
        help="a patch network's continuity at an extraordinary vertex",
    )
    parser.add_argument(
        "--local-refine",
        choices=["knots", "thb"],
        default="knots",
        help="how a freeform patch is refined (thb is experimental)",
    )
    parser.add_argument(
        "--freeform-max-control-points",
        type=int,
        default=400,
        help="the most poles of one freeform face outside the legacy mode",
    )
    parser.add_argument(
        "--layout-target-patch-size",
        type=float,
        default=0.1,
        help="a quad-layout patch's target size relative to the diagonal",
    )
    parser.add_argument(
        "--layout-max-valence-deviation",
        type=int,
        default=2,
        help="the largest departure from valence 4 of a layout vertex (1 to 3)",
    )
    parser.add_argument(
        "--resegment-rounds",
        type=int,
        default=2,
        help="re-segmentation rounds of a residual region before demotion",
    )
    parser.add_argument(
        "--no-layout-snap-to-features",
        action="store_true",
        help="do not snap layout edges to feature edges",
    )
    parser.add_argument(
        "--no-validate-by-surface-sampling",
        action="store_true",
        help="accept a freeform face without sampling its surface against the skin",
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
        refine_segmentation=args.refine,
        max_refine_iterations=args.max_refine_iterations,
        refine_epsilon=args.refine_epsilon,
        normal_weight=args.normal_weight,
        allow_merge=not args.no_refine_merge,
        allow_split=not args.no_refine_split,
        ransac_seed=args.ransac_seed,
        refine_tolerance=args.refine_tolerance,
        freeform_general=args.freeform_general,
        freeform_tolerance=args.freeform_tolerance,
        freeform_smoothing=args.freeform_smoothing,
        max_freeform_poles=args.max_freeform_poles,
        max_freeform_iterations=args.max_freeform_iterations,
        max_freeform_splits=args.max_freeform_splits,
        sew_tolerance=args.sew_tolerance,
        extend_surfaces=args.extend_surfaces,
        snap_junctions=args.snap_junctions,
        face_deviation=args.face_deviation,
        freeform_mode=args.freeform_mode,
        freeform_parameterization=args.freeform_parameterization,
        extraordinary_vertex=args.extraordinary_vertex,
        local_refine=args.local_refine,
        freeform_g1_tolerance_deg=args.g1_tolerance,
        freeform_max_control_points=args.freeform_max_control_points,
        layout_target_patch_size=args.layout_target_patch_size,
        layout_max_valence_deviation=args.layout_max_valence_deviation,
        resegment_rounds=args.resegment_rounds,
        layout_snap_to_features=not args.no_layout_snap_to_features,
        validate_by_surface_sampling=not args.no_validate_by_surface_sampling,
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
            f"{'valid' if report['valid'] else 'NOT valid'}, "
            f"{'closed' if report['closed'] else 'NOT closed'}"
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
        if report.get("refined"):
            print(
                f"  refined: {report['num_refine_splits']} split(s), "
                f"{report['num_refine_merges']} merge(s) in "
                f"{report['refine_iterations']} iteration(s); explained area "
                f"{report['initial_explained_area_fraction']:.6g} -> "
                f"{report['explained_area_fraction']:.6g}; "
                + (
                    "refined patches kept"
                    if report["refinement_kept"]
                    else "original patches kept"
                )
            )
        if "num_freeform_general" in report:
            print(
                f"  general freeform: {report['num_freeform_general']} face(s), "
                f"{report['num_freeform_refinements']} knot insertion round(s), "
                f"{report['avg_freeform_poles']:.6g} pole(s) on average"
            )
        if "num_extended_edges" in report:
            print(f"  extended: {report['num_extended_edges']} surface(s)")
        if "num_snapped_junctions" in report:
            print(f"  snapped:  {report['num_snapped_junctions']} junction(s)")
        if "hausdorff_max" in report:
            print(
                f"  face deviation: Hausdorff {report['hausdorff_max']:.6g}, "
                f"RMS {report['hausdorff_rms']:.6g}"
            )
        if "freeform_mode" in report:
            print(
                f"  freeform ({report['freeform_mode']}): "
                f"{report['num_freeform_adaptive']} adaptive face(s), "
                f"{report['num_network_patches']} network patch(es), "
                f"{report['num_escalations']} escalation(s), facet area "
                f"{report['facet_area_fraction']:.4g}"
            )
        print(
            f"  volume ratio {report['volume_ratio']:.6g}, max deviation "
            f"{report['max_deviation']:.6g} (tolerance {report['tolerance']:.6g})"
        )
    return 0 if report["valid"] else 2
