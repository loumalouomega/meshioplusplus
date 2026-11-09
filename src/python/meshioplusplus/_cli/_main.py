import argparse
import importlib
import sys
from sys import version_info
from typing import NamedTuple

from ..__about__ import __version__


class _Verb(NamedTuple):
    """One subcommand: its name, the ``_cli`` module that implements it, the
    module's handler and argument builder, and its ``--help`` line."""

    name: str
    module: str
    handler: str
    help: str
    add_args: str = "add_args"
    aliases: tuple = ()


# Every verb but ``completion``, in ``--help`` order. The modules are imported
# only for the verb a command line runs (roadmap 5.12.1): building all of them
# loads about 270 modules, a ``--help`` or a single verb needs none of them or
# one. ``test_cli.py`` checks each lazily built subparser against the full one.
_VERBS = (
    _Verb("convert", "_convert", "convert", "Convert mesh files", aliases=("c",)),
    _Verb("info", "_info", "info", "Print mesh info", aliases=("i",)),
    _Verb("compress", "_compress", "compress", "Compress mesh file"),
    _Verb("decompress", "_decompress", "decompress", "Decompress mesh file"),
    _Verb("ascii", "_ascii", "ascii", "Convert to ASCII", aliases=("a",)),
    _Verb("binary", "_binary", "binary", "Convert to binary", aliases=("b",)),
    _Verb(
        "quality", "_quality", "quality", "Print mesh quality metrics", aliases=("q",)
    ),
    _Verb(
        "curvature",
        "_curvature",
        "curvature_cmd",
        "Per-vertex mean and Gaussian curvature of a surface",
    ),
    _Verb(
        "normals",
        "_normals",
        "normals_cmd",
        "Point and cell normals of a surface, optionally split at creases",
    ),
    _Verb(
        "check",
        "_check",
        "check_cmd",
        "Quality gate: exit 1 when cells break thresholds (2: could not check)",
    ),
    _Verb(
        "check-tetrahedra",
        "_check_tetrahedra",
        "check_tetrahedra_cmd",
        (
            "Inverted, degenerate and overlapping tetrahedra: exit 1 when any "
            "(2: could not check)"
        ),
    ),
    _Verb(
        "check-conformity",
        "_check_conformity",
        "check_conformity_cmd",
        (
            "Hanging nodes, T-junctions and partial overlaps: exit 1 when the "
            "mesh is not conforming (2: could not check)"
        ),
    ),
    _Verb(
        "find-intersections",
        "_find_intersections",
        "find_intersections_cmd",
        (
            "Self-intersections of a surface, or where two surfaces meet: exit 1 "
            "when any are found (2: could not check)"
        ),
    ),
    _Verb(
        "thickness",
        "_thickness",
        "thickness_cmd",
        "Wall thickness or clearance per facet, by a ray along its normal",
    ),
    _Verb(
        "deviation",
        "_deviation",
        "deviation_cmd",
        "Signed deviation of a scan from a reference surface or fitted primitive, "
        "with its statistics",
    ),
    _Verb(
        "deviation-regions",
        "_deviation_regions",
        "deviation_regions_cmd",
        "Connected excess and missing-material regions beyond a deviation tolerance",
    ),
    _Verb(
        "exaggerate",
        "_exaggerate",
        "exaggerate_cmd",
        "Points displaced along their normals by a scaled deviation, for a figure",
    ),
    _Verb(
        "colorize",
        "_colorize",
        "colorize_cmd",
        "Bake a field into red/green/blue colours through a colormap (PLY, glTF)",
    ),
    _Verb(
        "feature-edges",
        "_feature_edges",
        "feature_edges_cmd",
        "Sharp, open and non-manifold edges of a surface as a line mesh",
    ),
    _Verb(
        "segment",
        "_segment",
        "segment_cmd",
        "Split a surface into smooth patches bounded by its feature edges",
    ),
    _Verb(
        "select",
        "_select",
        "select_cmd",
        (
            "Write a named point, cell or side region chosen by geometry, data, "
            "quality, component or a grown patch"
        ),
    ),
    _Verb(
        "fit-primitives",
        "_fit_primitives",
        "fit_primitives_cmd",
        "Fit a plane, sphere, cylinder, cone or torus to each patch",
    ),
    _Verb(
        "refine-patches",
        "_refine_patches",
        "refine_patches_cmd",
        (
            "Refine segment:patch until one surface explains each patch "
            "(variational shape approximation)"
        ),
    ),
    _Verb(
        "iga-refine",
        "_iga",
        "refine_cmd",
        "Refine a NURBS model exactly (knot insertion, degree elevation)",
        add_args="add_refine_args",
    ),
    _Verb(
        "iga-sample",
        "_iga",
        "sample_cmd",
        "Sample a NURBS model into Lagrange cells",
        add_args="add_sample_args",
    ),
    _Verb(
        "iga-evaluate",
        "_iga",
        "evaluate_cmd",
        "Evaluate a patch of a NURBS model at parametric coordinates (JSON)",
        add_args="add_evaluate_args",
    ),
    _Verb(
        "to-cad",
        "_to_cad",
        "to_cad_cmd",
        (
            "Fitted patches as a STEP solid with analytic, freeform, fillet and "
            "thread faces (analytic only without OpenCASCADE)"
        ),
    ),
    _Verb(
        "region-adjacency",
        "_region_adjacency",
        "region_adjacency_cmd",
        "Write facets shared by selected regions or cell blocks",
    ),
    _Verb(
        "find-interface",
        "_find_interface",
        "find_interface_cmd",
        "Find conforming or proximity-matched interface facets",
    ),
    _Verb(
        "contact-pairs",
        "_contact_pairs",
        "contact_pairs_cmd",
        "Project slave-region nodes onto master facets",
    ),
    _Verb(
        "split-interface",
        "_split_interface",
        "split_interface_cmd",
        "Split node fans along a Side region",
    ),
    _Verb(
        "hausdorff",
        "_hausdorff",
        "hausdorff_cmd",
        "Hausdorff distance between two surfaces (--max: nonzero exit above it)",
    ),
    _Verb(
        "coordinates",
        "_coordinates",
        "coordinates_cmd",
        "Positions, vector and tensor fields in cylindrical/spherical components",
    ),
    _Verb(
        "metric",
        "_metric",
        "metric_cmd",
        "Intersect, normalise, clamp and grade metric fields",
    ),
    _Verb(
        "untangle",
        "_untangle",
        "untangle_cmd",
        "Move nodes until no tetrahedron or hexahedron is inverted",
    ),
    _Verb(
        "morph",
        "_morph",
        "morph_cmd",
        "Move the boundary by a displacement, carried into the mesh",
    ),
    _Verb(
        "eigenmodes",
        "_eigenmodes",
        "eigenmodes_cmd",
        "Lowest eigenpairs of a surface's Laplace-Beltrami operator",
    ),
    _Verb(
        "register",
        "_register",
        "register_cmd",
        "Rigid registration of one mesh onto another (ICP)",
    ),
    _Verb(
        "register-landmarks",
        "_register",
        "register_landmarks_cmd",
        "The rigid motion laying picked points onto their targets",
        "add_landmark_args",
    ),
    _Verb(
        "auto-align",
        "_register",
        "auto_align_cmd",
        "Pose a scan on its largest flat face, or a cylinder on its axis",
        "add_auto_align_args",
    ),
    _Verb(
        "symmetry-plane",
        "_register",
        "symmetry_plane_cmd",
        "The plane a mesh is most nearly symmetric about",
        "add_symmetry_args",
    ),
    _Verb(
        "point-normals",
        "_point_normals",
        "point_normals_cmd",
        "Normals of a point cloud by local plane fits",
    ),
    _Verb(
        "reconstruct",
        "_reconstruct",
        "reconstruct_cmd",
        "A triangle surface through an oriented point cloud",
    ),
    _Verb(
        "hull",
        "_hull",
        "hull_cmd",
        "Convex hull of the points (3-D surface or 2-D loop)",
        add_args="add_hull_args",
    ),
    _Verb(
        "delaunay",
        "_hull",
        "delaunay_cmd",
        "Delaunay triangulation of the points in the x, y plane",
        add_args="add_delaunay_args",
    ),
    _Verb(
        "geodesic",
        "_geodesic",
        "geodesic_cmd",
        "Surface distance from source points (heat method)",
    ),
    _Verb(
        "parameterize",
        "_parameterize",
        "parameterize_cmd",
        "A UV atlas of a triangle surface (texture coordinates)",
    ),
    _Verb(
        "graph",
        "_graph",
        "graph_cmd",
        "Write a mesh's node or cell graph as a METIS or Scotch graph file",
    ),
    _Verb(
        "periodic",
        "_periodic",
        "periodic_cmd",
        "Match the nodes of two regions a transform maps onto each other",
    ),
    _Verb(
        "repair",
        "_repair",
        "repair_cmd",
        "Fix a surface's orientation, holes and pinched vertices",
    ),
    _Verb(
        "shrinkwrap",
        "_shrinkwrap",
        "shrinkwrap_cmd",
        "Project a mesh's points onto a target surface",
    ),
    _Verb(
        "sobolev-deform",
        "_sobolev",
        "sobolev_cmd",
        "Filter a displacement field over the mesh's P1 operators and apply it",
    ),
    _Verb(
        "extract-surface",
        "_extract_surface",
        "extract_surface_cmd",
        "Extract the boundary surface/edges",
        aliases=("surface",),
    ),
    _Verb(
        "reorder",
        "_reorder",
        "reorder_cmd",
        "Renumber nodes/elements (RCM / Morton / Hilbert)",
    ),
    _Verb(
        "diff",
        "_diff",
        "diff_cmd",
        "Compare two meshes (nonzero exit code if different)",
    ),
    _Verb(
        "merge",
        "_merge",
        "merge_cmd",
        "Merge two or more meshes into one (optional welding)",
    ),
    _Verb(
        "transform",
        "_transform",
        "transform_cmd",
        "Affine transform (translate / scale / rotate / matrix / units)",
    ),
    _Verb("clean", "_clean", "clean_cmd", "Weld / prune / de-dup a mesh in one pass"),
    _Verb(
        "crop",
        "_crop",
        "crop_cmd",
        "Subset a mesh by a bounding box, half-space or cell_data predicate",
    ),
    _Verb(
        "slice",
        "_slice",
        "slice_cmd",
        "Planar cross-section of a mesh (volume -> surface, surface -> lines)",
    ),
    _Verb(
        "fit-curves",
        "_fit_curves",
        "fit_curves_cmd",
        "Lines, arcs and circles from chained sections (a sketch)",
    ),
    _Verb(
        "sections",
        "_slice",
        "sections_cmd",
        "Chained line sections of a surface at parallel planes",
        "add_sections_args",
    ),
    _Verb(
        "isosurface",
        "_isosurface",
        "isosurface_cmd",
        "Level set of a scalar point_data field (volume -> surface, surface -> lines)",
    ),
    _Verb(
        "generate",
        "_generate",
        "generate_cmd",
        (
            "Build a mesh from parameters: rectangle, box, disk, spheres, "
            "cylinder, reference cells"
        ),
    ),
    _Verb(
        "voxelize",
        "_voxelize",
        "voxelize_cmd",
        "Regular hexahedron grid around a mesh (whole box / surface / interior)",
    ),
    _Verb(
        "octree",
        "_octree",
        "octree_cmd",
        "Octree / quadtree around a skin, fine only where the skin is",
    ),
    _Verb(
        "octree-mesh",
        "_octree_mesh",
        "octree_mesh_cmd",
        (
            "Volume mesh of a skin's inside on a balanced octree (castellated "
            "hexahedra, conforming tetrahedra or all-hexahedra)"
        ),
    ),
    _Verb(
        "mesh-cad-surface",
        "_mesh_cad",
        "mesh_cad_surface_cmd",
        (
            "Conforming triangle mesh of a CAD model's faces (CAD JSON, STEP, "
            "IGES) at a target edge size"
        ),
    ),
    _Verb(
        "project-to-cad",
        "_project_to_cad",
        "project_to_cad_cmd",
        "Move a mesh's points back onto the CAD faces and edges its cells name",
    ),
    _Verb(
        "boundary-layer",
        "_boundary_layer",
        "boundary_layer_cmd",
        "Prismatic boundary layers grown from a triangle and quad surface",
    ),
    _Verb(
        "extrude",
        "_extrude",
        "extrude_cmd",
        "Sweep a surface into a volume (or a curve into a surface) along a vector",
    ),
    _Verb(
        "revolve",
        "_revolve",
        "revolve_cmd",
        "Sweep a surface (or a curve) around an axis",
    ),
    _Verb(
        "mirror",
        "_mirror",
        "mirror_cmd",
        "Reflect a mesh across a plane, every cell right side out",
    ),
    _Verb(
        "tetrahedralize",
        "_tetrahedralize",
        "tetrahedralize_cmd",
        "Fill a closed triangle surface with tetrahedra, its triangles kept (Netgen)",
    ),
    _Verb(
        "remesh-quads",
        "_remesh_quads",
        "remesh_quads_cmd",
        "Field-aligned quads from a triangle surface (QuadriFlow, Instant Meshes)",
    ),
    _Verb(
        "cut",
        "_cut",
        "cut_cmd",
        "Cut a mesh by level sets: the pieces either side or the interface (CutCells)",
    ),
    _Verb(
        "cut-quadrature",
        "_cut",
        "cut_quadrature_cmd",
        "Quadrature rules on the part of a mesh that level sets select (CutCells)",
        add_args="add_quadrature_args",
    ),
    _Verb(
        "grid-sample",
        "_grid_transfer",
        "grid_sample_cmd",
        "Sample a mesh's point data onto a regular grid (for a CNN)",
        add_args="add_sample_args",
    ),
    _Verb(
        "grid-scatter",
        "_grid_transfer",
        "grid_scatter_cmd",
        "Write a grid's fields back onto a mesh (grid-sample's inverse)",
        add_args="add_scatter_args",
    ),
    _Verb(
        "grid-resample",
        "_grid_transfer",
        "grid_resample_cmd",
        "Resample a grid onto a finer or coarser lattice over the same box",
        add_args="add_resample_args",
    ),
    _Verb(
        "grid-spectrum",
        "_grid_transfer",
        "grid_spectrum_cmd",
        "Azimuthally averaged power spectrum of a field on a regular grid",
        add_args="add_spectrum_args",
    ),
    _Verb(
        "subsample",
        "_point_budget",
        "subsample_cmd",
        "Reduce a mesh to a point cloud under a token budget (farthest-point)",
    ),
    _Verb(
        "proximity-graph",
        "_proximity",
        "proximity_cmd",
        "Build a radius or k-nearest-neighbour graph over a mesh's points",
    ),
    _Verb(
        "predict",
        "_predict",
        "predict_cmd",
        "Run a trained PhysicsNeMo checkpoint on one mesh file",
    ),
    _Verb(
        "guard-fit",
        "_guard",
        "guard_fit_cmd",
        "Fit a geometry guardrail over a dataset manifest's split",
        add_args="add_fit_args",
    ),
    _Verb(
        "guard-check",
        "_guard",
        "guard_check_cmd",
        "Score a mesh against a fitted geometry guardrail",
        add_args="add_check_args",
    ),
    _Verb(
        "sdf",
        "_sdf",
        "sdf_cmd",
        "Signed distance field: generate a grid over a surface and fill it",
    ),
    _Verb(
        "split",
        "_split",
        "split_cmd",
        "Partition a mesh into multiple files (type / region / regions / component)",
    ),
    _Verb(
        "regions",
        "_regions",
        "regions_cmd",
        "List a mesh's named regions, or edit them (union, rename, retag, ...)",
    ),
    _Verb(
        "convert-cells",
        "_convert_cells",
        "convert_cells_cmd",
        "Convert the element representation (linearize / simplexify / elevate)",
    ),
    _Verb(
        "subdivide",
        "_subdivide",
        "subdivide_cmd",
        (
            "Polyhedrally refine: split every eligible 3D cell into one "
            "polyhedral child per face, connected to a new interior point"
        ),
    ),
    _Verb(
        "tessellate",
        "_tessellate",
        "tessellate_cmd",
        (
            "Isoparametric subdivision of a mesh's curved cells (quad9, quad8, "
            "triangle6, tetra10, hexahedron27) onto a refinement lattice, "
            "recording tessellate:* provenance"
        ),
    ),
    _Verb(
        "agglomerate",
        "_agglomerate",
        "agglomerate_cmd",
        (
            "Polyhedrally coarsen: merge groups of cells into single larger "
            "polyhedral cells via greedy seed-and-grow over the shared-face dual"
        ),
    ),
    _Verb(
        "tet-to-poly",
        "_tet_poly",
        "tet_to_poly_cmd",
        (
            "Convert a tetrahedral mesh into its polyhedral (median) dual: one "
            "polyhedron per vertex"
        ),
        add_args="add_tet_to_poly_args",
    ),
    _Verb(
        "poly-to-tet",
        "_tet_poly",
        "poly_to_tet_cmd",
        (
            "Decompose a polyhedral mesh into tetrahedra, keeping polyhedra that "
            "already are tetrahedra as one tetra"
        ),
        add_args="add_poly_to_tet_args",
    ),
    _Verb(
        "resample",
        "_resample",
        "resample_cmd",
        (
            "Resample a sequence onto new times (linear blend, nearest or "
            "previous), holding at most two steps in memory"
        ),
    ),
    _Verb(
        "refine",
        "_refine",
        "refine_cmd",
        (
            "Refine (subdivide cells into same-type children), every cell or a "
            "selected subset with a conforming closure"
        ),
    ),
    _Verb(
        "undo-green",
        "_undo_green",
        "undo_green_cmd",
        (
            "Restore a transitional (green) cell to its coarse parent, read "
            "verbatim from the coarse mesh"
        ),
    ),
    _Verb(
        "decimate",
        "_decimate",
        "decimate_cmd",
        "Reduce a surface mesh's face count by quadric edge collapse",
    ),
    _Verb(
        "decimate-volume",
        "_decimate_volume",
        "decimate_volume_cmd",
        "Reduce a tet mesh's cell count by quadric tet-edge collapse",
    ),
    _Verb(
        "adapt",
        "_mmg",
        "adapt_cmd",
        "Remesh to a size or a metric with MMG (a build with MESHIOPLUSPLUS_WITH_MMG)",
        add_args="add_adapt_args",
    ),
    _Verb(
        "discretize-levelset",
        "_mmg",
        "discretize_levelset_cmd",
        "Embed a level set conformally and remesh both sides with MMG",
        add_args="add_discretize_levelset_args",
    ),
    _Verb(
        "mesh-2d",
        "_mmg",
        "mesh_2d_cmd",
        "Triangulate the planar domain a line mesh bounds with MMG",
        add_args="add_mesh_2d_args",
    ),
    _Verb(
        "move",
        "_mmg",
        "move_cmd",
        "Move a mesh by a boundary displacement and remesh with MMG (needs ELAS)",
        add_args="add_move_args",
    ),
    _Verb(
        "remesh",
        "_remesh",
        "remesh_cmd",
        "Replace a surface's triangulation with a new, well-shaped one (ACVD)",
    ),
    _Verb(
        "remesh-volume",
        "_remesh_volume",
        "remesh_volume_cmd",
        "Retetrahedralize a volume (or closed surface) by isosurface stuffing",
    ),
    _Verb(
        "optimize-volume",
        "_optimize_volume",
        "optimize_volume_cmd",
        "ODT-remesh a tetrahedral mesh (relocate vertices + flip connectivity)",
    ),
    _Verb(
        "partition",
        "_partition",
        "partition_cmd",
        "Decompose into N balanced parts (SFC / KaHIP) for domain decomposition",
    ),
    _Verb(
        "pipeline",
        "_pipeline",
        "pipeline_cmd",
        "Run a settings.json operation chain (read -> ops -> write)",
    ),
    _Verb(
        "smooth",
        "_smooth",
        "smooth_cmd",
        "Relax node positions to improve element shape (Laplacian / Taubin)",
    ),
    _Verb(
        "interpolate",
        "_interpolate",
        "interpolate_cmd",
        (
            "Sample data arrays from a source mesh onto a target mesh (nearest / "
            "barycentric)"
        ),
    ),
    _Verb(
        "conservative-interpolate",
        "_conservative_interpolate",
        "conservative_interpolate_cmd",
        (
            "Mass-preserving cross-mesh field transfer (overlap-measure weighted,"
            " not pointwise)"
        ),
    ),
    _Verb(
        "stats",
        "_stats",
        "stats_cmd",
        "Print geometric statistics (bbox, area, volume, ...)",
    ),
    _Verb(
        "view",
        "_view",
        "view_cmd",
        "Open a mesh in an interactive viewer (desktop or browser)",
    ),
    _Verb(
        "screenshot",
        "_view",
        "screenshot_cmd",
        "Render a mesh to a PNG without opening a window",
        add_args="add_screenshot_args",
    ),
    _Verb(
        "snapshot",
        "_snapshot",
        "snapshot_cmd",
        (
            "Draw a mesh in this terminal (OUT = -) or to "
            ".png/.txt/.ansi/.html/.cast, with no display or GPU"
        ),
    ),
    _Verb(
        "tui",
        "_tui",
        "tui_cmd",
        "Orbit, zoom and pan a mesh inside this terminal, with no display or GPU",
    ),
    # Nested groups: `data <verb>` and `dataset <verb>` (the hand-editable
    # DatasetManifest JSON, doc/datasets.md). Their inner parsers each call
    # set_defaults(func=...), which overrides the outer default, so the
    # `args.func(args)` dispatch reaches the right handler unchanged.
    _Verb(
        "data",
        "_data",
        "data_cmd",
        "Inspect / rename / average / compute on mesh data arrays",
    ),
    _Verb(
        "dataset",
        "_dataset",
        "dataset_cmd",
        "Catalogue / split / tag / annotate a dataset manifest",
    ),
)

_COMPLETION_HELP = "Print a bash, zsh or fish completion script for this CLI"

_BY_NAME = {alias: v for v in _VERBS for alias in (v.name, *v.aliases)}

# ``_make_parser(only=_NO_VERB)`` registers every verb and fills in none.
_NO_VERB = object()


def _make_parser(only=None):
    """The top-level parser with every verb registered, so ``--help`` and the
    "invalid choice" error list them all; the arguments and handler are filled
    in for the verb ``only`` (a :class:`_Verb`), for none when it is
    ``_NO_VERB``, and for every verb, ``completion`` included, when it is
    ``None``."""
    parent_parser = argparse.ArgumentParser(
        description="Mesh input/output tools.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parent_parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=_get_version_text(),
        help="display version information",
    )

    subparsers = parent_parser.add_subparsers(title="subcommands", dest="command")
    subparsers.required = True

    for verb in _VERBS:
        parser = subparsers.add_parser(
            verb.name, help=verb.help, aliases=list(verb.aliases)
        )
        if only is None or verb is only:
            module = importlib.import_module(f".{verb.module}", __package__)
            getattr(module, verb.add_args)(parser)
            parser.set_defaults(func=getattr(module, verb.handler))

    parser = subparsers.add_parser("completion", help=_COMPLETION_HELP)
    if only is None:
        from . import _completion

        _completion.add_args(parser)
        parser.set_defaults(
            func=lambda args: _completion.completion_cmd(args, parent_parser)
        )

    return parent_parser


def build_parser():
    """The whole CLI's argparse parser; shell completions are generated from it
    (the ``completion`` verb), so they follow the verbs as they change."""
    return _make_parser()


def _verb_of(argv):
    """The name of the verb a command line runs, or ``None``. The top-level
    options take no value, so the verb is the first argument that is not an
    option."""
    for arg in argv:
        if not arg.startswith("-"):
            return arg
    return None


def main(argv=None):
    name = _verb_of(sys.argv[1:] if argv is None else argv)
    if name == "completion":
        # The completion script walks every verb's options.
        parser = build_parser()
    else:
        # No verb, ``--help``, ``--version`` or an unknown verb need only the
        # verbs' names and help lines, which every parser registers.
        parser = _make_parser(only=_BY_NAME.get(name, _NO_VERB))
    args = parser.parse_args(argv)
    # A handler that only reports returns None and errors leave by exception,
    # so None is success; a check verb returns its own non-zero codes.
    code = args.func(args)
    return 0 if code is None else code


def _get_version_text():
    python_version = f"{version_info.major}.{version_info.minor}.{version_info.micro}"
    return "\n".join(
        [
            f"meshio++ {__version__} [Python {python_version}]",
            "Copyright (c) 2015-2021 Nico Schlömer et al. (as meshio)",
            "Copyright (c) 2025 Vicente Mataix Ferrándiz",
            "Copyright (c) 2026 the meshio++ contributors",
        ]
    )
