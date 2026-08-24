from .._helpers import _writer_map, read, reader_map, write
from .._remesh import remesh


def add_args(parser):
    parser.add_argument("infile", type=str, help="surface mesh file to be read from")
    parser.add_argument("outfile", type=str, help="remeshed mesh to be written to")
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
        "--num-clusters",
        type=int,
        required=True,
        help="number of clusters, i.e. output vertices aimed for "
        "(>= 4, <= the subdivided input's vertex count)",
    )
    parser.add_argument(
        "--subdivide",
        type=int,
        default=None,
        help="uniform refine passes before clustering; unset picks the "
        "smallest count reaching --subsample-ratio items/cluster "
        "(capped at --max-subdivide); 0 disables subdivision",
    )
    parser.add_argument(
        "--subsample-ratio",
        type=float,
        default=10.0,
        help="items per cluster targeted by automatic subdivision (default: 10)",
    )
    parser.add_argument(
        "--max-subdivide",
        type=int,
        default=4,
        help="ceiling on automatic subdivision, ignored when --subdivide is "
        "set explicitly (default: 4)",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=100,
        help="maximum energy-minimisation sweeps per pass (default: 100)",
    )
    parser.add_argument(
        "--repair-passes",
        type=int,
        default=10,
        help="'split disconnected clusters, minimise again' passes; 0 skips "
        "repair (default: 10)",
    )
    parser.add_argument(
        "--metric",
        choices=["isotropic", "quadric", "anisotropic"],
        default="isotropic",
        help="isotropic (area-weighted centroidal distance, default), "
        "quadric (Garland-Heckbert error, preserves sharp edges/corners), "
        "or anisotropic (clusters shaped by a local curvature tensor -- "
        "see --max-anisotropy)",
    )
    parser.add_argument(
        "--max-anisotropy",
        type=float,
        default=4.0,
        help="under --metric anisotropic, the maximum ratio between the "
        "two in-plane target edge lengths (long axis / short axis) a "
        "per-vertex curvature tensor may request; 1.0 recovers the "
        "isotropic shape exactly (default: 4.0, a measured value -- see "
        "kRemeshDefaultMaxAnisotropy in remesh.hpp); an error to set away "
        "from the default under any other metric",
    )
    parser.add_argument(
        "--gradation",
        type=float,
        default=0.0,
        dest="gradation",
        help="curvature-gradation exponent gamma in the item weight "
        "area * kappa**gamma; 0 (default) disables gradation and reproduces "
        "plain area weighting",
    )
    parser.add_argument(
        "--no-preserve-boundary",
        dest="preserve_boundary",
        action="store_false",
        help="do not pin the input's open boundary or emit a boundary line "
        "block (a no-op on a closed mesh either way)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="do not print the run summary",
    )
    parser.set_defaults(preserve_boundary=True)


def remesh_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = remesh(
        mesh,
        args.num_clusters,
        subdivide=args.subdivide,
        subsample_ratio=args.subsample_ratio,
        max_subdivide=args.max_subdivide,
        max_iterations=args.iterations,
        max_repair_passes=args.repair_passes,
        metric=args.metric,
        gradation=args.gradation,
        preserve_boundary=args.preserve_boundary,
        max_anisotropy=args.max_anisotropy,
        return_report=True,
    )
    if not args.quiet:
        print(
            f"remeshed to {report['num_clusters']} clusters "
            f"({report['num_iterations']} iterations, "
            f"{report['subdivide_applied']} subdivide passes, "
            f"{report['num_isolated_clusters']} isolated clusters, "
            f"{report['num_non_manifold_vertices']} non-manifold vertices)"
        )
    write(args.outfile, out, file_format=args.output_format)
    return 0
