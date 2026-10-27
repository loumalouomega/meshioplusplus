import json

from .._helpers import _writer_map, read, reader_map, write
from .._project_to_cad import project_to_cad


def add_args(parser):
    parser.add_argument(
        "infile",
        type=str,
        help="surface mesh whose cells name their CAD face (cad:face) and edge "
        "(cad:edge), as mesh-cad-surface writes them",
    )
    parser.add_argument(
        "model", type=str, help="the CAD model (CAD JSON, STEP or IGES)"
    )
    parser.add_argument("outfile", type=str, help="projected mesh to be written to")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input mesh format",
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
        "--face-array",
        type=str,
        default="cad:face",
        help="the cell data naming each cell's face (default cad:face)",
    )
    parser.add_argument(
        "--edge-array",
        type=str,
        default="cad:edge",
        help="the cell data naming each line's edge (default cad:edge)",
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=0.0,
        help="leave a point farther than this from where it would go (0: no limit)",
    )
    parser.add_argument(
        "--trim-samples",
        type=int,
        default=16,
        help="samples per knot span when testing a face's trims (default 16)",
    )
    parser.add_argument("--json", action="store_true", help="print the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="do not print the summary"
    )


def project_to_cad_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    model = read(args.model)
    out, report = project_to_cad(
        mesh,
        model,
        face_array=args.face_array,
        edge_array=args.edge_array,
        max_distance=args.max_distance,
        trim_samples=args.trim_samples,
        return_report=True,
    )
    if args.json:
        print(json.dumps(report))
    elif not args.quiet:
        print(
            f"projected {report['num_projected']} point(s) "
            f"({report['num_on_edges']} on edges, {report['num_on_curves']} onto "
            f"their curves, {report['num_corners']} at corners), largest move "
            f"{report['max_move']:.6g}; "
            f"{report['num_ambiguous']} ambiguous, {report['num_too_far']} too far, "
            f"{report['num_failed']} failed, {report['num_outside_trim']} outside "
            "their trims"
        )
    write(args.outfile, out, file_format=args.output_format)
