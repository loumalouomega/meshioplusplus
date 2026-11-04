"""The ``geodesic`` verb: the surface distance from source points by the heat
method (``doc/geodesic.md``)."""

from .._geodesic import geodesic_distance
from .._helpers import _reader_names, read, write
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="triangle surface mesh")
    parser.add_argument("outfile", type=str, help="output mesh with geodesic:distance")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=_reader_names(),
        help="input file format",
        default=None,
    )
    parser.add_argument("--output-format", "-o", type=str, default=None)
    parser.add_argument(
        "--source",
        type=str,
        default="",
        metavar="ID[,ID...]",
        help="source point ids, 0-based, comma-separated",
    )
    parser.add_argument(
        "--source-region",
        type=str,
        default="",
        help="a region whose points are sources",
    )
    parser.add_argument(
        "--time-factor",
        type=float,
        default=1.0,
        help="heat time = factor * h^2, h the mean edge length (default 1)",
    )
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--max-iterations", type=int, default=0)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def geodesic_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    try:
        ids = [int(v) for v in args.source.split(",") if v.strip()]
    except ValueError:
        raise SystemExit(
            "geodesic: --source expects comma-separated integers"
        ) from None
    out, report = geodesic_distance(
        mesh,
        ids,
        source_region=args.source_region or None,
        time_factor=args.time_factor,
        tolerance=args.tolerance,
        max_iterations=args.max_iterations,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json({"output": args.outfile, **report})
        return 0
    print(
        f"geodesic distance from {report['num_sources']} source(s): max "
        f"{report['max_distance']:.6g}, {report['num_unreached']} unreached point(s)"
    )
    return 0
