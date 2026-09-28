from .._hausdorff import hausdorff_distance
from .._helpers import read, reader_map
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile_a", type=str, help="first mesh file")
    parser.add_argument("infile_b", type=str, help="second mesh file")
    for side in ("a", "b"):
        parser.add_argument(
            f"--input-format-{side}",
            type=str,
            choices=sorted(list(reader_map.keys())),
            help=f"input format of the {'first' if side == 'a' else 'second'} file",
            default=None,
        )
        parser.add_argument(
            f"--region-{side}",
            type=str,
            default="",
            help="restrict to this named cell region of surface cells",
        )
    parser.add_argument(
        "--face-samples",
        type=int,
        default=0,
        metavar="S",
        help="also sample the centroids of the S*S sub-triangles of every triangle "
        "(default: 0, vertices only -- a lower bound)",
    )
    parser.add_argument(
        "--max",
        type=float,
        default=None,
        metavar="D",
        help="exit with status 1 when the distance exceeds D",
    )
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")


def hausdorff_cmd(args):
    a = read(args.infile_a, file_format=args.input_format_a)
    b = read(args.infile_b, file_format=args.input_format_b)
    r = hausdorff_distance(
        a,
        b,
        face_samples=args.face_samples,
        region_a=args.region_a,
        region_b=args.region_b,
    )
    failed = args.max is not None and r["distance"] > args.max
    if args.json:
        out = dict(r)
        if args.max is not None:
            out["max"] = args.max
            out["passed"] = not failed
        emit_json(out)
    else:
        print(f"Hausdorff distance: {r['distance']:.9g}")
        print(
            f"  A -> B: max {r['a_to_b']:.9g}, mean {r['mean_a_to_b']:.9g}, "
            f"rms {r['rms_a_to_b']:.9g} ({r['num_samples_a']} samples)"
        )
        print(
            f"  B -> A: max {r['b_to_a']:.9g}, mean {r['mean_b_to_a']:.9g}, "
            f"rms {r['rms_b_to_a']:.9g} ({r['num_samples_b']} samples)"
        )
        if args.max is not None:
            print(f"  {'FAIL' if failed else 'pass'}: limit {args.max:g}")
    return 1 if failed else 0
