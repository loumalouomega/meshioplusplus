from .._decimate import decimate
from .._helpers import read, write
from ._common import add_decimate_args, decimate_summary

_HELP = {
    "in": "surface mesh file to be read from",
    "ratio": "fraction of the (triangulated) faces to KEEP, in (0, 1]",
    "count": "absolute face count to stop at (within one collapse)",
    "max_error": "collapse only while the cheapest quadric error is at most this",
    "placement": (
        "where the surviving vertex goes: the quadric minimizer (midpoint "
        "when ill-conditioned), the edge midpoint, or the cheaper endpoint "
        "(default: optimal)"
    ),
    "features": "allow feature vertices (sharp corners/creases) to collapse",
    "angle": "degrees between face normals above which a vertex is a feature "
    "(default: 30)",
}

_BOUNDARY = (
    ("--no-preserve-boundary",),
    {
        "dest": "preserve_boundary",
        "action": "store_false",
        "help": "allow boundary vertices to collapse (the outline may change)",
    },
)


def add_args(parser):
    add_decimate_args(parser, "faces", _HELP, _BOUNDARY)


def decimate_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = decimate(
        mesh,
        ratio=args.ratio,
        target_faces=args.target_faces,
        max_error=args.max_error,
        placement=args.placement,
        preserve_boundary=args.preserve_boundary,
        preserve_features=args.preserve_features,
        feature_angle=args.feature_angle,
        return_report=True,
    )
    if not args.quiet:
        print(decimate_summary(out, report, "faces", "faces_removed"))
    write(args.outfile, out, file_format=args.output_format)
    return 0
