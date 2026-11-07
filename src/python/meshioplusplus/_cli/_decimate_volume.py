from .._decimate_volume import decimate_volume
from .._helpers import read, write
from ._common import add_decimate_args, decimate_summary

_HELP = {
    "in": "tet mesh file to be read from",
    "ratio": "fraction of the tets to KEEP, in (0, 1]",
    "count": "absolute tet count to stop at (within one collapse)",
    "max_error": "collapse only while the cheapest boundary-touching quadric "
    "error is at most this",
    "placement": (
        "where the surviving vertex goes: the quadric minimizer (midpoint "
        "when ill-conditioned or purely interior), the edge midpoint, or "
        "the cheaper endpoint (default: optimal)"
    ),
    "features": "allow boundary feature vertices (sharp corners/creases) to collapse",
    "angle": "degrees between boundary-triangle normals above which a vertex "
    "is a feature (default: 30)",
}

_BOUNDARY = (
    ("--preserve-boundary",),
    {
        "action": "store_true",
        "help": "pin every boundary vertex outright, reproducing `decimate`'s own "
        "default instead of letting boundary vertices participate",
    },
)


def add_args(parser):
    add_decimate_args(parser, "cells", _HELP, _BOUNDARY)


def decimate_volume_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = decimate_volume(
        mesh,
        ratio=args.ratio,
        target_cells=args.target_cells,
        max_error=args.max_error,
        placement=args.placement,
        preserve_boundary=args.preserve_boundary,
        preserve_features=args.preserve_features,
        feature_angle=args.feature_angle,
        return_report=True,
    )
    if not args.quiet:
        print(decimate_summary(out, report, "tets", "tets_removed"))
    write(args.outfile, out, file_format=args.output_format)
    return 0
