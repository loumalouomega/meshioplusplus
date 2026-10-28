"""Shared pieces of the ``extrude``, ``revolve`` and ``mirror`` verbs."""

import json

from .._helpers import _writer_map, reader_map


def add_io_args(parser, out_help):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help=out_help)
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


def add_report_args(parser):
    parser.add_argument("--json", action="store_true", help="print the report as JSON")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="do not print the summary"
    )


def vec3(text, verb, flag):
    try:
        v = [float(x) for x in text.split(",")]
    except ValueError:
        v = []
    if len(v) != 3:
        raise ValueError(f"{verb}: --{flag} expects 'x,y,z'")
    return v


def print_sweep(args, report, revolve):
    if args.json:
        print(json.dumps(report))
        return
    if args.quiet:
        return
    text = (
        f"{report['num_cells']} cells ({report['num_hexahedra']} hexahedra, "
        f"{report['num_wedges']} wedges, {report['num_pyramids']} pyramids, "
        f"{report['num_tetrahedra']} tetrahedra, {report['num_quads']} quads, "
        f"{report['num_triangles']} triangles), {report['num_flipped']} profile "
        "cells turned"
    )
    if revolve:
        text += f", {report['num_on_axis']} points on the axis"
        if report["full"]:
            text += ", a full turn"
    print(f"{text}; least scaled Jacobian {report['min_scaled_jacobian']:.3f}")
    for w in report["warnings"]:
        print(f"warning: {w}")
