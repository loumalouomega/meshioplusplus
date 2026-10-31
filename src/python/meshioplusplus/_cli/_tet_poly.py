from .._helpers import _writer_map, read, reader_map, write
from .._tet_poly import poly_to_tet, tet_to_poly
from ._json import emit_json


def _add_io_args(parser, out_help):
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


def add_tet_to_poly_args(parser):
    _add_io_args(parser, "polyhedral dual mesh to be written to")
    parser.add_argument(
        "--feature-angle",
        type=float,
        default=1.0,
        metavar="DEG",
        help="boundary triangles whose normals differ by more start a new dual "
        "boundary face (default: 1)",
    )
    parser.add_argument(
        "--point-data-to-cell-data",
        action="store_true",
        help="also copy each vertex's point_data onto its dual cells as cell_data",
    )
    parser.add_argument(
        "--record-vertex-ids",
        action="store_true",
        help="attach a dual:vertex cell_data array with each cell's source vertex",
    )
    parser.add_argument("--json", action="store_true", help="emit a summary as JSON")


def add_poly_to_tet_args(parser):
    _add_io_args(parser, "tetrahedral mesh to be written to")
    parser.add_argument(
        "--record-parent-ids",
        action="store_true",
        help="attach a convert:parent_cell cell_data array with each tet's input cell",
    )
    parser.add_argument("--json", action="store_true", help="emit a summary as JSON")


def tet_to_poly_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = tet_to_poly(
        mesh,
        feature_angle=args.feature_angle,
        point_data_to_cell_data=args.point_data_to_cell_data,
        record_vertex_ids=args.record_vertex_ids,
        return_report=True,
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(
            {
                "cells_in": sum(len(c) for c in mesh.cells),
                "cells_out": sum(len(c) for c in out.cells),
                "num_new_points": report["num_new_points"],
                "num_boundary_faces": report["num_boundary_faces"],
                "num_split_vertices": report["num_split_vertices"],
            }
        )
    return 0


def poly_to_tet_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    out, report = poly_to_tet(
        mesh, record_parent_ids=args.record_parent_ids, return_report=True
    )
    write(args.outfile, out, file_format=args.output_format)
    if args.json:
        emit_json(
            {
                "cells_in": sum(len(c) for c in mesh.cells),
                "cells_out": sum(len(c) for c in out.cells),
                "num_exact_tets": report["num_exact_tets"],
                "num_fanned_cells": report["num_fanned_cells"],
            }
        )
    return 0
