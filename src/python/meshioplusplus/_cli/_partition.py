import os

from .._helpers import _writer_map, read, reader_map, write
from .._partition import partition, partition_labels


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument(
        "outpattern",
        type=str,
        help=(
            "output path pattern containing '{part}' (e.g. 'out_{part}.vtu'), "
            "expanded once per piece; or a '.pvtu'/'.pvtp' path, written as one "
            "index over every piece (halo layers included); with --labels-only, "
            "a single plain path"
        ),
    )
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
        "--nparts",
        "-n",
        type=int,
        required=True,
        help="number of parts to decompose into (>= 1)",
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=["auto", "sfc", "kahip"],
        default="auto",
        help=(
            "partitioning backend: 'sfc' (Hilbert curve cut, always available), "
            "'kahip' (KaHIP kaffpa on the dual graph; needs a KaHIP-enabled "
            "build and fails by name otherwise), or 'auto' (kahip when "
            "available, else sfc; default)"
        ),
    )
    parser.add_argument(
        "--imbalance",
        type=float,
        default=0.03,
        help="kahip only: allowed imbalance fraction in (0, 1) (default: 0.03)",
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["fast", "eco", "strong"],
        default="eco",
        help="kahip only: preconfiguration (default: eco)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="kahip only: random seed (deterministic per seed; default: 0)",
    )
    parser.add_argument(
        "--weights",
        type=str,
        default=None,
        help="name of a scalar cell_data array of per-cell weights",
    )
    parser.add_argument(
        "--ghost-layers",
        type=int,
        default=0,
        help="grow each piece by N shared-node BFS layers of other parts' "
        "cells (a halo), tagged partition:ghost (0 = owned)",
    )
    parser.add_argument(
        "--record-ids",
        action="store_true",
        help=(
            "attach partition:original_point_id / partition:original_cell_id "
            "arrays of the original input indices to every piece"
        ),
    )
    parser.add_argument(
        "--labels-only",
        action="store_true",
        help=(
            "write the input mesh once with the assignment attached as the "
            "Int64 'partition:part' cell_data instead of writing pieces"
        ),
    )


def _parallel_index_kind(outpattern, output_format):
    """``"pvtu"`` / ``"pvtp"`` when the output is a parallel index, else ``None``."""
    if output_format is not None:
        return output_format if output_format in ("pvtu", "pvtp") else None
    return {".pvtu": "pvtu", ".pvtp": "pvtp"}.get(
        os.path.splitext(outpattern)[1].lower()
    )


def partition_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)

    if args.labels_only:
        labels = partition_labels(
            mesh,
            args.nparts,
            method=args.method,
            imbalance=args.imbalance,
            mode=args.mode,
            seed=args.seed,
            weights=args.weights,
        )
        mesh.cell_data["partition:part"] = labels
        write(args.outpattern, mesh, file_format=args.output_format)
        return 0

    index_kind = _parallel_index_kind(args.outpattern, args.output_format)
    if "{part}" not in args.outpattern and index_kind is None:
        raise ValueError(
            "partition: the output pattern must contain '{part}' "
            "(e.g. 'out_{part}.vtu'), or be a .pvtu/.pvtp index (one file per "
            "part plus an index), or pass --labels-only"
        )
    pieces = partition(
        mesh,
        args.nparts,
        method=args.method,
        imbalance=args.imbalance,
        mode=args.mode,
        seed=args.seed,
        record_ids=args.record_ids,
        ghost_layers=args.ghost_layers,
        weights=args.weights,
    )
    if "{part}" not in args.outpattern:
        # One index over every piece, so the halo layers survive as vtkGhostType.
        from .. import pvtp, pvtu

        (pvtu if index_kind == "pvtu" else pvtp).write_pieces(args.outpattern, pieces)
        return 0
    for part_id, piece in enumerate(pieces):
        write(
            args.outpattern.format(part=part_id),
            piece,
            file_format=args.output_format,
        )
    return 0
