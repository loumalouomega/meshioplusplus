import argparse
from sys import version_info

from ..__about__ import __version__
from . import (
    _ascii,
    _binary,
    _clean,
    _compress,
    _convert,
    _convert_cells,
    _crop,
    _data,
    _decompress,
    _diff,
    _extract_surface,
    _info,
    _merge,
    _quality,
    _refine,
    _reorder,
    _split,
    _stats,
    _transform,
)


def main(argv=None):
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

    parser = subparsers.add_parser("convert", help="Convert mesh files", aliases=["c"])
    _convert.add_args(parser)
    parser.set_defaults(func=_convert.convert)

    parser = subparsers.add_parser("info", help="Print mesh info", aliases=["i"])
    _info.add_args(parser)
    parser.set_defaults(func=_info.info)

    parser = subparsers.add_parser("compress", help="Compress mesh file")
    _compress.add_args(parser)
    parser.set_defaults(func=_compress.compress)

    parser = subparsers.add_parser("decompress", help="Decompress mesh file")
    _decompress.add_args(parser)
    parser.set_defaults(func=_decompress.decompress)

    parser = subparsers.add_parser("ascii", help="Convert to ASCII", aliases=["a"])
    _ascii.add_args(parser)
    parser.set_defaults(func=_ascii.ascii)

    parser = subparsers.add_parser("binary", help="Convert to binary", aliases=["b"])
    _binary.add_args(parser)
    parser.set_defaults(func=_binary.binary)

    parser = subparsers.add_parser(
        "quality", help="Print mesh quality metrics", aliases=["q"]
    )
    _quality.add_args(parser)
    parser.set_defaults(func=_quality.quality)

    parser = subparsers.add_parser(
        "extract-surface",
        help="Extract the boundary surface/edges",
        aliases=["surface"],
    )
    _extract_surface.add_args(parser)
    parser.set_defaults(func=_extract_surface.extract_surface_cmd)

    parser = subparsers.add_parser(
        "reorder",
        help="Renumber nodes/elements (RCM / Morton / Hilbert)",
    )
    _reorder.add_args(parser)
    parser.set_defaults(func=_reorder.reorder_cmd)

    parser = subparsers.add_parser(
        "diff",
        help="Compare two meshes (nonzero exit code if different)",
    )
    _diff.add_args(parser)
    parser.set_defaults(func=_diff.diff_cmd)

    parser = subparsers.add_parser(
        "merge",
        help="Merge two or more meshes into one (optional welding)",
    )
    _merge.add_args(parser)
    parser.set_defaults(func=_merge.merge_cmd)

    parser = subparsers.add_parser(
        "transform",
        help="Affine transform (translate / scale / rotate / matrix / units)",
    )
    _transform.add_args(parser)
    parser.set_defaults(func=_transform.transform_cmd)

    parser = subparsers.add_parser(
        "clean",
        help="Weld / prune / de-dup a mesh in one pass",
    )
    _clean.add_args(parser)
    parser.set_defaults(func=_clean.clean_cmd)

    parser = subparsers.add_parser(
        "crop",
        help="Subset a mesh by a bounding box or half-space",
    )
    _crop.add_args(parser)
    parser.set_defaults(func=_crop.crop_cmd)

    parser = subparsers.add_parser(
        "split",
        help="Partition a mesh into multiple files (type / region / component)",
    )
    _split.add_args(parser)
    parser.set_defaults(func=_split.split_cmd)

    parser = subparsers.add_parser(
        "convert-cells",
        help="Convert the element representation (linearize / simplexify / elevate)",
    )
    _convert_cells.add_args(parser)
    parser.set_defaults(func=_convert_cells.convert_cells_cmd)

    parser = subparsers.add_parser(
        "refine",
        help="Uniformly refine (subdivide every cell into same-type children)",
    )
    _refine.add_args(parser)
    parser.set_defaults(func=_refine.refine_cmd)

    parser = subparsers.add_parser(
        "stats",
        help="Print geometric statistics (bbox, area, volume, ...)",
    )
    _stats.add_args(parser)
    parser.set_defaults(func=_stats.stats_cmd)

    # Nested group: `meshioplusplus data <verb>`. The inner parsers each call
    # set_defaults(func=...), which overrides the outer default, so the
    # `args.func(args)` dispatch below reaches the right handler unchanged.
    parser = subparsers.add_parser(
        "data",
        help="Inspect / rename / average / compute on mesh data arrays",
    )
    _data.add_args(parser)
    parser.set_defaults(func=_data.data_cmd)

    args = parent_parser.parse_args(argv)

    return args.func(args)


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
