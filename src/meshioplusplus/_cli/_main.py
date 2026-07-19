import argparse
from sys import version_info

from ..__about__ import __version__
from . import (
    _ascii,
    _binary,
    _compress,
    _convert,
    _decompress,
    _diff,
    _extract_surface,
    _info,
    _quality,
    _reorder,
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
