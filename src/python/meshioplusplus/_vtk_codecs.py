"""The block codecs the C++ VTK XML writers accept.

``vti``, ``vtm``, ``vtp``, ``vtr``, ``vts``, ``vtu``, ``pvd`` and the parallel
formats (:mod:`meshioplusplus._parallel_xml`) all route a write to the core only
when it can produce the requested codec; anything else (a custom
``header_type``, a codec this build lacks) goes to the Python reference, so a
request never silently degrades to a different codec.
"""

from . import _core

# Python ``compression`` argument -> the core's codec token.
CPP_CODECS = {
    None: "none",
    "zlib": "zlib",
    "lz4": "lz4",
    "zstd": "zstd",
    "lzma": "lzma",
}


def cpp_codec_ok(compression, header_type):
    """Whether the C++ writer can encode ``compression`` with ``header_type``."""
    if header_type is not None or compression not in CPP_CODECS:
        return False
    if compression == "lz4" and not getattr(_core, "__has_lz4__", False):
        return False
    if compression == "zstd" and not getattr(_core, "__has_zstd__", False):
        return False
    if compression == "lzma" and not getattr(_core, "__has_lzma__", False):
        return False
    return True
