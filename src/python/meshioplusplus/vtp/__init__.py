from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import vtk_xml_shims
from ._vtp import read as _py_read  # noqa: F401
from ._vtp import write as _py_write  # noqa: F401

read, write = vtk_xml_shims(
    __name__,
    "vtp",
    """Read a VTP (VTK XML PolyData) file.

    Uses the C++ core for ascii and inline binary (uncompressed or zlib) files,
    falling back to the reference Python reader for anything it doesn't handle
    (lzma, appended data, triangle strips, multi-piece).
    """,
    """Write a VTP (VTK XML PolyData) file.

    Uses the C++ core for the cases it supports (ASCII and binary, the latter
    uncompressed or zlib-compressed, written to a real file path) and otherwise
    falls back to the reference Python writer. The fallback also catches any
    limitation hit by the C++ path (e.g. cell types PolyData cannot hold raise
    WriteError in both implementations).
    """,
)


register_format("vtp", [".vtp"], read, {"vtp": write})

__all__ = ["read", "write"]
