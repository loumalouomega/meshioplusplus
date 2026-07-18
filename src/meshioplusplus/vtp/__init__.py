from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._vtp import read as _py_read
from ._vtp import write as _py_write


def read(filename):
    """Read a VTP (VTK XML PolyData) file.

    Uses the C++ core for ascii and inline binary (uncompressed or zlib) files,
    falling back to the reference Python reader for anything it doesn't handle
    (lzma, appended data, triangle strips, multi-piece).
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.vtp_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    """Write a VTP (VTK XML PolyData) file.

    Uses the C++ core for the cases it supports (ASCII and binary, the latter
    uncompressed or zlib-compressed, written to a real file path) and otherwise
    falls back to the reference Python writer. The fallback also catches any
    limitation hit by the C++ path (e.g. cell types PolyData cannot hold raise
    WriteError in both implementations).
    """
    cpp_compression_ok = compression is None or compression == "zlib"
    if header_type is None and cpp_compression_ok and not is_buffer(filename, "w"):
        try:
            _core.vtp_write(str(filename), mesh, binary, compression == "zlib")
            return
        except Exception:
            pass
    return _py_write(
        filename, mesh, binary=binary, compression=compression, header_type=header_type
    )


register_format("vtp", [".vtp"], read, {"vtp": write})

__all__ = ["read", "write"]
