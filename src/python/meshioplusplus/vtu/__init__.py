from .. import _core
from .._exceptions import WriteError
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._vtu import read as _py_read
from ._vtu import write as _py_write


def read(filename, points_only=False, arrays=None):
    """Read a VTU file.

    Uses the C++ core for ascii and inline binary (uncompressed or zlib) files,
    falling back to the reference Python reader for anything it doesn't handle
    (lzma, appended/raw binary, multi-piece).
    """
    # points_only/arrays reach the C++ reader, which skips the unwanted
    # <DataArray>/section bodies outright. The Python fallback below has no
    # selective support, so _helpers.read trims its result instead -- same
    # answer, just without the saving.
    if not is_buffer(filename, "r"):
        try:
            return _core.vtu_read(str(filename), points_only=points_only, arrays=arrays)
        except Exception as exc:
            if not core_declined(exc, "vtu", "read", filename):
                raise
    return _py_read(filename)


def write(
    filename, mesh, binary=True, compression="zlib", header_type=None, appended=False
):
    """Write a VTU file.

    Uses the C++ core for the cases it supports (ASCII and binary, the latter
    uncompressed or zlib-compressed, polyhedra included, written to a real file
    path) and otherwise falls back to the reference Python writer. The
    fallback also catches any limitation hit by the C++ path, so behaviour is
    identical to the pure-Python implementation.

    ``appended=True`` (since v16.20.0) writes every piece array into one raw
    ``<AppendedData>`` section -- binary, with no base64 -- compressed with
    ``compression`` as usual. It needs ``binary=True``.
    """
    if appended and not binary:
        raise WriteError(
            "VTU: appended data is binary; appended=True needs binary=True"
        )
    # zstd/lz4 reach the C++ writer when this build has them; otherwise the
    # Python reference path handles it (needing the `codecs` extra), so a
    # request never silently degrades to a different codec.
    _CPP_CODECS = {None: "none", "zlib": "zlib", "lz4": "lz4", "zstd": "zstd"}
    cpp_compression_ok = compression in _CPP_CODECS
    if compression == "lz4" and not getattr(_core, "__has_lz4__", False):
        cpp_compression_ok = False
    if compression == "zstd" and not getattr(_core, "__has_zstd__", False):
        cpp_compression_ok = False
    if header_type is None and cpp_compression_ok and not is_buffer(filename, "w"):
        try:
            if appended:
                _core.vtu_write_appended(str(filename), mesh, _CPP_CODECS[compression])
            else:
                _core.vtu_write_codec(
                    str(filename), mesh, binary, _CPP_CODECS[compression]
                )
            return
        except Exception as exc:
            if not core_declined(exc, "vtu", "write", filename):
                raise
    return _py_write(
        filename,
        mesh,
        binary=binary,
        compression=compression,
        header_type=header_type,
        appended=appended,
    )


register_format("vtu", [".vtu"], read, {"vtu": write})

__all__ = ["read", "write"]
