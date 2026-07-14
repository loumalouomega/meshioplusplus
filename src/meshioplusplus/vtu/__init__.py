from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._vtu import read as _py_read
from ._vtu import write as _py_write


def _has_polyhedron(mesh):
    return any(c.type.startswith("polyhedron") for c in mesh.cells)


def read(filename):
    """Read a VTU file.

    Uses the C++ core for ascii and inline binary (uncompressed or zlib) files,
    falling back to the reference Python reader for anything it doesn't handle
    (lzma, appended/raw binary, polyhedron, multi-piece).
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.vtu_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    """Write a VTU file.

    Uses the C++ core for the cases it supports (ASCII and binary, the latter
    uncompressed or zlib-compressed, for non-polyhedron meshes written to a real
    file path) and otherwise falls back to the reference Python writer. The
    fallback also catches any limitation hit by the C++ path, so behaviour is
    identical to the pure-Python implementation.
    """
    cpp_compression_ok = compression is None or compression == "zlib"
    if (
        header_type is None
        and cpp_compression_ok
        and not is_buffer(filename, "w")
        and not _has_polyhedron(mesh)
    ):
        try:
            _core.vtu_write(str(filename), mesh, binary, compression == "zlib")
            return
        except Exception:
            pass
    return _py_write(
        filename, mesh, binary=binary, compression=compression, header_type=header_type
    )


register_format("vtu", [".vtu"], read, {"vtu": write})

__all__ = ["read", "write"]
