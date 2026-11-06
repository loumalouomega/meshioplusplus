from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from .._vtk_codecs import CPP_CODECS, cpp_codec_ok
from ._vts import read as _py_read
from ._vts import write as _py_write


def read(filename, points_only=False, arrays=None):
    """Read a VTK XML StructuredGrid file."""
    if not is_buffer(filename, "r"):
        try:
            return _core.vts_read(str(filename), points_only=points_only, arrays=arrays)
        except Exception as exc:
            if not core_declined(exc, "vts", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    """Write a mesh as VTK XML StructuredGrid. The mesh must be a dense lattice."""
    if cpp_codec_ok(compression, header_type) and not is_buffer(filename, "w"):
        try:
            _core.vts_write_codec(str(filename), mesh, binary, CPP_CODECS[compression])
            return
        except Exception as exc:
            if not core_declined(exc, "vts", "write", filename):
                raise
    return _py_write(
        filename, mesh, binary=binary, compression=compression, header_type=header_type
    )


register_format("vts", [".vts"], read, {"vts": write})

__all__ = ["read", "write"]
