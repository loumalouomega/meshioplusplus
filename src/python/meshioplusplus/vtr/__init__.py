from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from .._vtk_codecs import CPP_CODECS, cpp_codec_ok
from ._vtr import read as _py_read
from ._vtr import write as _py_write


def read(filename, points_only=False, arrays=None):
    """Read a VTK XML RectilinearGrid file. Fully general: coordinates need
    not be uniformly spaced (a graded mesh reads correctly)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.vtr_read(str(filename), points_only=points_only, arrays=arrays)
        except Exception as exc:
            if not core_declined(exc, "vtr", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    """Write a mesh as VTK XML RectilinearGrid.

    Unlike the reader, the mesh must be a UNIFORM dense lattice today: a
    genuinely graded mesh cannot be written as ``.vtr`` yet.
    """
    if cpp_codec_ok(compression, header_type) and not is_buffer(filename, "w"):
        try:
            _core.vtr_write_codec(str(filename), mesh, binary, CPP_CODECS[compression])
            return
        except Exception as exc:
            if not core_declined(exc, "vtr", "write", filename):
                raise
    return _py_write(
        filename, mesh, binary=binary, compression=compression, header_type=header_type
    )


register_format("vtr", [".vtr"], read, {"vtr": write})

__all__ = ["read", "write"]
