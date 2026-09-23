from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._nastran import read as _py_read
from ._nastran import write as _py_write


def read(filename):
    """Read a Nastran or OptiStruct bulk-data file.

    Uses the C++ core; the reference Python reader is the fallback (and the
    only path for file-like buffers). HyperMesh components and OptiStruct
    SET cards become ``mesh.regions``.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.nastran_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "nastran", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh, point_format="fixed-large", cell_format="fixed-small"):
    """Write a Nastran bulk-data file.

    Uses the C++ core for the default fixed-large/fixed-small layout;
    otherwise falls back to the Python writer.
    """
    if (
        point_format == "fixed-large"
        and cell_format == "fixed-small"
        and not is_buffer(filename, "w")
    ):
        try:
            _core.nastran_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "nastran", "write", filename):
                raise
    return _py_write(filename, mesh, point_format=point_format, cell_format=cell_format)


register_format("nastran", [".bdf", ".fem", ".nas"], read, {"nastran": write})

__all__ = ["read", "write"]
