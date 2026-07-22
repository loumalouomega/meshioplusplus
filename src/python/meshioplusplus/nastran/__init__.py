from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._nastran import read as _py_read
from ._nastran import write as _py_write


def read(filename):
    """Read a Nastran bulk-data file.

    Uses the C++ core for files written by this library (recognised via a
    sentinel comment); all other files use the reference Python reader.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.nastran_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, point_format="fixed-large", cell_format="fixed-small"):
    """Write a Nastran bulk-data file.

    Uses the C++ core for the default fixed-large/fixed-small layout on meshes
    without nastran:ref data; otherwise falls back to the Python writer.
    """
    if (
        point_format == "fixed-large"
        and cell_format == "fixed-small"
        and not is_buffer(filename, "w")
        and "nastran:ref" not in mesh.point_data
        and "nastran:ref" not in mesh.cell_data
    ):
        try:
            _core.nastran_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, point_format=point_format, cell_format=cell_format)


register_format("nastran", [".bdf", ".fem", ".nas"], read, {"nastran": write})

__all__ = ["read", "write"]
