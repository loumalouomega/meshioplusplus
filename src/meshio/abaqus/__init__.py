from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._abaqus import read as _py_read
from ._abaqus import write as _py_write


def read(filename):
    """Read an Abaqus .inp file.

    Uses the C++ core for files limited to *NODE and *ELEMENT; falls back to the
    Python reader for *NSET / *ELSET / *INCLUDE and anything else.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.abaqus_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, float_fmt=".16e", translate_cell_names=True):
    """Write an Abaqus .inp file.

    Uses the C++ core for meshes without point_sets/cell_sets; otherwise falls
    back to the Python writer.
    """
    if (
        float_fmt == ".16e"
        and translate_cell_names
        and not is_buffer(filename, "w")
        and not mesh.point_sets
        and not mesh.cell_sets
    ):
        try:
            _core.abaqus_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(
        filename, mesh, float_fmt=float_fmt, translate_cell_names=translate_cell_names
    )


register_format("abaqus", [".inp"], read, {"abaqus": write})

__all__ = ["read", "write"]
